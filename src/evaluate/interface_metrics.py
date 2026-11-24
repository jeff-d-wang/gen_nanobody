#!/usr/bin/env python3
"""B3: interface and geometry metrics on all three structure sets, per candidate.

For each of the 84 RFantibody candidates, the same measurements on:

  generated  A's ProteinMPNN output, the designed sequence on the RFdiffusion
             backbone in the docked pose it was generated in
  rf2        A's RF2 prediction of that candidate (a PIPELINE predictor, so it
             can never be the independent final scorer; see non-negotiable #1)
  boltz      B's Boltz prediction from sequence alone, no template, no pose

Three contact counts are kept separate, as B3 requires, each with its denominator
stated in the output: hotspot contacts over the 4 hotspots A conditioned on,
full-epitope contacts over C's 5 A epitope, and the recovered fraction of each.

TWO THINGS MAKE THIS LESS OBVIOUS THAN IT LOOKS.

1. The generated complex has NO SIDE CHAINS. A's ProteinMPNN output carries
   exactly N, CA, C, O; RF2 and Boltz carry full side chains. A heavy-atom 5 A
   interface is therefore not comparable across the three sets: the generated
   complex would look artificially sparse. So every contact count is computed
   twice, once on backbone atoms only (comparable across all three) and once on
   all heavy atoms (comparable between rf2 and boltz, and matching C's epitope
   definition and fold_metrics.py). `has_sidechains` flags which rows support
   the heavy-atom reading. Compare like with like.

2. A's files continue antigen numbering after the nanobody. Chain T runs
   117..183 in a 116 + 67 complex, while Boltz numbers chain B 1..67. The
   hotspot lists and C's epitope sets are in 1..N T numbering, so A's residues
   are shifted back by `t_offset` before any set comparison. The offset is read
   from the structure, then checked against the antigen length from the manifest
   and against contiguity, because a silently wrong offset would land every
   lookup on the wrong residue and still produce plausible numbers.

Buried area uses Shrake-Rupley from Biopython, as
(SASA_nanobody + SASA_antigen - SASA_complex) / 2, hydrogens stripped so RF2's
explicit hydrogens do not inflate it. It is only meaningful with side chains, so
it is left blank for the generated complex rather than reported as zero.

Rosetta InterfaceAnalyzer is deliberately absent. PLAN.md makes it droppable, and
contacts, clashes and buried area cover the same ground without an academic
licence or a build.

  python3 src/evaluate/interface_metrics.py --demo
  python3 src/evaluate/interface_metrics.py
"""
import argparse
import copy
import csv
import glob
import json
import os
import warnings

import numpy as np

import fold_metrics as fm

BACKBONE = ("N", "CA", "C", "O")
CLASH = 2.0          # heavy-atom separation below this is a steric clash
NB_AG = {"generated": ("H", "T"), "rf2": ("H", "T"), "boltz": ("A", "B")}


def candidate_paths(run_dir, pdb_id, arm, idx):
    """(generated, rf2) structure paths for one candidate, or (None, None)."""
    base = os.path.join(run_dir, "%s_%s" % (pdb_id, arm))
    gen = os.path.join(base, "proteinmpnn_attempt_1", "outputs",
                       "nb_%d_dldesign_0.pdb" % idx)
    rf2 = os.path.join(base, "rf2_attempt_1", "outputs",
                       "nb_%d_dldesign_0_best.pdb" % idx)
    return (gen if os.path.exists(gen) else None,
            rf2 if os.path.exists(rf2) else None)


def boltz_path(boltz_dir, tag, job):
    hits = sorted(glob.glob(os.path.join(boltz_dir, tag, job, "**", "*model_0.pdb"),
                            recursive=True))
    return hits[0] if hits else None


def antigen_offset(struct, ag_chain, n_expected):
    """How much to subtract from this file's antigen ids to reach 1..N numbering.

    Verified rather than assumed: the antigen must have the expected number of
    residues and occupy a contiguous span. A's HLT files give offset 116 for a
    116-residue nanobody; Boltz gives 0.
    """
    ids = sorted(struct.get(ag_chain) or {})
    if not ids:
        raise SystemExit("no antigen chain %s" % ag_chain)
    if len(ids) != n_expected:
        raise SystemExit("antigen has %d residues, manifest says %d"
                         % (len(ids), n_expected))
    if ids[-1] - ids[0] + 1 != len(ids):
        raise SystemExit("antigen numbering is not contiguous: %d..%d over %d residues"
                         % (ids[0], ids[-1], len(ids)))
    return ids[0] - 1


def contacts_T(path, nb, ag, n_ag, cutoff, atoms):
    """Contacted antigen residues, renumbered to 1..N, plus the offset used."""
    st = fm.read_heavy(path, nb + ag, atoms=atoms)
    off = antigen_offset(st, ag, n_ag)
    return {r - off for r in fm.contact_antigen_residues(st, nb, ag, cutoff)}, off


def count_clashes(path, nb, ag, cutoff=CLASH):
    """Inter-chain heavy-atom pairs closer than `cutoff`."""
    st = fm.read_heavy(path, nb + ag)
    a = np.concatenate(list(st[nb].values()))
    n = 0
    for xyz in st[ag].values():
        d2 = ((xyz[:, None, :] - a[None, :, :]) ** 2).sum(-1)
        n += int((d2 < cutoff * cutoff).sum())
    return n


def buried_area(path, nb, ag):
    """Interface buried SASA in A^2, or None if it cannot be computed."""
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
    except ImportError:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = PDBParser(QUIET=True).get_structure("x", path)[0]
    for ch in list(model):
        if ch.id not in (nb, ag):
            model.detach_child(ch.id)
    # hydrogens would change the numbers between RF2 (explicit H) and the rest
    for ch in model:
        for res in ch:
            for at in [a for a in res if a.element == "H"]:
                res.detach_child(at.get_id())
    if len(list(model)) != 2:
        return None
    sr = ShrakeRupley()
    sr.compute(model, level="C")
    together = sum(c.sasa for c in model)
    alone = 0.0
    for keep in (nb, ag):
        solo = copy.deepcopy(model)
        for ch in list(solo):
            if ch.id != keep:
                solo.detach_child(ch.id)
        sr.compute(solo, level="C")
        alone += sum(c.sasa for c in solo)
    return (alone - together) / 2.0


def row_for(kind, path, nb, ag, n_ag, hotspots, epitope, cutoff):
    ep_bb, off = contacts_T(path, nb, ag, n_ag, cutoff, BACKBONE)
    ep_hv, _ = contacts_T(path, nb, ag, n_ag, cutoff, None)
    # Decided from the atoms present, not from `kind`: a backbone-only file has
    # exactly N, CA, C, O, so >5 heavy atoms per residue means side chains exist.
    st = fm.read_heavy(path, nb, atoms=None)
    side = sum(len(v) for v in st[nb].values()) / max(len(st[nb]), 1) > 5.0
    r = dict(structure_set=kind, has_sidechains=int(side), t_offset=off,
             cutoff=cutoff, clash_cutoff=CLASH,
             n_hotspot=len(hotspots), n_epitope=len(epitope),
             hs_contacted_bb=len(hotspots & ep_bb),
             hs_contacted_heavy=len(hotspots & ep_hv),
             ep_contacted_bb=len(epitope & ep_bb),
             ep_contacted_heavy=len(epitope & ep_hv),
             n_ag_contacts_bb=len(ep_bb), n_ag_contacts_heavy=len(ep_hv),
             n_clashes=count_clashes(path, nb, ag))
    r["hs_recovery_bb"] = round(r["hs_contacted_bb"] / len(hotspots), 4) if hotspots else ""
    r["hs_recovery_heavy"] = round(r["hs_contacted_heavy"] / len(hotspots), 4) if hotspots else ""
    r["ep_recovery_bb"] = round(r["ep_contacted_bb"] / len(epitope), 4) if epitope else ""
    r["ep_recovery_heavy"] = round(r["ep_contacted_heavy"] / len(epitope), 4) if epitope else ""
    bsa = buried_area(path, nb, ag) if side else None
    r["bsa_A2"] = "" if bsa is None else round(bsa, 1)
    return r


FIELDS = ["target", "pdb_id", "arm", "backbone_index", "candidate_id", "structure_set",
          "has_sidechains", "t_offset", "cutoff", "clash_cutoff",
          "n_hotspot", "hs_contacted_bb", "hs_recovery_bb",
          "hs_contacted_heavy", "hs_recovery_heavy",
          "n_epitope", "ep_contacted_bb", "ep_recovery_bb",
          "ep_contacted_heavy", "ep_recovery_heavy",
          "n_ag_contacts_bb", "n_ag_contacts_heavy", "n_clashes", "bsa_A2",
          "native_included"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true")
    p.add_argument("--run", default="runs/rfab_benchmark_1ce3344e4822")
    p.add_argument("--boltz", default="runs/boltz_out5")
    p.add_argument("--inputs", default="runs/boltz_in")
    p.add_argument("--stability", default="results/native_stability.csv")
    p.add_argument("--cutoff", type=float, default=fm.CUTOFF)
    p.add_argument("--out", default="results/interface_metrics.csv")
    a = p.parse_args()
    if a.demo:
        return demo()

    spec = json.load(open(os.path.join(a.run, "inputs", "targets.json")))
    spec = {d["pdb_id"]: d for d in (spec.values() if isinstance(spec, dict) else spec)}
    inc = {r["target"]: r["stable"] for r in csv.DictReader(open(a.stability))} \
        if os.path.exists(a.stability) else {}

    rows, miss = [], []
    for man in sorted(glob.glob(os.path.join(a.inputs, "*", "manifest.csv"))):
        tag = os.path.basename(os.path.dirname(man))
        for j in csv.DictReader(open(man)):
            if j["kind"] != "rfab":
                continue
            pdb_id, arm = tag.replace("_H_T", ""), j["arm"]
            idx = int(j["backbone"])
            n_ag = int(j["len_antigen"])
            sp = spec.get(pdb_id, {})
            hot = {int(h[1:]) for h in sp.get("hotspots", [])}
            epi = set(sp.get("epitope_T", []))
            gen, rf2 = candidate_paths(a.run, pdb_id, arm, idx)
            bz = boltz_path(a.boltz, tag, j["job"])
            for kind, path in (("generated", gen), ("rf2", rf2), ("boltz", bz)):
                if not path:
                    miss.append((pdb_id, arm, idx, kind))
                    continue
                nb, ag = NB_AG[kind]
                r = row_for(kind, path, nb, ag, n_ag, hot, epi, a.cutoff)
                r.update(target=tag, pdb_id=pdb_id, arm=arm, backbone_index=idx,
                         candidate_id=j["candidate_id"],
                         native_included=inc.get(pdb_id, ""))
                rows.append(r)
        print("  %s done (%d rows so far)" % (tag, len(rows)), flush=True)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print("\nwrote %s (%d rows)" % (a.out, len(rows)))
    if miss:
        print("missing structures: %d  e.g. %s" % (len(miss), miss[:3]))
    for kind in ("generated", "rf2", "boltz"):
        g = [r for r in rows if r["structure_set"] == kind]
        if g:
            print("  %-10s n=%3d  sidechains=%d  median hotspot recovery (bb)=%.3f"
                  % (kind, len(g), sum(r["has_sidechains"] for r in g),
                     float(np.median([r["hs_recovery_bb"] or 0 for r in g]))))


def demo():
    import tempfile
    # a 2-residue "nanobody" and a 3-residue "antigen" numbered A's way, so the
    # offset has to be found and removed before the hotspot set means anything
    def atom(i, name, ch, res, x):
        return ("ATOM  %5d  %-3s ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s"
                % (i, name, ch, res, x, 0.0, 0.0, name[0]))
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "t.pdb")
        ln, i = [], 1
        for res in (1, 2):                      # chain H at x=0
            for nm in BACKBONE:
                ln.append(atom(i, nm, "H", res, 0.0)); i += 1
        # antigen numbered 3,4,5 -> offset 2 -> T indices 1,2,3
        for k, res in enumerate((3, 4, 5)):
            for nm in BACKBONE:
                ln.append(atom(i, nm, "T", res, 3.0 if k == 0 else 50.0)); i += 1
        open(f, "w").write("\n".join(ln) + "\nEND\n")

        st = fm.read_heavy(f, "HT")
        assert antigen_offset(st, "T", 3) == 2, antigen_offset(st, "T", 3)
        # only the first antigen residue is within 5 A, and after the offset it is T1
        got, off = contacts_T(f, "H", "T", 3, 5.0, BACKBONE)
        assert got == {1}, got
        assert off == 2
        # the offset must be applied: without it the hit would read as residue 3
        assert 1 in got and 3 not in got
        # wrong expected length must fail loudly rather than shift the offset
        try:
            antigen_offset(st, "T", 4); raise AssertionError("should have failed")
        except SystemExit:
            pass
        # backbone-only input is detected as having no side chains
        r = row_for("generated", f, "H", "T", 3, {1}, {1, 2}, 5.0)
        assert r["has_sidechains"] == 0, r
        assert r["hs_contacted_bb"] == 1 and r["hs_recovery_bb"] == 1.0, r
        assert r["ep_contacted_bb"] == 1 and r["ep_recovery_bb"] == 0.5, r
        assert r["bsa_A2"] == "", r      # not reported without side chains
        assert r["n_clashes"] == 0, r
    print("demo ok")


if __name__ == "__main__":
    main()
