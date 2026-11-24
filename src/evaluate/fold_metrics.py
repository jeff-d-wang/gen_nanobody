#!/usr/bin/env python3
"""Metrics for the end-to-end target: pLDDT, interchain PAE, hotspot-contact recovery.

Chain T residue numbering is identical across prepped_v2/, prepped_mpnn/ and
data/processed/hlt/, and the antigen sequence is unchanged between the native
and every design. So antigen residue i in a folded model IS antigen residue i
in the native structure, and contact recovery is a set comparison on T indices.
No structural alignment, no mapping tables.

  validate  reproduce C's 5A epitope set from the HLT file  (run this first)
  report    native + designs -> results/end_to_end_metrics.csv and Table 2

  python3 src/evaluate/fold_metrics.py validate --target pdb_00008qf4
  python3 src/evaluate/fold_metrics.py report --target pdb_00008qf4 \
      --boltz runs/boltz_out --inputs runs/boltz_in

Self-check:  python3 src/evaluate/fold_metrics.py --demo
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np

CUTOFF = 5.0          # matches EPITOPE_CUTOFF in src/data/04_prep_structures_v2.py
NANOBODY_CHAIN = "A"  # in Boltz output; the native HLT file uses H
ANTIGEN_CHAIN = "B"   # in Boltz output; the native HLT file uses T


def read_heavy(path, chains):
    """{chain: {resseq: (n_atoms, 3) float array}} of heavy-atom coordinates.

    Hydrogens are excluded so the cutoff means the same thing as C's, which was
    computed on heavy atoms only.
    """
    if not os.path.exists(path):
        raise SystemExit("no structure at %s" % path)
    out = {c: {} for c in chains}
    with open(path) as f:
        for ln in f:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            ch = ln[21]
            if ch not in out:
                continue
            el = ln[76:78].strip() or ln[12:16].strip()[:1]
            if el.upper() in ("H", "D"):
                continue
            try:
                res = int(ln[22:26])
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
            out[ch].setdefault(res, []).append(xyz)
    return {c: {r: np.asarray(v, dtype=float) for r, v in d.items()} for c, d in out.items()}


def contact_antigen_residues(struct, nb_chain, ag_chain, cutoff=CUTOFF):
    """Antigen residue ids within `cutoff` A (heavy atom) of any nanobody residue."""
    nb = struct.get(nb_chain) or {}
    ag = struct.get(ag_chain) or {}
    if not nb or not ag:
        raise SystemExit("missing chain: nanobody %s has %d residues, antigen %s has %d"
                         % (nb_chain, len(nb), ag_chain, len(ag)))
    nb_xyz = np.concatenate(list(nb.values()))
    hit = set()
    for res, xyz in ag.items():
        # (n_ag_atoms, n_nb_atoms) distances for this residue only; keeps memory flat
        d2 = ((xyz[:, None, :] - nb_xyz[None, :, :]) ** 2).sum(-1)
        if d2.min() <= cutoff * cutoff:
            hit.add(res)
    return hit


def native_epitope(hlt_path, cutoff=CUTOFF):
    return contact_antigen_residues(read_heavy(hlt_path, "HT"), "H", "T", cutoff)


def charlie_epitope(pdb_id, proc="data/processed"):
    """C's epitope for this target, expressed as chain-T indices.

    epitope_residues_v2.csv is keyed by ORIGINAL chain and residue number, so it
    has to come through antigen_residue_mapping_v2.csv to reach T indices.
    Chain ids collide across namespaces (9dc3's original chains are T, U, e), so
    PDB_ID is always part of the key.
    """
    key = {}
    with open(os.path.join(proc, "antigen_residue_mapping_v2.csv")) as f:
        for r in csv.DictReader(f):
            if r["PDB_ID"] == pdb_id:
                key[(r["original_chain"], r["original_res_id"],
                     r["original_insertion_code"])] = int(r["output_res_id"])
    out, missing = set(), []
    with open(os.path.join(proc, "epitope_residues_v2.csv")) as f:
        for r in csv.DictReader(f):
            if r["PDB_ID"] != pdb_id:
                continue
            k = (r["original_chain"], r["original_res_id"], r["original_insertion_code"])
            if k in key:
                out.add(key[k])
            else:
                missing.append(k)
    if missing:
        raise SystemExit("%s: %d epitope residues absent from the mapping table: %s"
                         % (pdb_id, len(missing), missing[:3]))
    return out


def boltz_outputs(job_dir):
    """(confidence dict, plddt array, pae matrix, structure path) for one Boltz job."""
    def one(pattern, what):
        hits = sorted(glob.glob(os.path.join(job_dir, "**", pattern), recursive=True))
        if not hits:
            raise SystemExit("%s: no %s (looked for %s)" % (job_dir, what, pattern))
        return hits[0]

    conf = json.load(open(one("confidence_*model_0.json", "confidence json")))
    plddt = np.load(one("plddt_*model_0.npz", "plddt npz"))
    plddt = plddt[plddt.files[0]]
    pae = np.load(one("pae_*model_0.npz", "pae npz"))
    pae = pae[pae.files[0]]
    struct = sorted(glob.glob(os.path.join(job_dir, "**", "*model_0.pdb"), recursive=True))
    if not struct:
        raise SystemExit("%s: no *model_0.pdb - run boltz with --output_format pdb" % job_dir)
    return conf, np.asarray(plddt).ravel(), np.asarray(pae), struct[0]


def metrics_for(job_dir, n_nb, cdr_idx, native_ep, cutoff=CUTOFF):
    conf, plddt, pae, struct_path = boltz_outputs(job_dir)
    n_tot = len(plddt)
    if n_tot != pae.shape[0]:
        raise SystemExit("%s: plddt has %d entries, pae is %s" % (job_dir, n_tot, pae.shape))
    if n_tot <= n_nb:
        raise SystemExit("%s: %d tokens but the nanobody alone is %d residues; "
                         "is the antigen chain missing?" % (job_dir, n_tot, n_nb))
    # Boltz reports pLDDT on 0-1; the report quotes the usual 0-100 scale.
    if plddt.max() <= 1.0:
        plddt = plddt * 100.0

    nb = slice(0, n_nb)
    cdr_mask = np.zeros(n_nb, dtype=bool)
    cdr_mask[sorted(cdr_idx)] = True

    inter = np.concatenate([pae[:n_nb, n_nb:].ravel(), pae[n_nb:, :n_nb].ravel()])
    st = read_heavy(struct_path, NANOBODY_CHAIN + ANTIGEN_CHAIN)
    got = contact_antigen_residues(st, NANOBODY_CHAIN, ANTIGEN_CHAIN, cutoff)

    return dict(
        plddt_mean="%.1f" % plddt.mean(),
        plddt_cdr="%.1f" % plddt[nb][cdr_mask].mean(),
        plddt_framework="%.1f" % plddt[nb][~cdr_mask].mean(),
        pae_interchain="%.2f" % inter.mean(),
        pae_interchain_min="%.2f" % inter.min(),
        iptm="%.3f" % conf["iptm"] if "iptm" in conf else "",
        ptm="%.3f" % conf["ptm"] if "ptm" in conf else "",
        n_native_epitope=len(native_ep),
        n_design_contacts=len(got),
        n_recovered=len(native_ep & got),
        # fraction of the NATIVE epitope contacted. Not A's 4-residue hotspot set:
        # those are a separate, smaller target recorded in A's inputs/targets.json.
        epitope_recovery="%.3f" % (len(native_ep & got) / len(native_ep)) if native_ep else "",
    )


def cmd_validate(a):
    tag = a.target if a.target.endswith("_H_T") else a.target + "_H_T"
    pdb_id = tag.replace("_H_T", "")
    mine = native_epitope(os.path.join(a.hlt, tag + ".pdb"), a.cutoff)
    theirs = charlie_epitope(pdb_id, a.processed)
    only_mine, only_theirs = sorted(mine - theirs), sorted(theirs - mine)
    print("%s at %.1f A" % (pdb_id, a.cutoff))
    print("  mine   %3d chain-T residues" % len(mine))
    print("  C's    %3d chain-T residues" % len(theirs))
    print("  shared %3d" % len(mine & theirs))
    if only_mine:
        print("  only mine:   %s" % only_mine)
    if only_theirs:
        print("  only C's:    %s" % only_theirs)
    if mine == theirs:
        print("MATCH: our contact code reproduces C's epitope exactly")
    else:
        raise SystemExit("MISMATCH - do not report metrics until this is understood")


def cmd_report(a):
    tag = a.target if a.target.endswith("_H_T") else a.target + "_H_T"
    man_path = os.path.join(a.inputs, tag, "manifest.csv")
    if not os.path.exists(man_path):
        raise SystemExit("no manifest at %s - run make_fold_inputs.py first" % man_path)
    jobs = list(csv.DictReader(open(man_path)))
    if jobs and "cdr_idx" not in jobs[0]:
        raise SystemExit("%s predates per-job CDR indices - regenerate it with "
                         "make_fold_inputs.py" % man_path)

    # The epitope is the one thing shared by every job: chain-T numbering is identical
    # across A's inputs and B's, so the native epitope is the common reference set.
    native_ep = native_epitope(os.path.join(a.hlt, tag + ".pdb"), a.cutoff)

    rows = []
    for j in jobs:
        job_dir = os.path.join(a.boltz, tag, j["job"])
        if not os.path.isdir(job_dir):
            print("skip %-20s (no output at %s)" % (j["job"], job_dir))
            continue
        # Per job, never from the native HLT: an RFantibody candidate is a 116-aa
        # chain on a different framework with its own CDR labels, so native indices
        # would point at the wrong residues (and run off the end of the chain).
        n_nb = int(j["len_nanobody"])
        cdr_idx = [int(i) for i in j["cdr_idx"].split(";") if i]
        m = metrics_for(job_dir, n_nb, cdr_idx, native_ep, a.cutoff)
        m.update(job=j["job"], kind=j["kind"], sample=j["sample"], arm=j["arm"],
                 candidate_id=j["candidate_id"], len_nanobody=n_nb,
                 mpnn_global_score=j["global_score"], n_cdr_mutations=j["n_mutations"],
                 target=tag, predictor="boltz", role="final scorer, excluded from "
                 "generation and selection")
        rows.append(m)
    if not rows:
        raise SystemExit("no Boltz outputs found under %s" % os.path.join(a.boltz, tag))

    order = ["target", "job", "kind", "arm", "candidate_id", "sample", "mpnn_global_score",
             "len_nanobody", "n_cdr_mutations",
             "plddt_mean", "plddt_framework", "plddt_cdr", "pae_interchain",
             "pae_interchain_min", "iptm", "ptm", "n_native_epitope", "n_design_contacts",
             "n_recovered", "epitope_recovery", "predictor", "role"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=order)
        w.writeheader()
        w.writerows(rows)

    hdr = ["job", "kind", "pLDDT", "fw", "CDR", "iPAE", "ipTM", "epi_recov"]
    print("\nTable 2: end-to-end metrics for %s (Boltz, no template)" % tag)
    print("%-20s %-7s %6s %6s %6s %7s %6s %9s" % tuple(hdr))
    for r in rows:
        print("%-20s %-7s %6s %6s %6s %7s %6s %9s"
              % (r["job"], r["kind"], r["plddt_mean"], r["plddt_framework"], r["plddt_cdr"],
                 r["pae_interchain"], r["iptm"], r["epitope_recovery"]))
    print("\nnative epitope: %d chain-T residues at %.1f A" % (len(native_ep), a.cutoff))
    print("wrote %s (%d rows)" % (a.out, len(rows)))


def demo():
    """Contacts on hand-placed coordinates, where the answer is known by construction."""
    import tempfile

    def atom(serial, chain, res, xyz, el="C"):
        return ("ATOM  %5d  CA  ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n"
                % (serial, chain, res, xyz[0], xyz[1], xyz[2], el))

    # nanobody atom at the origin; antigen residues at 3, 4.5, 6 and 20 A
    lines = [atom(1, "H", 1, (0, 0, 0))]
    for i, d in enumerate([3.0, 4.5, 6.0, 20.0], start=1):
        lines.append(atom(10 + i, "T", i, (d, 0, 0)))
    lines.append(atom(99, "T", 99, (1.0, 0, 0), el="H"))   # hydrogen, must be ignored
    lines.append("END\n")

    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x.pdb")
        open(p, "w").writelines(lines)
        st = read_heavy(p, "HT")
        assert len(st["T"]) == 4, "hydrogen-only residue must be dropped: %s" % sorted(st["T"])
        got = contact_antigen_residues(st, "H", "T", 5.0)
        assert got == {1, 2}, got                    # 3.0 and 4.5 are within 5 A
        assert contact_antigen_residues(st, "H", "T", 6.5) == {1, 2, 3}

    # recovery arithmetic
    native, design = {1, 2, 3, 4}, {2, 3, 9}
    assert len(native & design) / len(native) == 0.5

    # interchain PAE must take both off-diagonal blocks, not just one
    n_nb, n = 2, 5
    pae = np.zeros((n, n))
    pae[:n_nb, n_nb:] = 10.0
    pae[n_nb:, :n_nb] = 20.0
    inter = np.concatenate([pae[:n_nb, n_nb:].ravel(), pae[n_nb:, :n_nb].ravel()])
    assert abs(inter.mean() - 15.0) < 1e-9, inter.mean()

    print("demo ok")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true", help="run self-check and exit")
    p.add_argument("--hlt", default="data/processed/hlt")
    p.add_argument("--processed", default="data/processed")
    p.add_argument("--cutoff", type=float, default=CUTOFF)
    sub = p.add_subparsers(dest="cmd")

    v = sub.add_parser("validate", help="check our contacts against C's epitope set")
    v.add_argument("--target", required=True)
    v.set_defaults(func=cmd_validate)

    r = sub.add_parser("report", help="build Table 2 from Boltz outputs")
    r.add_argument("--target", required=True)
    r.add_argument("--boltz", default="runs/boltz_out")
    r.add_argument("--inputs", default="runs/boltz_in")
    r.add_argument("--out", default="results/end_to_end_metrics.csv")
    r.set_defaults(func=cmd_report)

    a = p.parse_args()
    if a.demo:
        return demo()
    if not a.cmd:
        p.error("need a subcommand (validate, report) or --demo")
    a.func(a)


if __name__ == "__main__":
    main()
