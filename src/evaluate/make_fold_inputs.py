#!/usr/bin/env python3
"""Build Boltz input YAMLs: native, B's MPNN redesigns, and A's RFantibody candidates.

Boltz is the excluded final scorer (PLAN.md), so it folds from SEQUENCE ONLY,
with no template and no designed backbone. That is also what the proposal
requires: "predict complexes without supplying the generated CDR backbone or
binding pose as a template".

The nanobody sequence comes from ProteinMPNN's output FASTA (designed H chain
only, no antigen). The antigen comes from parsed.jsonl, so it is byte-identical
across every fold. That identity is what makes hotspot-contact recovery a plain
set comparison on chain-T indices, with no structural alignment anywhere.

A's candidates are a different molecule, not a CDR variant of the native nanobody:
a de novo loop design (H1 7, H2 6, H3 8) on the upstream h-NbBCII10 framework, 116 aa
against the native's 126. So there is no native to diff against and no shared numbering,
and the native HLT files must not be used for their CDR indices. Each candidate carries
its own REMARK PDBinfo-LABEL lines; those are what `--from-rfab` reads, and the manifest
now records per-job CDR indices so one Table 2 can hold both kinds of nanobody.

What makes the two comparable at all: the chain-T antigen is byte-identical between A's
run inputs and B's parsed.jsonl on all 21 targets, numbered 1..N (verified, see PLAN.md).
That is re-checked here, because every contact metric downstream rests on it.

  # B's MPNN redesigns of the native backbone (native + top-N by global_score)
  python3 src/evaluate/make_fold_inputs.py --target pdb_00008qf4 \
      --run runs/mpnn_seed37 --hlt data/processed/hlt --out runs/boltz_in

  # all 84 of A's candidates, merged into the same per-target dirs
  python3 src/evaluate/make_fold_inputs.py --from-rfab runs/rfab_benchmark_1ce3344e4822

Self-check:  python3 src/evaluate/make_fold_inputs.py --demo
"""
import argparse
import csv
import json
import os
import re
import sys

NANOBODY_CHAIN = "A"   # Boltz chain id for the designed VHH (our chain H)
ANTIGEN_CHAIN = "B"    # Boltz chain id for the target  (our chain T)

# One schema for both builders, so they can merge into one manifest per target.
ROW_FIELDS = ["job", "target", "kind", "sample", "global_score", "arm", "backbone",
              "candidate_id", "nanobody_chain", "antigen_chain", "cdr_idx",
              "n_cdr_positions", "n_mutations", "len_nanobody", "len_antigen",
              "msa", "nanobody_seq"]


def read_jsonl(path):
    with open(path) as f:
        return [json.loads(ln) for ln in f if ln.strip()]


def load_designs(fa_path, n_top):
    """(native_seq, [(sample, global_score, seq)]) - best n_top by global_score.

    Lower global_score is better. Record 0 is the native; the rest are samples.
    """
    from Bio import SeqIO

    recs = list(SeqIO.parse(fa_path, "fasta"))
    if len(recs) < 2:
        raise SystemExit("%s has %d records, expected native + samples" % (fa_path, len(recs)))
    native = str(recs[0].seq).replace("/", "")

    out = []
    for r in recs[1:]:
        seq = str(r.seq).replace("/", "")
        m_s = re.search(r"sample=(\d+)", r.description)
        m_g = re.search(r"global_score=([0-9.]+)", r.description)
        if not (m_s and m_g):
            raise SystemExit("cannot parse sample/global_score from %r" % r.description[:80])
        out.append((int(m_s.group(1)), float(m_g.group(1)), seq))
    out.sort(key=lambda x: x[1])
    return native, out[:n_top]


def yaml_for(nanobody, antigen, msa, ag_msa=None):
    """Boltz v1 input.

    msa="server" omits the key (pair with `boltz --use_msa_server`), "empty" is
    single-sequence for both chains, and "cached" points the antigen at a prefetched
    MSA while keeping the nanobody single-sequence. "cached" exists because the 21
    antigens repeat across 84 designs, so one fetch each replaces 105, and because
    server hits for a designed VHH are framework noise either way.
    """
    if msa == "cached" and not ag_msa:
        raise SystemExit("msa=cached needs an antigen MSA path")

    def block(cid, seq, m):
        s = "  - protein:\n      id: %s\n      sequence: %s\n" % (cid, seq)
        if m:
            s += "      msa: %s\n" % m
        return s

    nb_m = "empty" if msa in ("empty", "cached") else None
    ag_m = "empty" if msa == "empty" else (ag_msa if msa == "cached" else None)
    return "version: 1\nsequences:\n" + block(NANOBODY_CHAIN, nanobody, nb_m) \
        + block(ANTIGEN_CHAIN, antigen, ag_m)


def write_manifest(dest, rows):
    """Merge `rows` into dest/manifest.csv by job name.

    Merging rather than overwriting is what lets the MPNN and RFantibody builders both
    contribute jobs for one target, so `fold_metrics.py report` emits a single table
    covering the native, B's redesigns and A's candidates under the same scorer.
    Native sorts first: it is the control, and the runner folds in manifest order.
    """
    path = os.path.join(dest, "manifest.csv")
    keep = {}
    if os.path.exists(path):
        old = list(csv.DictReader(open(path)))
        if old and set(old[0]) != set(ROW_FIELDS):
            print("note: %s has an older schema, replacing it" % path)
        else:
            keep = {r["job"]: r for r in old}
    keep.update({r["job"]: r for r in rows})
    def order(j):
        r = keep[j]
        return (r["kind"] != "native", float(r["global_score"] or 0), j)

    out = [keep[j] for j in sorted(keep, key=order)]

    # the identity that the whole metric design rests on
    assert len({str(r["len_antigen"]) for r in out}) == 1, \
        "antigen length varies across folds in %s" % dest
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_FIELDS)
        w.writeheader()
        w.writerows(out)
    return out


def row(**kw):
    r = {k: "" for k in ROW_FIELDS}
    r.update(nanobody_chain=NANOBODY_CHAIN, antigen_chain=ANTIGEN_CHAIN)
    r.update(kw)
    bad = set(r) - set(ROW_FIELDS)
    assert not bad, "unknown manifest fields %s" % bad
    return r


def build(target, run_dir, hlt_dir, out_dir, n_top, msa, msa_dir=None):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mpnn_cdr import hlt_cdr_indices

    tag = target if target.endswith("_H_T") else target + "_H_T"
    parsed = {e["name"]: e for e in read_jsonl(os.path.join(run_dir, "parsed.jsonl"))}
    if tag not in parsed:
        raise SystemExit("%s not in %s/parsed.jsonl (have %d targets)"
                         % (tag, run_dir, len(parsed)))
    seq_h = parsed[tag]["seq_chain_H"]
    seq_t = parsed[tag]["seq_chain_T"]

    native, designs = load_designs(os.path.join(run_dir, "seqs", tag + ".fa"), n_top)
    if native != seq_h:
        raise SystemExit("%s: FASTA native H does not match parsed.jsonl seq_chain_H" % tag)

    # Designs must differ from the native ONLY inside the CDRs. If a design has
    # drifted in the framework the mask was wrong, and every metric computed
    # downstream would describe a different experiment than the one we claim.
    cdrs = hlt_cdr_indices(os.path.join(hlt_dir, tag + ".pdb"), len(seq_h))
    free = {i for v in cdrs.values() for i in v}
    for sample, score, seq in designs:
        if len(seq) != len(seq_h):
            raise SystemExit("%s sample=%d: length %d != native %d"
                             % (tag, sample, len(seq), len(seq_h)))
        drift = [i for i in range(len(seq)) if i not in free and seq[i] != seq_h[i]]
        if drift:
            raise SystemExit("%s sample=%d: %d framework positions differ (first at %d)"
                             % (tag, sample, len(drift), drift[0]))

    dest = os.path.join(out_dir, tag)
    os.makedirs(dest, exist_ok=True)
    jobs = [("native", 0, None, native)] + \
           [("design_s%d" % s, s, g, q) for s, g, q in designs]
    cdr_idx = ";".join(str(i) for i in sorted(free))

    rows = []
    for name, sample, score, seq in jobs:
        with open(os.path.join(dest, name + ".yaml"), "w") as f:
            f.write(yaml_for(seq, seq_t, msa, msa_path(msa_dir, tag)))
        n_mut = sum(1 for i in range(len(seq)) if seq[i] != seq_h[i])
        rows.append(row(job=name, target=tag, kind="native" if sample == 0 else "design",
                        sample=sample, global_score="" if score is None else "%.4f" % score,
                        cdr_idx=cdr_idx, n_cdr_positions=len(free), n_mutations=n_mut,
                        len_nanobody=len(seq), len_antigen=len(seq_t),
                        msa=msa, nanobody_seq=seq))

    merged = write_manifest(dest, rows)
    print("%-20s %-8s %8s %10s %s" % ("job", "kind", "sample", "glob_score", "CDR muts"))
    for r in rows:
        print("%-20s %-8s %8s %10s %s/%s"
              % (r["job"], r["kind"], r["sample"], r["global_score"] or "-",
                 r["n_mutations"], r["n_cdr_positions"]))
    print("\nwrote %d YAMLs to %s (manifest now holds %d jobs)"
          % (len(rows), dest, len(merged)))
    print("antigen chain %s is byte-identical across all folds (%d aa)"
          % (ANTIGEN_CHAIN, len(seq_t)))
    print("msa mode: %s" % msa)
    return rows


def msa_path(msa_dir, tag):
    """Where the antigen MSA lives, as Boltz will see it on the folding host.

    Deliberately relative: the YAMLs are built here and rsynced to EC2, so an
    absolute local path would not resolve there. run_boltz.sh cds to the input root.
    ponytail: flat <pdb_id>.csv naming, fine while one antigen == one target.

    .csv, not .a3m: with --use_msa_server Boltz writes the MERGED msa as
    boltz_results_<job>/msa/<job>_<chain_idx>.csv (header key,sequence). The raw
    a3m files it leaves in msa/<job>_unpaired_tmp_env/ are the separate MMseqs2
    databases (uniref, bfd/env), unmerged, so concatenating them is NOT the same
    input. boltz/main.py dispatches on suffix and parses .a3m or .csv, so the
    harvested CSV is reused directly. Verified Nov 24 on 8qf4: cached CSV folded
    with zero server queries, iptm 0.940 against the known interface.
    """
    return None if not msa_dir else "%s/%s.csv" % (msa_dir.rstrip("/"), tag.replace("_H_T", ""))


def pdb_chain_seq(path, chain):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from make_hlt import AA3, chain_residues

    lines = open(path).read().splitlines(True)
    res = chain_residues(lines, chain)
    if not res:
        raise SystemExit("%s has no chain %s" % (path, chain))
    return "".join(AA3.get(r, "X") for r, _ in res)


def load_rfab(run_dir):
    """[(pdb_id, arm, label, seq_h, cdr_idx)] for A's accepted candidates.

    Candidate identity is (pdb_id, arm, backbone index). The accepted stage directory
    is read out of each state.json rather than inferred from a path, because A's
    handoff warns that the accepted attempt is not always attempt 1 and an interrupted
    attempt was left on disk. Globbing the run tree would count abandoned files.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mpnn_cdr import hlt_cdr_indices

    name = os.path.basename(os.path.normpath(run_dir))

    def local(p):
        # state.json stores absolute paths from the machine that ran the generation
        return os.path.join(run_dir, p.split(name + "/", 1)[1]) if name + "/" in p else p

    fasta, hdr = {}, None
    for ln in open(os.path.join(run_dir, "designed_sequences.fasta")):
        if ln.startswith(">"):
            hdr = ln[1:].strip()
            fasta[hdr] = ""
        elif hdr:
            fasta[hdr] += ln.strip()

    out, frameworks = [], set()
    for d in sorted(x for x in os.listdir(run_dir) if x.startswith("pdb_")):
        sp = os.path.join(run_dir, d, "state.json")
        if not os.path.isfile(sp):
            continue
        st = json.load(open(sp))
        for o in st["stages"]["proteinmpnn"]["outputs"]:
            f = local(o["file"])
            label = os.path.basename(f)[:-len(".pdb")]
            seq_h = pdb_chain_seq(f, "H")
            cdrs = hlt_cdr_indices(f, len(seq_h))      # the candidate's OWN labels
            idx = sorted({i for v in cdrs.values() for i in v})
            cid = "%s_%s_%s" % (st["pdb_id"], st["arm"], label)
            if fasta.get(cid) != seq_h:
                raise SystemExit("%s: designed_sequences.fasta does not match chain H of %s"
                                 % (cid, f))
            frameworks.add("".join(c for i, c in enumerate(seq_h) if i not in idx))
            out.append((st["pdb_id"], st["arm"], label, seq_h, idx))

    if not out:
        raise SystemExit("no accepted ProteinMPNN outputs under %s" % run_dir)
    # Fixed-framework invariant, checkable without any alignment: RFantibody designs
    # only H1/H2/H3, so stripping each candidate's own CDR positions must leave one
    # shared framework string. A's audit checks this upstream; a corrupted or
    # mis-paired handoff would show up right here.
    if len(frameworks) != 1:
        raise SystemExit("%d distinct framework sequences across %d candidates - the "
                         "framework should be fixed" % (len(frameworks), len(out)))
    print("%d candidates, framework fixed at %d aa across all of them"
          % (len(out), len(next(iter(frameworks)))))
    return out


def build_rfab(run_dir, mpnn_run, out_dir, msa, msa_dir):
    """One YAML per accepted RFantibody candidate, merged into the per-target dirs."""
    cands = load_rfab(run_dir)

    pj = os.path.join(mpnn_run, "parsed.jsonl")
    ref = {}
    if os.path.exists(pj):
        ref = {e["name"].replace("_H_T", ""): e["seq_chain_T"] for e in read_jsonl(pj)}
    else:
        print("note: no %s, skipping the antigen cross-check" % pj)

    by_target = {}
    for pid, arm, label, seq_h, idx in cands:
        by_target.setdefault(pid, []).append((arm, label, seq_h, idx))

    total = 0
    for pid in sorted(by_target):
        tag = pid + "_H_T"
        seq_t = pdb_chain_seq(os.path.join(run_dir, "inputs", "targets", pid + "_T.pdb"), "T")
        # Every contact metric assumes antigen residue i is the same residue in A's
        # inputs and in B's. Check it, never assume it (PLAN.md, Nov 16).
        if pid in ref and ref[pid] != seq_t:
            raise SystemExit("%s: antigen differs between A's run inputs and %s. Contact "
                             "recovery on shared T indices is invalid until this is "
                             "resolved." % (pid, pj))

        dest = os.path.join(out_dir, tag)
        os.makedirs(dest, exist_ok=True)
        rows = []
        for arm, label, seq_h, idx in sorted(by_target[pid]):
            bb = re.search(r"nb_(\d+)", label)
            job = "rfab_%s_nb%s" % (arm, bb.group(1) if bb else label)
            with open(os.path.join(dest, job + ".yaml"), "w") as f:
                f.write(yaml_for(seq_h, seq_t, msa, msa_path(msa_dir, tag)))
            rows.append(row(job=job, target=tag, kind="rfab", arm=arm,
                            backbone=bb.group(1) if bb else "",
                            candidate_id="%s_%s_%s" % (pid, arm, label),
                            cdr_idx=";".join(str(i) for i in idx),
                            n_cdr_positions=len(idx), len_nanobody=len(seq_h),
                            len_antigen=len(seq_t), msa=msa, nanobody_seq=seq_h))
        merged = write_manifest(dest, rows)
        total += len(rows)
        print("%-14s %2d rfab jobs, antigen %3d aa, manifest %2d jobs"
              % (pid, len(rows), len(seq_t), len(merged)))

    print("\nwrote %d candidate YAMLs across %d targets under %s"
          % (total, len(by_target), out_dir))
    print("msa mode: %s" % msa)
    return total


def demo():
    """Asserts the invariants that would silently invalidate every fold metric."""
    import tempfile

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from make_hlt import write_hlt

    seq_h = "AAACDRONECCCCCDRTWOGGGGGGCDRTHREEWWW"
    seq_t = "YYYYYMMMM"
    spans = {"H1": "CDRONE", "H2": "CDRTWO", "H3": "CDRTHREE"}
    cdrs = {k: list(range(seq_h.index(v), seq_h.index(v) + len(v))) for k, v in spans.items()}
    free = sorted(i for v in cdrs.values() for i in v)

    with tempfile.TemporaryDirectory() as d:
        run, hlt, out = os.path.join(d, "run"), os.path.join(d, "hlt"), os.path.join(d, "out")
        os.makedirs(os.path.join(run, "seqs"))
        os.makedirs(hlt)
        tag = "T1_H_T"

        with open(os.path.join(run, "parsed.jsonl"), "w") as f:
            f.write(json.dumps({"name": tag, "seq_chain_H": seq_h, "seq_chain_T": seq_t}) + "\n")

        lines = ["ATOM      1  CA  ALA H   1       0.0   0.0   0.0\n", "END\n"]
        write_hlt(lines, cdrs, os.path.join(hlt, tag + ".pdb"))

        good = seq_h[:free[0]] + "Q" + seq_h[free[0] + 1:]          # one CDR mutation
        with open(os.path.join(run, "seqs", tag + ".fa"), "w") as f:
            f.write(">%s, score=1.0, global_score=1.40\n%s\n" % (tag, seq_h))
            f.write(">T=0.1, sample=1, score=0.5, global_score=1.20\n%s\n" % good)
            f.write(">T=0.1, sample=2, score=0.6, global_score=1.10\n%s\n" % seq_h)

        rows = build(tag, run, hlt, out, 2, "empty")

        assert rows[0]["kind"] == "native" and rows[0]["n_mutations"] == 0
        # best global_score must come first among designs
        assert rows[1]["sample"] == 2 and rows[1]["global_score"] == "1.1000", rows[1]
        assert rows[2]["sample"] == 1 and rows[2]["n_mutations"] == 1
        assert len({r["len_antigen"] for r in rows}) == 1

        y = open(os.path.join(out, tag, "native.yaml")).read()
        assert "id: %s\n      sequence: %s" % (NANOBODY_CHAIN, seq_h) in y, y
        assert "id: %s\n      sequence: %s" % (ANTIGEN_CHAIN, seq_t) in y, y
        assert y.count("msa: empty") == 2, "single-seq mode must mark both chains"
        assert y.index(NANOBODY_CHAIN + "\n") < y.index(ANTIGEN_CHAIN + "\n"), \
            "nanobody chain must come first so Boltz chain order matches H,T"

        # per-job CDR indices: fold_metrics.py reads these instead of the native HLT,
        # because A's generated chains are a different length with their own labels
        assert rows[0]["cdr_idx"] == ";".join(str(i) for i in free), rows[0]["cdr_idx"]
        man = list(csv.DictReader(open(os.path.join(out, tag, "manifest.csv"))))
        assert [r["job"] for r in man] == [r["job"] for r in rows], man
        assert list(man[0]) == ROW_FIELDS, "manifest schema drifted"

        # cached MSA: antigen reads a prefetched MSA, designed nanobody stays single-seq
        y2 = yaml_for("AAA", "CCC", "cached", msa_path("msa", tag))
        assert "msa: msa/T1.csv" in y2, y2
        assert y2.count("msa: empty") == 1, y2
        assert msa_path("", tag) is None

        # framework drift must be refused, not folded
        with open(os.path.join(run, "seqs", tag + ".fa"), "w") as f:
            f.write(">%s, score=1.0, global_score=1.40\n%s\n" % (tag, seq_h))
            f.write(">T=0.1, sample=1, score=0.5, global_score=1.20\n%s\n" % ("Q" + seq_h[1:]))
        try:
            build(tag, run, hlt, out, 1, "empty")
            raise AssertionError("framework drift must raise")
        except SystemExit as e:
            assert "framework positions differ" in str(e), e

    print("demo ok")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true", help="run self-check and exit")
    p.add_argument("--target", help="e.g. pdb_00008qf4")
    p.add_argument("--from-rfab", metavar="RUNDIR",
                   help="A's RFantibody run dir; folds every accepted candidate")
    p.add_argument("--run", default="runs/mpnn_seed37", help="ProteinMPNN run dir")
    p.add_argument("--hlt", default="data/processed/hlt")
    p.add_argument("--out", default="runs/boltz_in")
    p.add_argument("-n", "--n-top", type=int, default=4, help="designs to fold besides native")
    p.add_argument("--msa", choices=("server", "empty", "cached"), default="server",
                   help="server: pair with boltz --use_msa_server. empty: single-sequence. "
                        "cached: antigen reads --msa-dir/<pdb_id>.csv, nanobody single-seq")
    p.add_argument("--msa-dir", default="msa",
                   help="antigen MSA dir AS SEEN ON THE FOLDING HOST, relative to the "
                        "input root (default: msa)")
    a = p.parse_args()

    if a.demo:
        return demo()
    if a.from_rfab:
        return build_rfab(a.from_rfab, a.run, a.out, a.msa, a.msa_dir)
    if not a.target:
        p.error("need --target or --from-rfab")
    build(a.target, a.run, a.hlt, a.out, a.n_top, a.msa, a.msa_dir)


if __name__ == "__main__":
    main()
