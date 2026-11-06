#!/usr/bin/env python3
"""ProteinMPNN CDR redesign: build fixed-position masks, then score CDR recovery.

Both subcommands read ProteinMPNN's OWN parsed jsonl, so residue indices can
never drift from what the model actually saw (parsing the PDB separately is how
you get a silently wrong mask):

  mask     parsed.jsonl + --hlt DIR          -> fixed_positions.jsonl
  recover  parsed.jsonl + --hlt DIR + seqs/  -> cdr_recovery.csv

CDRs come from the `REMARK PDBinfo-LABEL: <abs_idx> H1|H2|H3` lines that
RFantibody's scripts/util/chothia2HLT.py writes, so A's RFdiffusion and B's
recovery scoring consume one artifact with one Chothia definition
(H1 26-32, H2 52-56, H3 95-102).

Self-check:  python3 mpnn_cdr.py demo
"""
import argparse
import csv
import glob
import json
import os
import re
from collections import defaultdict

from Bio import SeqIO

HLT_LABELS = {"H1": "cdr1", "H2": "cdr2", "H3": "cdr3"}


def find_once(hay, needle, what):
    """Index of `needle` in `hay`, erroring unless it occurs exactly once."""
    i = hay.find(needle)
    if i < 0:
        raise ValueError("%s: %r not found" % (what, needle))
    if hay.find(needle, i + 1) >= 0:
        raise ValueError("%s: %r occurs more than once, cannot place it" % (what, needle))
    return i


def read_jsonl(path):
    with open(path) as f:
        return [json.loads(ln) for ln in f if ln.strip()]


def hlt_cdr_indices(pdb_path, n_heavy):
    """{cdrN: [0-based H-chain indices]} from HLT REMARK PDBinfo-LABEL lines.

    HLT stores 1-indexed ABSOLUTE residue indices over chain order H,L,T, so for
    the heavy chain absolute == per-chain. That is asserted, not assumed: if the
    heavy chain is not first the indices would silently point into the antigen.
    """
    if not os.path.exists(pdb_path):
        raise SystemExit("no HLT file at %s - run chothia2HLT.py on C's Chothia PDBs" % pdb_path)
    out = defaultdict(list)
    with open(pdb_path) as f:
        for ln in f:
            if not ln.startswith("REMARK PDBinfo-LABEL:"):
                continue
            parts = ln.split(":", 1)[1].split()
            if len(parts) >= 2 and parts[1] in HLT_LABELS:
                out[HLT_LABELS[parts[1]]].append(int(parts[0]) - 1)
    if not out:
        raise SystemExit("%s has no H1/H2/H3 PDBinfo-LABEL remarks - is it HLT-formatted?" % pdb_path)
    for cdr, idx in out.items():
        if max(idx) >= n_heavy:
            raise SystemExit(
                "%s: %s index %d exceeds the %d-residue H chain. HLT requires chain "
                "order H,L,T - check chothia2HLT.py output." % (pdb_path, cdr, max(idx) + 1, n_heavy))
    return {k: sorted(v) for k, v in out.items()}


def chains_of(entry):
    return sorted(k[len("seq_chain_"):] for k in entry if k.startswith("seq_chain_"))


def build_mask(entries, heavy, hlt_dir):
    """{name: {chain: [1-indexed fixed positions]}} - everything fixed but the CDRs."""
    out = {}
    for entry in entries:
        name = entry["name"]
        key = "seq_chain_%s" % heavy
        if key not in entry:
            raise SystemExit("%s has no chain %s (chains: %s)" % (name, heavy, chains_of(entry)))
        seq_h = entry[key]
        cdrs = hlt_cdr_indices(os.path.join(hlt_dir, name + ".pdb"), len(seq_h))
        designed = {i for idx in cdrs.values() for i in idx}
        fixed = {}
        for c in chains_of(entry):
            n = len(entry["seq_chain_%s" % c])
            keep = range(n) if c != heavy else (i for i in range(n) if i not in designed)
            fixed[c] = [i + 1 for i in keep]
        out[name] = fixed
    return out


def recovery_rows(entries, seqs_dir, heavy, hlt_dir):
    parsed = {e["name"]: e for e in entries}
    rows = []
    for fa in sorted(glob.glob(os.path.join(seqs_dir, "*.fa"))):
        name = os.path.splitext(os.path.basename(fa))[0]
        if name not in parsed:
            continue
        recs = list(SeqIO.parse(fa, "fasta"))   # native first, then samples
        if len(recs) < 2:
            continue
        seq_h = parsed[name]["seq_chain_%s" % heavy]
        native = str(recs[0].seq).replace("/", "")
        off = find_once(native, seq_h, "%s native H chain" % name)
        cdrs = hlt_cdr_indices(os.path.join(hlt_dir, name + ".pdb"), len(seq_h))
        for rec in recs[1:]:
            m = re.search(r"sample=(\d+)", rec.description)
            sample = int(m.group(1)) if m else len(rows)
            des_h = str(rec.seq).replace("/", "")[off:off + len(seq_h)]
            if len(des_h) != len(seq_h):
                raise ValueError("%s sample=%s length mismatch" % (name, sample))
            pooled_n = pooled_id = 0
            for cdr in sorted(cdrs):
                idx = cdrs[cdr]
                nat = "".join(seq_h[i] for i in idx)
                des = "".join(des_h[i] for i in idx)
                nid = sum(a == b for a, b in zip(nat, des))
                pooled_n += len(idx)
                pooled_id += nid
                rows.append(dict(target=name, sample=sample, cdr=cdr, n_pos=len(idx),
                                 n_identical=nid, recovery=nid / len(idx),
                                 native=nat, designed=des))
            rows.append(dict(target=name, sample=sample, cdr="all", n_pos=pooled_n,
                             n_identical=pooled_id, recovery=pooled_id / pooled_n,
                             native="", designed=""))
    return rows


def cmd_mask(a):
    entries = read_jsonl(a.parsed)
    mask = build_mask(entries, a.heavy, a.hlt)
    if not mask:
        raise SystemExit("no targets in %s" % a.parsed)
    with open(a.out, "w") as f:
        f.write(json.dumps(mask) + "\n")
    n_des = sum(len(e["seq_chain_%s" % a.heavy]) - len(mask[e["name"]][a.heavy]) for e in entries)
    print("wrote %s: %d targets, %d designed positions" % (a.out, len(mask), n_des))


def cmd_recover(a):
    rows = recovery_rows(read_jsonl(a.parsed), a.seqs, a.heavy, a.hlt)
    if not rows:
        raise SystemExit("no sequences matched in %s" % a.seqs)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for cdr in ("cdr1", "cdr2", "cdr3", "all"):
        v = [r["recovery"] for r in rows if r["cdr"] == cdr]
        if v:
            print("%-5s n=%3d  mean recovery %.3f" % (cdr, len(v), sum(v) / len(v)))
    print("wrote %s (%d rows)" % (a.out, len(rows)))


def demo(_a=None):
    """Asserts the things that would silently corrupt the experiment."""
    import tempfile
    seq_h = "AAACDRONECCCCCDRTWOGGGGGGCDRTHREEWWW"
    seq_t = "YYYYY"
    entry = {"name": "T1", "seq_chain_H": seq_h, "seq_chain_T": seq_t}
    spans = {"cdr1": "CDRONE", "cdr2": "CDRTWO", "cdr3": "CDRTHREE"}
    want = {k: list(range(seq_h.index(v), seq_h.index(v) + len(v))) for k, v in spans.items()}

    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "T1.pdb"), "w") as f:
            f.write("ATOM      1  N   ALA H   1       0.0   0.0   0.0\n")
            for label, cdr in (("H1", "cdr1"), ("H2", "cdr2"), ("H3", "cdr3")):
                for i in want[cdr]:
                    f.write("REMARK PDBinfo-LABEL: %4d %s\n" % (i + 1, label))
        assert hlt_cdr_indices(os.path.join(d, "T1.pdb"), len(seq_h)) == want

        # remarks pointing past the H chain must fail loudly, not silently
        # design the antigen
        with open(os.path.join(d, "BAD.pdb"), "w") as f:
            f.write("REMARK PDBinfo-LABEL: %4d H3\n" % (len(seq_h) + 3))
        try:
            hlt_cdr_indices(os.path.join(d, "BAD.pdb"), len(seq_h))
            raise AssertionError("out-of-range HLT index must raise")
        except SystemExit:
            pass

        mask = build_mask([entry], "H", d)["T1"]
        free = {i for v in want.values() for i in v}
        assert mask["H"] == [i + 1 for i in range(len(seq_h)) if i not in free]
        assert len(mask["T"]) == len(seq_t), "antigen must be fully fixed"
        assert min(mask["H"]) >= 1, "must be 1-indexed"

        # recovery: the offset must follow the H chain even when it is second
        mut = seq_h.replace("CDRTHREE", "CDRTHREF")          # 7/8 on cdr3
        with open(os.path.join(d, "T1.fa"), "w") as f:
            f.write(">T1, score=0\n%s/%s\n" % (seq_t, seq_h))
            f.write(">T=0.1, sample=1, score=0\n%s/%s\n" % (seq_t, mut))
        got = {r["cdr"]: r for r in recovery_rows([entry], d, "H", d)}

    assert got["cdr1"]["recovery"] == 1.0
    assert got["cdr3"]["n_identical"] == 7 and got["cdr3"]["n_pos"] == 8
    assert got["all"]["n_identical"] == 19 and got["all"]["n_pos"] == 20
    print("demo ok")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--heavy", default="H", help="nanobody chain id")
    p.add_argument("--hlt", help="dir of HLT-formatted PDBs")
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("mask", help="build fixed_positions jsonl")
    m.add_argument("--parsed", required=True)
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_mask)

    r = sub.add_parser("recover", help="score CDR recovery")
    r.add_argument("--parsed", required=True)
    r.add_argument("--seqs", required=True, help="ProteinMPNN out_folder/seqs")
    r.add_argument("--out", default="results/cdr_recovery.csv")
    r.set_defaults(func=cmd_recover)

    sub.add_parser("demo", help="run self-check").set_defaults(func=demo)

    a = p.parse_args()
    if a.cmd != "demo" and not a.hlt:
        p.error("--hlt DIR is required")
    a.func(a)
