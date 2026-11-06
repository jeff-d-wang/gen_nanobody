#!/usr/bin/env python3
"""Add Chothia CDR remarks to C's prepped structures -> HLT files for ProteinMPNN.

C's prepped_mpnn PDBs are IMGT-numbered (max resnum 128, CDR3 insertions at
111A/112A) and carry no CDR annotation. A's RFantibody applies CHOTHIA spans
(H1 26-32, H2 52-56, H3 95-102). Scoring recovery on IMGT spans while A
generates on Chothia would mean Figure 3 described a stage the pipeline never
ran, so we renumber to Chothia ourselves with ANARCI and write the remarks.

Coordinates are copied byte for byte. Only `REMARK PDBinfo-LABEL` lines are
added, because C's files already fix a real ProteinMPNN bug (the parser
misorders insertion codes in 8qf4) and rewriting atoms would undo that.

  python3 src/evaluate/make_hlt.py --pdbs data/processed/prepped_mpnn \
      --out data/processed/hlt [--targets data/processed/targets.csv]

Needs ANARCI, in its own env. Installing it into a conda base env fails: the
solve tries to remove conda's own `requests`/`tqdm` and aborts.

  conda create -y -n nbnum -c bioconda -c conda-forge python=3.10 anarci hmmer biopython
  conda run -n nbnum python src/evaluate/make_hlt.py --pdbs ... --out ...

Only this script needs that env. The rest of the pipeline runs on base.
"""
import argparse
import csv
import glob
import os
import sys

# Chothia CDR spans, matching RFantibody's scripts/util/chothia2HLT.py.
# Testing the integer part alone is deliberate: it captures Chothia's insertion
# codes (31A/31B, 52A-52C, 100A-100K) without enumerating them.
CHOTHIA_CDR = {"H1": (26, 32), "H2": (52, 56), "H3": (95, 102)}

# Chothia-numbered VH anchors. A numbering shift moves ALL of them, so >=2
# mismatches means the numbering is wrong and every CDR index is wrong with it.
# Exactly one mismatch is a real substitution and only worth reporting: 8yj8
# genuinely carries R at 103 instead of the near-invariant Trp.
ANCHORS = {22: "C", 36: "W", 92: "C", 103: "W"}
MAX_BAD_ANCHORS = 1

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}


def chain_residues(lines, chain):
    """[(resname, resseq+icode)] for `chain`, in file order, one per residue."""
    out, seen = [], set()
    for ln in lines:
        if ln.startswith("ATOM") and ln[21] == chain and ln[12:16].strip() == "CA":
            key = ln[22:27]
            if key not in seen:
                seen.add(key)
                out.append((ln[17:20].strip(), key))
    return out


def chothia_map(seq, name):
    """[(pos, icode)] aligned 1:1 with `seq`, or None where ANARCI saw no domain."""
    try:
        from anarci import run_anarci
    except ImportError:
        raise SystemExit(
            "ANARCI not importable. It needs its own env (installing into a conda "
            "base env aborts on conda's own deps):\n"
            "  conda create -y -n nbnum -c bioconda -c conda-forge "
            "python=3.10 anarci hmmer biopython\n"
            "  conda run -n nbnum python %s ..." % os.path.basename(__file__))
    _, numbered, details, _ = run_anarci([(name, seq)], scheme="chothia")
    if not numbered[0]:
        raise SystemExit("%s: ANARCI found no antibody domain" % name)
    if len(numbered[0]) > 1:
        raise SystemExit("%s: ANARCI found %d domains, expected 1 VHH"
                         % (name, len(numbered[0])))
    ctype = details[0][0]["chain_type"]
    if ctype != "H":
        raise SystemExit("%s: ANARCI typed this chain %r, not H" % (name, ctype))

    numbering, start, end = numbered[0][0]
    out = [None] * len(seq)
    i = start
    for (pos, icode), aa in numbering:
        if aa == "-":
            continue
        if aa != seq[i]:
            raise SystemExit("%s: ANARCI residue %r at Chothia %d%s does not match "
                             "sequence residue %r at index %d"
                             % (name, aa, pos, icode.strip(), seq[i], i))
        out[i] = (pos, icode.strip())
        i += 1
    if i != end + 1:
        raise SystemExit("%s: consumed %d residues, ANARCI reported %d"
                         % (name, i - start, end - start + 1))

    got = {pc[0]: seq[k] for k, pc in enumerate(out) if pc and pc[0] in ANCHORS}
    bad = {p: (got[p], a) for p, a in ANCHORS.items() if p in got and got[p] != a}
    if len(bad) > MAX_BAD_ANCHORS:
        raise SystemExit("%s: %d Chothia anchors wrong, numbering is not Chothia "
                         "(pos: got/want) %s" % (name, len(bad), bad))
    missing = sorted(set(ANCHORS) - set(got))
    return out, missing, bad


def cdr_indices(cmap):
    """{H1|H2|H3: [0-based seq indices]} from a Chothia position map."""
    out = {}
    for cdr, (lo, hi) in CHOTHIA_CDR.items():
        out[cdr] = [i for i, pc in enumerate(cmap) if pc and lo <= pc[0] <= hi]
    return out


def write_hlt(lines, cdrs, dest):
    """Copy the PDB, inserting PDBinfo-LABEL remarks before END.

    HLT indices are 1-based and ABSOLUTE over chain order H,L,T. The heavy
    chain is first here, so absolute == per-chain; mpnn_cdr.py re-checks that
    against the parsed H length rather than trusting it.
    """
    remarks = ["REMARK PDBinfo-LABEL: %4d %s\n" % (i + 1, cdr)
               for cdr in ("H1", "H2", "H3") for i in cdrs[cdr]]
    # Keep the original END record rather than writing a fresh one: PDB lines are
    # padded to 80 columns and synthesizing it would be the one byte we changed.
    body = [ln for ln in lines if not ln.startswith("END")]
    end = [ln for ln in lines if ln.startswith("END")] or ["END\n"]
    with open(dest, "w") as f:
        f.writelines(body + remarks + end)


def demo():
    """Round-trip the remarks through mpnn_cdr.py, the only consumer that matters.

    Asserts the producer and the consumer agree on HLT indexing. If they ever
    drift, every CDR in Figure 3 is off by that drift and nothing else notices.
    """
    import tempfile

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mpnn_cdr import hlt_cdr_indices

    def atom(serial, chain, resnum):
        return ("ATOM  %5d  CA  ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n"
                % (serial, chain, resnum, 0.0, 0.0, 0.0))

    n_h, n_t = 120, 30
    lines = ["CRYST1\n"]
    lines += [atom(i + 1, "H", i + 1) for i in range(n_h)]
    lines += ["TER\n"]
    lines += [atom(n_h + i + 1, "T", i + 1) for i in range(n_t)]
    lines += ["TER\n", "END\n"]

    res = chain_residues(lines, "H")
    assert len(res) == n_h, "chain_residues lost residues: %d" % len(res)
    assert len(chain_residues(lines, "T")) == n_t, "antigen residues miscounted"

    # Chothia 1..120 with no insertions, so the spans land on known indices.
    cmap = [(i + 1, "") for i in range(n_h)]
    cdrs = cdr_indices(cmap)
    assert cdrs["H1"] == list(range(25, 32)), cdrs["H1"]     # 26-32 -> 0-based 25-31
    assert cdrs["H2"] == list(range(51, 56)), cdrs["H2"]     # 52-56
    assert cdrs["H3"] == list(range(94, 102)), cdrs["H3"]    # 95-102

    # an insertion inside a span must be picked up by the integer-part test
    ins = list(cmap)
    ins[100] = (100, "A")
    assert 100 in cdr_indices(ins)["H3"], "Chothia 100A must count as CDR3"

    with tempfile.TemporaryDirectory() as d:
        dest = os.path.join(d, "T1.pdb")
        write_hlt(lines, cdrs, dest)
        out = open(dest).readlines()
        assert out[-1] == "END\n" and out.count("END\n") == 1, "END must be last and single"
        assert sum(1 for ln in out if ln.startswith("ATOM")) == n_h + n_t, "atoms changed"
        got = hlt_cdr_indices(dest, n_h)

    assert got == {"cdr1": cdrs["H1"], "cdr2": cdrs["H2"], "cdr3": cdrs["H3"]}, got

    # remarks pointing past the H chain must be rejected, not silently applied
    # to the antigen
    with tempfile.TemporaryDirectory() as d:
        dest = os.path.join(d, "BAD.pdb")
        write_hlt(lines, {"H1": [0], "H2": [1], "H3": [n_h + 5]}, dest)
        try:
            hlt_cdr_indices(dest, n_h)
            raise AssertionError("out-of-range remark must raise")
        except SystemExit:
            pass

    print("demo ok")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true", help="run self-check and exit")
    p.add_argument("--pdbs", help="dir of C's prepped H/T PDBs")
    p.add_argument("--out", help="dir to write HLT files into")
    p.add_argument("--targets", help="targets.csv, to cross-check against C's IMGT spans")
    p.add_argument("--heavy", default="H")
    a = p.parse_args()

    if a.demo:
        return demo()
    if not (a.pdbs and a.out):
        p.error("--pdbs and --out are required")

    imgt = {}
    if a.targets:
        for r in csv.DictReader(open(a.targets)):
            imgt[r["PDB_ID"]] = (r["CDRH1"], r["CDRH2"], r["CDRH3"])

    os.makedirs(a.out, exist_ok=True)
    pdbs = sorted(glob.glob(os.path.join(a.pdbs, "*.pdb")))
    if not pdbs:
        raise SystemExit("no PDBs in %s" % a.pdbs)

    print("%-16s %4s %5s %5s %5s  %s" % ("target", "nH", "H1", "H2", "H3", "notes"))
    tot = 0
    for pdb in pdbs:
        tag = os.path.splitext(os.path.basename(pdb))[0]
        lines = open(pdb).readlines()
        res = chain_residues(lines, a.heavy)
        if not res:
            raise SystemExit("%s: no chain %s" % (tag, a.heavy))
        seq = "".join(AA3.get(rn, "X") for rn, _ in res)

        cmap, missing, bad = chothia_map(seq, tag)
        cdrs = cdr_indices(cmap)
        for cdr, idx in cdrs.items():
            if not idx:
                raise SystemExit("%s: %s came back empty" % (tag, cdr))

        notes = []
        if missing:
            notes.append("anchors absent from model: %s" % missing)
        for p, (g, w) in sorted(bad.items()):
            notes.append("substitution at Chothia %d: %s not %s" % (p, g, w))
        # C annotated these independently under IMGT. Chothia H3 should sit
        # inside C's IMGT CDRH3; a disjoint span means one of us is wrong.
        pid = tag.replace("_H_T", "")
        if pid in imgt:
            i3 = seq.find(imgt[pid][2])
            if i3 < 0:
                notes.append("C's IMGT CDRH3 not in this structure")
            else:
                span = set(range(i3, i3 + len(imgt[pid][2])))
                ov = len(span & set(cdrs["H3"])) / len(cdrs["H3"])
                notes.append("H3 overlap with C's IMGT %.0f%%" % (100 * ov))

        write_hlt(lines, cdrs, os.path.join(a.out, tag + ".pdb"))
        n = sum(len(v) for v in cdrs.values())
        tot += n
        print("%-16s %4d %5d %5d %5d  %s"
              % (pid, len(seq), len(cdrs["H1"]), len(cdrs["H2"]), len(cdrs["H3"]),
                 "; ".join(notes)))

    print("\nwrote %d HLT files to %s" % (len(pdbs), a.out))
    print("mean designed positions per target: %.1f" % (tot / len(pdbs)))
    print("next: src/evaluate/run_mpnn_cdr.sh %s" % a.out)


if __name__ == "__main__":
    main()
