#!/usr/bin/env python3
"""Build RFdiffusion target inputs from C's benchmark.

Two things RFantibody needs that are not on main:

  1. A target-only PDB. C's prepped_v2 files carry the native nanobody in chain
     H alongside the antigen in chain T. RFdiffusion generates a new H, so the
     native one must not be in the target file.
  2. ppi.hotspot_res in the numbering of that target file. C's
     epitope_residues_v2.csv is in ORIGINAL antigen numbering (chain A, residue
     211); the crop renumbered the merged antigen to T1..Tn. The translation
     lives in antigen_residue_mapping_v2.csv (output_res_id).

Hotspots are a subset, not the whole epitope: a 25-residue epitope passed whole
over-constrains the dock, and RFantibody's own examples use 3-6. Selection is
deterministic -- seed with the most buried epitope residue (smallest
min_distance_to_nanobody_A), then greedily add the next most buried residue
whose CA is at least HOTSPOT_MIN_SEP from every residue already chosen, so the
hotspots spread across the interface instead of clustering in one pocket. The
full epitope is still written to the epitope_all_T column.

stdlib only, so it runs under any python3 on PSC without an env.
"""
import argparse
import csv
import math
import os

HOTSPOT_MIN_SEP = 8.0


def parse_t_ca(pdb_path):
    coords = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[21] != "T" or line[12:16].strip() != "CA":
                continue
            coords[int(line[22:26])] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def dist(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def pick_hotspots(cand, coords, n):
    """cand is epitope residue ids in T numbering, most buried first."""
    chosen = [cand[0]]
    for resid in cand[1:]:
        if len(chosen) >= n:
            break
        if resid not in coords:
            continue
        seps = [dist(coords[resid], coords[c]) for c in chosen if c in coords]
        if not seps or min(seps) >= HOTSPOT_MIN_SEP:
            chosen.append(resid)
    # If the separation rule could not fill n slots, top up by burial order
    # rather than returning fewer hotspots than asked for.
    for resid in cand[1:]:
        if len(chosen) >= n:
            break
        if resid not in chosen:
            chosen.append(resid)
    return sorted(chosen[:n])


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed", default="data/processed")
    ap.add_argument("--out-dir", default="data/processed/rfab_inputs")
    ap.add_argument("--n-hotspots", type=int, default=4)
    args = ap.parse_args()

    proc = args.processed
    epi = read_rows(os.path.join(proc, "epitope_residues_v2.csv"))
    amap = read_rows(os.path.join(proc, "antigen_residue_mapping_v2.csv"))
    prep = {r["PDB_ID"]: r for r in read_rows(os.path.join(proc, "prepped_v2.csv"))}
    qc = {r["PDB_ID"]: r for r in read_rows(os.path.join(proc, "structural_qc_final.csv"))}

    def key(r):
        return (r["PDB_ID"], r["original_chain"], r["original_res_id"],
                r["original_insertion_code"] or "")

    to_out = {key(r): int(r["output_res_id"]) for r in amap}

    by_pdb = {}
    unmapped = 0
    for r in epi:
        out_id = to_out.get(key(r))
        if out_id is None:
            unmapped += 1
            continue
        by_pdb.setdefault(r["PDB_ID"], []).append(
            (float(r["min_distance_to_nanobody_A"]), out_id))

    os.makedirs(os.path.join(args.out_dir, "targets"), exist_ok=True)
    rows = []
    for pdb_id in sorted(by_pdb):
        src = os.path.join(proc, "prepped_v2", pdb_id + "_H_T.pdb")
        coords = parse_t_ca(src)

        dst = os.path.join(args.out_dir, "targets", pdb_id + "_T.pdb")
        with open(src) as fi, open(dst, "w") as fo:
            for line in fi:
                if line.startswith(("ATOM", "HETATM")) and line[21] != "T":
                    continue
                if line.startswith(("CRYST1", "ATOM", "HETATM", "TER", "END")):
                    fo.write(line)

        cand = [rid for _, rid in sorted(by_pdb[pdb_id])]
        hs = pick_hotspots(cand, coords, args.n_hotspots)
        q = qc.get(pdb_id, {})
        rows.append({
            "PDB_ID": pdb_id,
            "target_pdb": os.path.relpath(dst, args.out_dir),
            "T_residues": len(coords),
            "native_H_residues": prep.get(pdb_id, {}).get("nanobody_residues", ""),
            "epitope_residues": len(cand),
            "hotspot_res": "[" + ",".join("T%d" % i for i in hs) + "]",
            "epitope_all_T": " ".join("T%d" % i for i in cand),
            "QC": q.get("QC", ""),
            "QC_issues": q.get("issues", ""),
        })

    rows.sort(key=lambda r: r["T_residues"])
    for i, r in enumerate(rows, start=1):
        r["size_rank"] = i

    cols = ["size_rank", "PDB_ID", "target_pdb", "T_residues", "native_H_residues",
            "epitope_residues", "hotspot_res", "epitope_all_T", "QC", "QC_issues"]
    out_csv = os.path.join(args.out_dir, "rfab_targets.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    print("unmapped epitope residues: %d" % unmapped)
    print("targets written: %d" % len(rows))
    print()
    print("%4s %-14s %5s %5s  %-24s %s" %
          ("rank", "PDB", "T", "epi", "hotspot_res", "QC"))
    for r in rows:
        print("%4d %-14s %5d %5d  %-24s %s" %
              (r["size_rank"], r["PDB_ID"], r["T_residues"],
               r["epitope_residues"], r["hotspot_res"], r["QC"]))
    print()
    print("wrote %s" % out_csv)
    if unmapped:
        raise SystemExit("STOP: %d epitope residues did not map to a T residue id"
                         % unmapped)


if __name__ == "__main__":
    main()
