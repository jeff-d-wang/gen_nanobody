#!/usr/bin/env python3
"""Generation status for the midterm report, read off runs/rfab/manifest.csv.

Writes results/generation_status.csv and prints the paragraph numbers. Counts
backbones on disk rather than trusting num_designs, so a partial run reports
what it actually produced. stdlib only.
"""
import argparse
import csv
import glob
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", default="results/generation_status.csv")
    args = ap.parse_args()

    man = os.path.join(args.repo, "runs", "rfab", "manifest.csv")
    if not os.path.exists(man):
        raise SystemExit("no manifest at %s -- nothing has run yet" % man)

    rows = list(csv.DictReader(open(man)))
    # Keep the last attempt per (target, arm); earlier ones were skipped or failed.
    latest = {}
    for r in rows:
        latest[(r["pdb_id"], r["arm"])] = r

    out_rows = []
    for (pid, arm), r in sorted(latest.items()):
        d = os.path.join(args.repo, "runs", "rfab", "%s_%s" % (pid, arm))
        on_disk = len([p for p in glob.glob(os.path.join(d, "*.pdb"))])
        out_rows.append({
            "pdb_id": pid,
            "arm": arm,
            "requested": int(r["num_designs"]),
            "backbones_on_disk": on_disk,
            "wall_seconds": int(r["wall_seconds"]),
            "exit_code": int(r["exit_code"]),
            "job_id": r["job_id"],
            "utc": r["utc"],
        })

    os.makedirs(os.path.dirname(os.path.join(args.repo, args.out)) or ".", exist_ok=True)
    cols = ["pdb_id", "arm", "requested", "backbones_on_disk", "wall_seconds",
            "exit_code", "job_id", "utc"]
    with open(os.path.join(args.repo, args.out), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)

    ok = [r for r in out_rows if r["exit_code"] == 0 and r["backbones_on_disk"] > 0]
    tot_bb = sum(r["backbones_on_disk"] for r in ok)
    tot_s = sum(r["wall_seconds"] for r in ok)
    targets = sorted({r["pdb_id"] for r in ok})
    arms = sorted({r["arm"] for r in ok})
    failed = [r for r in out_rows if r["exit_code"] != 0 or r["backbones_on_disk"] == 0]

    print("runs recorded          %d" % len(out_rows))
    print("runs with backbones    %d" % len(ok))
    print("targets with backbones %d / 21" % len(targets))
    print("arms represented       %s" % (", ".join(arms) or "none"))
    print("backbones on disk      %d" % tot_bb)
    print("GPU wall time          %.2f h (%d s)" % (tot_s / 3600.0, tot_s))
    if tot_bb:
        print("mean per backbone      %.1f s" % (tot_s / tot_bb))
    if failed:
        print("\nfailed or empty:")
        for r in failed:
            print("  %s %-8s exit %s, %d pdb" %
                  (r["pdb_id"], r["arm"], r["exit_code"], r["backbones_on_disk"]))
    print("\nwrote %s" % os.path.join(args.out))


if __name__ == "__main__":
    main()
