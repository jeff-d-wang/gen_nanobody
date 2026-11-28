#!/usr/bin/env python3
"""Independent re-check of the RQ2 feature-rank outputs.

Reproduces the reported macro/micro concordance from the per-comparison audit,
then adds three analyses the original run did not report: the dose-response
across the pre-specified KD separation thresholds, percentile bootstrap
intervals over groups, and a McNemar test between the paratope-side and
epitope-side contact counts on the comparisons both of them scored. Also checks
whether the paratope count is a proxy for sequence coverage.

The pre-declared protocol withholds p-values for a single concordance value
because comparisons inside a group are dependent. The bootstrap here resamples
groups, and the McNemar test is paired on identical comparisons, so neither
treats within-group comparisons as independent.

Usage:
  python3 rq2_recheck.py --root <nanobody_rq2_clean> --out <results dir>

stdlib only.
"""
import argparse
import csv
import math
import os
import random
import statistics as st
from collections import defaultdict

OUTCOME = {"correct": 1.0, "incorrect": 0.0, "tie": 0.5}
PRIMARY = "structure_paper_proxy_EXPLORATORY"
SENS = "antigen_only_mixed_assays_SENSITIVITY"
METRICS = ["vhh_residues", "ag_residues", "atom_contacts", "residue_contacts", "clashes"]


def spearman(xs, ys):
    n = len(xs)
    def rk(v):
        o = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[o[j + 1]] == v[o[i]]:
                j += 1
            a = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[o[k]] = a
            i = j + 1
        return r
    rx, ry = rk(xs), rk(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((p - mx) * (q - my) for p, q in zip(rx, ry))
    den = math.sqrt(sum((p - mx) ** 2 for p in rx) * sum((q - my) ** 2 for q in ry))
    return num / den if den else float("nan")


def mcnemar_exact(b, c):
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--draws", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=37)
    a = ap.parse_args()

    rep = os.path.join(a.root, "reports", "rq2_feature_rank_v1")
    groups = list(csv.DictReader(open(os.path.join(rep, "group_pairwise_metrics.csv"))))
    comps = list(csv.DictReader(open(os.path.join(rep, "informative_pair_comparisons.csv"))))
    scores = [r for r in csv.DictReader(open(os.path.join(a.root, "tables", "structure_scores.csv")))
              if r["mapping_status"] == "MAPPED_CONTACT"]
    os.makedirs(a.out, exist_ok=True)

    for r in groups:
        r["thr"] = float(r["kd_fold_threshold"])
        r["n"] = int(r["n_KD_informative_comparisons"])
        r["ties"] = int(r["score_ties"])
        r["c"] = float(r["concordance"]) if r["concordance"] else None

    rng = random.Random(a.seed)
    out = []
    for tier in (PRIMARY, SENS):
        for thr in (1.5, 3.0, 10.0):
            for m in METRICS:
                sub = [r for r in groups if r["tier"] == tier and r["thr"] == thr
                       and r["metric"] == m and r["n"] > 0]
                if not sub:
                    continue
                per = [r["c"] for r in sub]
                n = len(per)
                macro = st.mean(per)
                n_comp = sum(r["n"] for r in sub)
                micro = sum(r["c"] * r["n"] for r in sub) / n_comp
                lo = hi = float("nan")
                if n >= 2:
                    d = sorted(st.mean(rng.choice(per) for _ in range(n)) for _ in range(a.draws))
                    lo, hi = d[int(0.025 * a.draws)], d[int(0.975 * a.draws)]
                out.append({
                    "tier": tier, "kd_fold_threshold": thr, "metric": m,
                    "n_groups": n, "n_comparisons": n_comp,
                    "n_ties": sum(r["ties"] for r in sub),
                    "macro_concordance": round(macro, 4),
                    "micro_concordance": round(micro, 4),
                    "boot_ci_low": round(lo, 4), "boot_ci_high": round(hi, 4),
                })

    cols = ["tier", "kd_fold_threshold", "metric", "n_groups", "n_comparisons",
            "n_ties", "macro_concordance", "micro_concordance",
            "boot_ci_low", "boot_ci_high"]
    p1 = os.path.join(a.out, "rq2_recheck_concordance.csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)

    pv, pa = {}, {}
    for r in comps:
        if r["tier"] != PRIMARY:
            continue
        key = (r["antigen_id"], r["structure_study_proxy"], r["pair_i"], r["pair_j"])
        if r["metric"] == "vhh_residues":
            pv[key] = OUTCOME[r["outcome"]]
        elif r["metric"] == "ag_residues":
            pa[key] = OUTCOME[r["outcome"]]
    both = sorted(set(pv) & set(pa))
    b = sum(1 for k in both if pv[k] > pa[k])
    c = sum(1 for k in both if pv[k] < pa[k])
    p_mcnemar = mcnemar_exact(b, c)

    by_pair = defaultdict(list)
    for r in scores:
        by_pair[r["pair_id"]].append(r)

    def med(rs, col):
        v = [float(x[col]) for x in rs if x[col] not in ("", "nan")]
        return st.median(v) if v else None

    xs, ys = [], []
    for rs in by_pair.values():
        p, q = med(rs, "vhh_residues"), med(rs, "vhh_coverage")
        if p is not None and q is not None:
            xs.append(p)
            ys.append(q)
    rho_cov = spearman(xs, ys)
    covered = sum(1 for q in ys if q >= 0.95) / len(ys)

    xs2, ys2 = [], []
    for rs in by_pair.values():
        p, q = med(rs, "vhh_residues"), med(rs, "ag_residues")
        if p is not None and q is not None:
            xs2.append(p)
            ys2.append(q)
    rho_sides = spearman(xs2, ys2)

    p2 = os.path.join(a.out, "rq2_recheck_controls.csv")
    with open(p2, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["check", "value", "n", "note"])
        w.writerow(["mcnemar_vhh_beats_ag", b, len(both), "VHH correct where antigen wrong"])
        w.writerow(["mcnemar_ag_beats_vhh", c, len(both), "antigen correct where VHH wrong"])
        w.writerow(["mcnemar_exact_p", round(p_mcnemar, 4), len(both), "paired, primary tier, 3x"])
        w.writerow(["spearman_vhhres_vhhcoverage", round(rho_cov, 4), len(xs),
                    "near zero rules out a coverage artefact"])
        w.writerow(["frac_vhh_coverage_ge_0.95", round(covered, 4), len(ys), ""])
        w.writerow(["spearman_vhhres_agres", round(rho_sides, 4), len(xs2),
                    "the two sides measure largely the same extent"])

    print("reproduced primary tier, 3x separation:")
    print("%-20s %8s %8s %22s %7s" % ("metric", "macro", "micro", "95% CI (groups)", "ties"))
    for r in out:
        if r["tier"] == PRIMARY and r["kd_fold_threshold"] == 3.0:
            print("%-20s %8.3f %8.3f      [%.3f, %.3f] %7d" %
                  (r["metric"], r["macro_concordance"], r["micro_concordance"],
                   r["boot_ci_low"], r["boot_ci_high"], r["n_ties"]))
    print("\ndose response, primary tier, macro concordance:")
    print("%-20s %8s %8s %8s" % ("metric", "1.5x", "3.0x", "10x"))
    for m in METRICS:
        v = []
        for thr in (1.5, 3.0, 10.0):
            hit = [r for r in out if r["tier"] == PRIMARY
                   and r["kd_fold_threshold"] == thr and r["metric"] == m]
            v.append(hit[0]["macro_concordance"] if hit else float("nan"))
        print("%-20s %8.3f %8.3f %8.3f" % (m, *v))
    print("\nMcNemar, paratope vs epitope count on %d shared comparisons: %d vs %d, p = %.4f"
          % (len(both), b, c, p_mcnemar))
    print("Spearman(vhh_residues, vhh_coverage) = %.3f over %d pairs; %.1f%% cover >= 0.95"
          % (rho_cov, len(xs), 100 * covered))
    print("Spearman(vhh_residues, ag_residues)  = %.3f over %d pairs" % (rho_sides, len(xs2)))
    print("\nwrote %s\n      %s" % (p1, p2))


if __name__ == "__main__":
    main()
