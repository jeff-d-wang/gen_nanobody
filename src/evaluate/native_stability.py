#!/usr/bin/env python3
"""Native-control quality, and the target inclusion rule it defines.

A target is INCLUDED only if its native control passes both:

    min iptm over the 5 diffusion samples >= 0.7     (the fold is reliable)
    native epitope recovery              >= 0.5      (the fold is in the RIGHT PLACE)

Both are properties of the target and the predictor, never of any candidate, so
neither leaks anything about the designs. Full rationale in the Nov 24 addendum
of PLAN.md; the short version:

  - The iptm condition exists because at 1 diffusion sample the weak targets vary
    by 0.25-0.47 iptm between draws while the 8qf4 control varies by 0.023, so a
    single draw is not a statistic. Five criteria (min>=0.7, median>=0.7,
    min>=0.6, model_0>=0.7 AND spread<=0.10, median>=0.7 AND spread<=0.10) pick
    the same 9 targets, so that stratification is a property of the data.
    `model_0 >= 0.7` alone is the one to avoid: it admits 8rwf (model_0 0.722,
    min 0.236), 9h77 (0.726, 0.333) and 9ljj (0.830, 0.364).

  - The recovery condition exists because iptm and placement are decoupled IN
    BOTH DIRECTIONS. 9g3e folds confidently at the wrong interface: min iptm
    0.904, spread 0.045, and it contacts 1 of 18 annotated epitope residues.
    9g2g does the reverse, min iptm 0.481 with recovery 0.917. Confidence of any
    kind -- plDDT, iptm -- measures how sure the predictor is about an interface,
    not whether it is THE interface. Only epitope recovery measures placement.

Ordering, which non-negotiable #3 cares about: the iptm condition was frozen
before any design number was read. The recovery condition was added AFTER the
first report run, on native-only evidence. It is therefore not blind, and the
report must say so rather than present the final target count as the original
plan.

  python3 src/evaluate/native_stability.py --boltz runs/boltz_out5
  python3 src/evaluate/native_stability.py --demo
"""
import argparse
import csv
import glob
import json
import os
import statistics as st

import fold_metrics as fm

FLOOR_IPTM = 0.7
FLOOR_RECOVERY = 0.5


def sample_iptms(job_dir):
    """{sample_index: iptm} over every diffusion sample in one job directory."""
    out = {}
    for f in glob.glob(os.path.join(job_dir, "**", "confidence_*model_*.json"),
                       recursive=True):
        try:
            idx = int(f.split("model_")[1].split(".json")[0])
        except (IndexError, ValueError):
            continue
        out[idx] = json.load(open(f))["iptm"]
    return out


def native_recovery(boltz_dir, tag, hlt_dir, cutoff=fm.CUTOFF):
    """(n_epitope, recovered_fraction) for a target's native control, or None.

    Reuses fold_metrics' contact code rather than reimplementing it, so the
    epitope definition here cannot drift from the one `report` uses. That code
    is validated against C's independently derived epitope on all 21 targets.
    """
    hlt = os.path.join(hlt_dir, tag + ".pdb")
    struct = sorted(glob.glob(os.path.join(boltz_dir, tag, "native", "**",
                                           "*model_0.pdb"), recursive=True))
    if not os.path.exists(hlt) or not struct:
        return None
    ep = fm.native_epitope(hlt, cutoff)
    if not ep:
        return None
    got = fm.contact_antigen_residues(
        fm.read_heavy(struct[0], fm.NANOBODY_CHAIN + fm.ANTIGEN_CHAIN),
        fm.NANOBODY_CHAIN, fm.ANTIGEN_CHAIN, cutoff)
    return len(ep), len(ep & got) / len(ep)


def rows_for(boltz_dir, recovery_fn):
    """One row per target. recovery_fn(tag) -> (n_epitope, fraction) or None."""
    rows = []
    for d in sorted(glob.glob(os.path.join(boltz_dir, "*_H_T"))):
        tag = os.path.basename(d)
        v = sample_iptms(os.path.join(d, "native"))
        if not v:
            continue
        vals = sorted(v.values())
        rec = recovery_fn(tag)
        n_ep, recov = rec if rec else ("", None)
        iptm_ok = vals[0] >= FLOOR_IPTM
        rec_ok = recov is not None and recov >= FLOOR_RECOVERY
        rows.append(dict(
            target=tag.replace("_H_T", ""), n_samples=len(vals),
            model_0=round(v.get(0, float("nan")), 4),
            iptm_min=round(vals[0], 4), iptm_median=round(st.median(vals), 4),
            iptm_max=round(vals[-1], 4), spread=round(vals[-1] - vals[0], 4),
            n_epitope=n_ep,
            epitope_recovery="" if recov is None else round(recov, 4),
            iptm_ok=int(iptm_ok), recovery_ok=int(rec_ok),
            stable=int(iptm_ok and rec_ok)))
    return sorted(rows, key=lambda r: (-r["stable"], -r["iptm_min"]))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true", help="run self-check and exit")
    p.add_argument("--boltz", default="runs/boltz_out5")
    p.add_argument("--hlt", default="data/processed/hlt")
    p.add_argument("--cutoff", type=float, default=fm.CUTOFF)
    p.add_argument("--out", default="results/native_stability.csv")
    a = p.parse_args()
    if a.demo:
        return demo()

    rows = rows_for(a.boltz, lambda tag: native_recovery(a.boltz, tag, a.hlt, a.cutoff))
    if not rows:
        raise SystemExit("no native outputs under %s" % a.boltz)
    ns = {r["n_samples"] for r in rows}
    if len(ns) != 1:
        # unequal n would make `min` a different test per target
        raise SystemExit("sample count differs between targets: %s - refold so every "
                         "native has the same n before applying the rule" % sorted(ns))
    missing = [r["target"] for r in rows if r["epitope_recovery"] == ""]
    if missing:
        raise SystemExit("no epitope recovery for %s - check %s and the HLT files"
                         % (missing, a.boltz))

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fields = ["target", "n_samples", "model_0", "iptm_min", "iptm_median", "iptm_max",
              "spread", "n_epitope", "epitope_recovery", "iptm_ok", "recovery_ok",
              "stable"]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print("inclusion rule: min iptm over %d samples >= %.2f  AND  native epitope "
          "recovery >= %.2f\n" % (rows[0]["n_samples"], FLOOR_IPTM, FLOOR_RECOVERY))
    print("%-16s %8s %8s %8s %7s %9s  %-6s %-6s %s"
          % ("target", "min", "median", "spread", "n_epi", "recovery", "iptm", "recov",
             "INCLUDED"))
    for r in rows:
        print("%-16s %8.3f %8.3f %8.3f %7s %9.3f  %-6s %-6s %s"
              % (r["target"], r["iptm_min"], r["iptm_median"], r["spread"],
                 r["n_epitope"], r["epitope_recovery"],
                 "ok" if r["iptm_ok"] else "FAIL", "ok" if r["recovery_ok"] else "FAIL",
                 "yes" if r["stable"] else "no"))
    keep = [r["target"] for r in rows if r["stable"]]
    print("\n%d of %d targets included: %s" % (len(keep), len(rows), " ".join(keep)))
    print("excluded by iptm only:     %s"
          % " ".join(r["target"] for r in rows if not r["iptm_ok"] and r["recovery_ok"]))
    print("excluded by recovery only: %s"
          % " ".join(r["target"] for r in rows if r["iptm_ok"] and not r["recovery_ok"]))
    print("excluded by both:          %s"
          % " ".join(r["target"] for r in rows
                     if not r["iptm_ok"] and not r["recovery_ok"]))
    print("wrote %s" % a.out)


def demo():
    import tempfile
    # the rule has two independent conditions, so the self-check needs all four
    # corners, including the two the single-condition rules get wrong
    cases = {
        "pdb_good_H_T":      ([0.91, 0.88, 0.93, 0.90, 0.92], (20, 0.90)),  # include
        "pdb_cointoss_H_T":  ([0.72, 0.24, 0.41, 0.55, 0.30], (20, 0.90)),  # iptm fails
        "pdb_wrongplace_H_T": ([0.91, 0.90, 0.93, 0.90, 0.92], (18, 0.06)),  # 9g3e case
        "pdb_bothbad_H_T":   ([0.45, 0.18, 0.30, 0.22, 0.28], (20, 0.10)),  # both fail
    }
    with tempfile.TemporaryDirectory() as td:
        for tag, (vals, _) in cases.items():
            d = os.path.join(td, tag, "native", "boltz_results_native")
            os.makedirs(d)
            for i, x in enumerate(vals):
                json.dump({"iptm": x},
                          open(os.path.join(d, "confidence_native_model_%d.json" % i), "w"))
        rows = {r["target"]: r for r in rows_for(td, lambda t: cases[t][1])}

        assert rows["pdb_good"]["stable"] == 1, rows["pdb_good"]
        # a reliable fold in the wrong place must be EXCLUDED: this is the whole
        # point of the amendment, and the condition the iptm floor alone misses
        assert rows["pdb_wrongplace"]["iptm_ok"] == 1
        assert rows["pdb_wrongplace"]["recovery_ok"] == 0
        assert rows["pdb_wrongplace"]["stable"] == 0, rows["pdb_wrongplace"]
        # and a well-placed but unreliable fold is excluded too
        assert rows["pdb_cointoss"]["recovery_ok"] == 1
        assert rows["pdb_cointoss"]["iptm_ok"] == 0
        assert rows["pdb_cointoss"]["stable"] == 0
        assert rows["pdb_bothbad"]["stable"] == 0
        assert rows["pdb_good"]["iptm_median"] == 0.91
        assert abs(rows["pdb_cointoss"]["spread"] - 0.48) < 1e-9
        # included targets sort first
        assert [r["target"] for r in rows_for(td, lambda t: cases[t][1])][0] == "pdb_good"
    print("demo ok")


if __name__ == "__main__":
    main()
