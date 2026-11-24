#!/usr/bin/env python3
"""B5: paired hotspot-vs-unguided comparison, cluster bootstrap, failure taxonomy.

Two outputs:

  results/b5_paired.csv     one row per target: each arm's median, and the
                            paired difference. This is the unit of analysis.
  results/b5_failures.csv   one row per candidate: geometry measurements and
                            the failure classes it falls into.

THE UNIT OF ANALYSIS IS THE TARGET, NEVER THE CANDIDATE. 84 candidates are not
84 independent observations: four come from each target, two per arm, on the
same antigen and the same framework. Every comparison below summarises within
a target first, then treats the per-target differences as the sample.

CLUSTERING RULE, which B5 owes explicitly. `data/processed/targets.csv` carries
SAbDab's `ab_ag_cluster`, and all 21 targets have distinct values because C's
benchmark construction kept one representative per cluster, chosen by
resolution. So each target is its own antibody/antigen cluster and a bootstrap
over clusters is a bootstrap over targets. The script asserts this rather than
assuming it: if a future benchmark ever puts two targets in one cluster, the
resampling has to group them and the assertion will fail loudly.

SEQUENCE DUPLICATES, also owed: 83 distinct sequences among the 84 candidates.
The one duplicate is on 8qf4, where a hotspot candidate's designed sequence
equals an unguided candidate's, so that target's arm comparison is partly a
self-comparison. Recomputed here rather than quoted.

FAILURE TAXONOMY. Measured on B's Boltz prediction against A's generated
complex, which is the designed intent. All three use CA atoms only.

  framework_rmsd  superimpose on the non-CDR residues, RMSD on those. The
                  framework was meant to be held fixed, so this is how far the
                  prediction drifted from a backbone it was supposed to keep.
  cdr_rmsd        after that same framework superposition, RMSD on the CDR
                  residues. Loops that moved while the framework stayed put.
  pose_rmsd       superimpose on the ANTIGEN, then RMSD of the nanobody. How
                  far the nanobody sits from where it was designed to bind.

Thresholds are round numbers fixed in advance, not tuned, and they only label
failures for description. They never select candidates, so they are not the
kind of threshold non-negotiable #3 governs. State them with any count.

  python3 src/evaluate/b5_analysis.py --demo
  python3 src/evaluate/b5_analysis.py
"""
import argparse
import csv
import collections
import glob
import json
import os
import statistics as sm

import numpy as np

import fold_metrics as fm
import interface_metrics as im

FRAMEWORK_RMSD = 2.0     # A; the framework was supposed to be fixed
CDR_RMSD = 3.0           # A; after framework superposition
POSE_RMSD = 10.0         # A; antigen-superposed nanobody displacement
IPTM_CONFIDENT = 0.7     # the confident-but-displaced class
RECOVERY_MISS = 0.1


def superpose(P, Q):
    """Rotation+translation taking P onto Q, and the RMSD on those points.

    Kabsch. Returns (apply_fn, rmsd) so the transform fitted on the framework
    can be applied to the CDRs, which is the whole point of doing it this way.
    """
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    apply_fn = lambda X: (X - pc) @ R.T + qc
    return apply_fn, rmsd(apply_fn(P), Q)


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def ca_map(path, chain):
    return {r: v[0] for r, v in fm.read_heavy(path, chain, atoms=("CA",))[chain].items()}


def geometry(gen_path, bz_path, cdr_idx, n_ag):
    """framework / cdr / pose RMSD between the Boltz prediction and A's design."""
    g_nb, b_nb = ca_map(gen_path, "H"), ca_map(bz_path, "A")
    shared = sorted(set(g_nb) & set(b_nb))
    cdr_res = {i + 1 for i in cdr_idx}          # cdr_idx is 0-based into the chain
    fw = [r for r in shared if r not in cdr_res]
    cd = [r for r in shared if r in cdr_res]
    if len(fw) < 3 or not cd:
        return None
    G = np.array([g_nb[r] for r in fw]); B = np.array([b_nb[r] for r in fw])
    xf, fw_rmsd = superpose(B, G)
    cdr = rmsd(xf(np.array([b_nb[r] for r in cd])), np.array([g_nb[r] for r in cd]))

    # pose: superimpose on the antigen instead, then look at the nanobody
    g_ag_raw, b_ag = ca_map(gen_path, "T"), ca_map(bz_path, "B")
    off = min(g_ag_raw) - 1
    g_ag = {r - off: v for r, v in g_ag_raw.items()}
    ag = sorted(set(g_ag) & set(b_ag))
    if len(ag) < 3 or len(g_ag) != n_ag:
        return None
    xf2, _ = superpose(np.array([b_ag[r] for r in ag]), np.array([g_ag[r] for r in ag]))
    pose = rmsd(xf2(np.array([b_nb[r] for r in shared])),
                np.array([g_nb[r] for r in shared]))
    return dict(framework_rmsd=round(fw_rmsd, 2), cdr_rmsd=round(cdr, 2),
                pose_rmsd=round(pose, 2), n_framework=len(fw), n_cdr=len(cd))


def classify(row):
    """Failure classes for one candidate. A candidate can be in several."""
    out = []
    if row.get("framework_rmsd") is not None and row["framework_rmsd"] > FRAMEWORK_RMSD:
        out.append("framework_distortion")
    if row.get("cdr_rmsd") is not None and row["cdr_rmsd"] > CDR_RMSD:
        out.append("cdr_refolding")
    if row.get("pose_rmsd") is not None and row["pose_rmsd"] > POSE_RMSD:
        out.append("wrong_pose")
    if row.get("n_clashes") and row["n_clashes"] > 0:
        out.append("clashes")
    if (row.get("iptm") is not None and row.get("epitope_recovery") is not None
            and row["iptm"] >= IPTM_CONFIDENT and row["epitope_recovery"] < RECOVERY_MISS):
        out.append("confident_but_displaced")
    return out


def bootstrap(diffs, n=10000, seed=37):
    """Percentile CI for the mean paired difference, resampling TARGETS."""
    if not diffs:
        return None
    rng = np.random.default_rng(seed)
    d = np.asarray(diffs, dtype=float)
    means = d[rng.integers(0, len(d), size=(n, len(d)))].mean(axis=1)
    return (float(d.mean()), float(np.percentile(means, 2.5)),
            float(np.percentile(means, 97.5)),
            float((means <= 0).mean()))


def sign_test(diffs):
    from math import comb
    pos = sum(1 for x in diffs if x > 0)
    neg = sum(1 for x in diffs if x < 0)
    nt = pos + neg
    if nt == 0:
        return pos, neg, len(diffs) - nt, 1.0
    p = 2 * sum(comb(nt, k) for k in range(max(pos, neg), nt + 1)) / 2 ** nt
    return pos, neg, len(diffs) - nt, min(p, 1.0)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true")
    p.add_argument("--run", default="runs/rfab_benchmark_1ce3344e4822")
    p.add_argument("--boltz", default="runs/boltz_out5")
    p.add_argument("--inputs", default="runs/boltz_in")
    p.add_argument("--metrics", default="results/end_to_end_metrics.csv")
    p.add_argument("--interface", default="results/interface_metrics.csv")
    p.add_argument("--targets", default="data/processed/targets.csv")
    p.add_argument("--out-paired", default="results/b5_paired.csv")
    p.add_argument("--out-failures", default="results/b5_failures.csv")
    a = p.parse_args()
    if a.demo:
        return demo()

    # --- clustering rule: assert one cluster per target ----------------------
    tg = list(csv.DictReader(open(a.targets)))
    clusters = collections.Counter(r["ab_ag_cluster"] for r in tg)
    shared = {c: n for c, n in clusters.items() if n > 1}
    if shared:
        raise SystemExit("targets share an ab_ag_cluster %s - the bootstrap must "
                         "resample those together, not as separate targets" % shared)
    print("clustering rule: %d targets, %d distinct ab_ag_cluster values, so one "
          "cluster per target" % (len(tg), len(clusters)))

    ee = list(csv.DictReader(open(a.metrics)))
    iface = {(r["pdb_id"], r["arm"], r["backbone_index"]): r
             for r in csv.DictReader(open(a.interface)) if r["structure_set"] == "boltz"}

    # --- duplicate count, recomputed --------------------------------------
    seqs = []
    for man in sorted(glob.glob(os.path.join(a.inputs, "*", "manifest.csv"))):
        seqs += [(r["target"], r["arm"], r["nanobody_seq"])
                 for r in csv.DictReader(open(man)) if r["kind"] == "rfab"]
    uniq = len({s for _, _, s in seqs})
    cross = [t for t in {t for t, _, _ in seqs}
             if ({s for tt, ar, s in seqs if tt == t and ar == "hotspot"} &
                 {s for tt, ar, s in seqs if tt == t and ar == "unguided"})]
    print("sequences: %d candidates, %d distinct; cross-arm duplicates on %s"
          % (len(seqs), uniq, [t.replace("_H_T", "") for t in cross] or "none"))

    spec = json.load(open(os.path.join(a.run, "inputs", "targets.json")))
    spec = {d["pdb_id"]: d for d in (spec.values() if isinstance(spec, dict) else spec)}

    # --- per candidate: geometry + taxonomy -------------------------------
    fails = []
    for r in ee:
        if r["kind"] != "rfab":
            continue
        pid = r["target"].replace("_H_T", "")
        arm = r["arm"]
        idx = r["candidate_id"].rsplit("_nb_", 1)[1].split("_")[0]
        man = next((m for m in csv.DictReader(
            open(os.path.join(a.inputs, r["target"], "manifest.csv")))
            if m["job"] == r["job"]), None)
        if man is None:
            continue
        cdr_idx = [int(i) for i in man["cdr_idx"].split(";") if i]
        gen, _ = im.candidate_paths(a.run, pid, arm, int(man["backbone"]))
        bz = im.boltz_path(a.boltz, r["target"], r["job"])
        row = dict(target=pid, arm=arm, backbone_index=man["backbone"],
                   candidate_id=r["candidate_id"],
                   native_included=r["native_included"],
                   iptm=float(r["iptm"]) if r["iptm"] else None,
                   epitope_recovery=float(r["epitope_recovery"]) if r["epitope_recovery"] else None,
                   plddt_framework=r["plddt_framework"], plddt_cdr=r["plddt_cdr"])
        ik = (pid, arm, man["backbone"])
        row["n_clashes"] = int(iface[ik]["n_clashes"]) if ik in iface else None
        row["hs_recovery_bb"] = (float(iface[ik]["hs_recovery_bb"])
                                 if ik in iface and iface[ik]["hs_recovery_bb"] else None)
        if gen and bz:
            g = geometry(gen, bz, cdr_idx, int(man["len_antigen"]))
            if g:
                row.update(g)
        row["failures"] = ";".join(classify(row))
        fails.append(row)

    ff = ["target", "arm", "backbone_index", "candidate_id", "native_included",
          "iptm", "epitope_recovery", "hs_recovery_bb", "n_clashes",
          "plddt_framework", "plddt_cdr", "n_framework", "n_cdr",
          "framework_rmsd", "cdr_rmsd", "pose_rmsd", "failures"]
    os.makedirs(os.path.dirname(a.out_failures) or ".", exist_ok=True)
    with open(a.out_failures, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ff, extrasaction="ignore")
        w.writeheader(); w.writerows(fails)

    # --- paired per-target comparison -------------------------------------
    paired = []
    for t in sorted({r["target"] for r in fails}):
        g = [r for r in fails if r["target"] == t]
        row = dict(target=t, native_included=g[0]["native_included"])
        for metric in ("epitope_recovery", "hs_recovery_bb"):
            for arm in ("hotspot", "unguided"):
                v = [r[metric] for r in g if r["arm"] == arm and r[metric] is not None]
                row["%s_%s" % (metric, arm)] = round(sm.median(v), 4) if v else ""
            h, u = row["%s_hotspot" % metric], row["%s_unguided" % metric]
            row["%s_diff" % metric] = round(h - u, 4) if h != "" and u != "" else ""
        row["n_candidates"] = len(g)
        paired.append(row)
    pf = ["target", "native_included", "n_candidates",
          "epitope_recovery_hotspot", "epitope_recovery_unguided", "epitope_recovery_diff",
          "hs_recovery_bb_hotspot", "hs_recovery_bb_unguided", "hs_recovery_bb_diff"]
    with open(a.out_paired, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pf, extrasaction="ignore")
        w.writeheader(); w.writerows(paired)

    # --- report ------------------------------------------------------------
    print("\nwrote %s (%d targets) and %s (%d candidates)"
          % (a.out_paired, len(paired), a.out_failures, len(fails)))
    for scope, keep in (("8 included targets", lambda r: r["native_included"] == "1"),
                        ("all 21 targets", lambda r: True)):
        print("\n=== %s ===" % scope)
        for metric in ("epitope_recovery", "hs_recovery_bb"):
            d = [r["%s_diff" % metric] for r in paired
                 if keep(r) and r["%s_diff" % metric] != ""]
            if not d:
                continue
            pos, neg, tie, pv = sign_test(d)
            bs = bootstrap(d)
            print("%-18s n=%2d pairs  hotspot better %d, unguided %d, tied %d  sign p=%.4f"
                  % (metric, len(d), pos, neg, tie, pv))
            print("%-18s mean paired diff %+.4f  95%% CI [%+.4f, %+.4f]  "
                  "bootstrap P(diff<=0)=%.4f" % ("", bs[0], bs[1], bs[2], bs[3]))
            if neg == 0 and tie:
                # With no negative difference, a resampled mean can only fall at or
                # below zero when every draw happens to be a tie. The "p-value" is
                # then just (tie/n)^n and says nothing about the effect. Quote the
                # sign test instead.
                print("%-18s WARNING: no negative differences, so the bootstrap is "
                      "degenerate here:" % "")
                print("%-18s          P(diff<=0) == (%d/%d)^%d = %.6f, the chance of "
                      "drawing only ties." % ("", tie, len(d), len(d),
                                              (tie / len(d)) ** len(d)))
                print("%-18s          Report the sign test (p=%.4f); do not quote this "
                      "CI as significance." % ("", pv))

    print("\n=== failure taxonomy (Boltz prediction vs A's design; CA atoms) ===")
    print("thresholds: framework>%.1fA, cdr>%.1fA (after framework superposition), "
          "pose>%.1fA, clashes>0, iptm>=%.1f with recovery<%.1f"
          % (FRAMEWORK_RMSD, CDR_RMSD, POSE_RMSD, IPTM_CONFIDENT, RECOVERY_MISS))
    for scope, keep in (("8 included", lambda r: r["native_included"] == "1"),
                        ("all 84", lambda r: True)):
        g = [r for r in fails if keep(r)]
        c = collections.Counter()
        for r in g:
            for lab in (r["failures"].split(";") if r["failures"] else ["none"]):
                c[lab] += 1
        print("\n  %s (n=%d candidates)" % (scope, len(g)))
        for lab, n in c.most_common():
            print("    %-24s %3d  (%.0f%%)" % (lab, n, 100 * n / len(g)))
        for k in ("framework_rmsd", "cdr_rmsd", "pose_rmsd"):
            v = [r[k] for r in g if r.get(k) is not None]
            if v:
                print("    median %-18s %6.2f A   (max %.1f)" % (k, sm.median(v), max(v)))


def demo():
    # superposition must be exact on a rigid transform and recover a known shift
    rng = np.random.default_rng(0)
    P = rng.normal(size=(12, 3))
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    Q = P @ R.T + np.array([5.0, -2.0, 1.0])
    xf, r0 = superpose(P, Q)
    assert r0 < 1e-9, r0
    assert np.allclose(xf(P), Q, atol=1e-9)
    # a pure translation of one point set is fully removed by superposition
    assert superpose(P, P + 3.0)[1] < 1e-9
    # and a genuine deformation is not
    Pd = P.copy(); Pd[0] += 10.0
    assert superpose(Pd, P)[1] > 0.5

    assert abs(rmsd(np.zeros((2, 3)), np.array([[3.0, 0, 0], [0, 4.0, 0]]))
               - np.sqrt((9 + 16) / 2)) < 1e-9

    # taxonomy: each class fires on its own and they can co-occur
    assert classify(dict(framework_rmsd=5.0)) == ["framework_distortion"]
    assert classify(dict(cdr_rmsd=9.0)) == ["cdr_refolding"]
    assert classify(dict(pose_rmsd=30.0)) == ["wrong_pose"]
    assert classify(dict(n_clashes=2)) == ["clashes"]
    assert classify(dict(iptm=0.9, epitope_recovery=0.0)) == ["confident_but_displaced"]
    # confident AND placed is not a failure
    assert classify(dict(iptm=0.9, epitope_recovery=0.8)) == []
    # below every threshold is clean
    assert classify(dict(framework_rmsd=0.5, cdr_rmsd=1.0, pose_rmsd=2.0,
                         n_clashes=0, iptm=0.9, epitope_recovery=0.9)) == []
    assert set(classify(dict(framework_rmsd=5.0, pose_rmsd=30.0, n_clashes=1))) == \
        {"framework_distortion", "wrong_pose", "clashes"}
    # missing measurements must not invent a failure
    assert classify(dict(framework_rmsd=None, cdr_rmsd=None, pose_rmsd=None,
                         n_clashes=None, iptm=None, epitope_recovery=None)) == []

    # bootstrap: deterministic under a seed, and CI brackets the mean
    d = [0.1, 0.2, 0.0, 0.3, -0.05, 0.15, 0.2, 0.0]
    m, lo, hi, pneg = bootstrap(d, n=2000, seed=1)
    assert abs(m - sm.mean(d)) < 1e-12
    assert lo <= m <= hi, (lo, m, hi)
    assert bootstrap(d, n=2000, seed=1) == bootstrap(d, n=2000, seed=1)
    # an all-positive sample should rarely bootstrap to <=0
    assert bootstrap([0.2] * 8, n=2000, seed=1)[3] == 0.0

    pos, neg, tie, pv = sign_test([1, 1, 1, 1, 1, 0, 0, 0])
    assert (pos, neg, tie) == (5, 0, 3), (pos, neg, tie)
    assert abs(pv - 2 / 32) < 1e-12, pv
    print("demo ok")


if __name__ == "__main__":
    main()
