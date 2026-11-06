#!/usr/bin/env python3
"""Figure 3: ProteinMPNN CDR sequence recovery, per CDR, across targets.

One dot per TARGET (mean over that target's sampled sequences), not one per
sequence - designs from a single target are not independent observations
(proposal, statistical treatment). Box summarises the across-target spread.

  python3 src/evaluate/fig3_recovery.py [results/cdr_recovery.csv]
"""
import csv
import random
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else "results/cdr_recovery.csv"
OUT = "results/fig3_cdr_recovery"
GROUPS = {"cdr1": "CDR1", "cdr2": "CDR2", "cdr3": "CDR3", "all": "pooled"}
BLUE, INK, MUTED, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#dcdbd6"
CHANCE = 1 / 20  # uniform over the 20 amino acids

per_target = defaultdict(lambda: defaultdict(list))
with open(SRC, newline="") as f:
    for r in csv.DictReader(f):
        per_target[r["cdr"]][r["target"]].append(float(r["recovery"]))
if not per_target:
    raise SystemExit("no rows in %s - run src/evaluate/run_mpnn_cdr.sh first" % SRC)

data = [[sum(v) / len(v) for v in per_target[g].values()] for g in GROUPS]
n_targets = max(len(per_target[g]) for g in GROUPS)
rng = random.Random(0)

fig, ax = plt.subplots(figsize=(5.2, 3.4))
bp = ax.boxplot(data, widths=0.5, showfliers=False, medianprops=dict(color=INK, lw=1.6),
                whiskerprops=dict(color=MUTED, lw=1), capprops=dict(color=MUTED, lw=1),
                boxprops=dict(color=MUTED, lw=1), patch_artist=True)
for (g, box), vals in zip(zip(GROUPS, bp["boxes"]), data):
    pooled = g == "all"
    box.set_facecolor("#f0efec" if pooled else "#eef4fc")
    i = list(GROUPS).index(g) + 1
    ax.scatter([i + rng.uniform(-0.13, 0.13) for _ in vals], vals, s=22, zorder=3,
               color=MUTED if pooled else BLUE, alpha=0.85, linewidths=0.8,
               edgecolors="white")
    ax.text(i, 1.04, "%.2f" % (sum(vals) / len(vals)), ha="center", va="bottom",
            fontsize=9, color=INK)

ax.axhline(CHANCE, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
ax.text(0.55, CHANCE + 0.02, "chance (1/20)", fontsize=8, color=MUTED, va="bottom", ha="left")

ax.set_xticks(range(1, len(GROUPS) + 1), list(GROUPS.values()))
ax.set_ylim(0, 1.12)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_ylabel("sequence recovery vs native", fontsize=10, color=INK)
ax.set_title("CDR sequence recovery (n=%d targets, mean per target)" % n_targets,
             fontsize=10.5, color=INK, loc="left", pad=14)
ax.yaxis.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(GRID)
ax.tick_params(colors=MUTED, length=0, labelsize=9)

fig.tight_layout()
for ext in ("pdf", "png"):          # pdf for LaTeX, png to eyeball
    fig.savefig("%s.%s" % (OUT, ext), dpi=200, bbox_inches="tight")
print("wrote %s.pdf / .png  (%d targets)" % (OUT, n_targets))
for g, vals in zip(GROUPS, data):
    print("  %-6s mean %.3f  min %.3f  max %.3f" % (g, sum(vals) / len(vals), min(vals), max(vals)))
