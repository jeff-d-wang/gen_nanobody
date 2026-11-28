from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

df = pd.read_csv(OUT / "lysozyme_interface_audit.csv")
df = df[df["status"] == "MAPPED_CONTACT"].copy()

metrics = [
    "residue_contacts",
    "atom_contacts",
    "vhh_residues",
    "ag_residues",
]

required = ["pdb_id", "vhh_key", "kd_M", "pKD"] + metrics
if any(c not in df.columns for c in required):
    raise ValueError("Missing required interface fields")

df = df.drop_duplicates(["pdb_id", "vhh_key", "kd_M"])
for c in ["pKD"] + metrics:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["pKD"] + metrics)

# Primary analysis:
# One PDB/measurement per VHH, selected by a fixed,
# affinity-independent ordering, not by score.
primary = (
    df.sort_values(["vhh_key", "pdb_id"])
      .drop_duplicates("vhh_key", keep="first")
      .copy()
)

# Sensitivity analysis:
# Median affinity and structural descriptors across
# available records for each unique VHH.
median_data = (
    df.groupby("vhh_key", as_index=False)
      [["pKD"] + metrics].median()
)

# Additional sensitivity:
# For VHHs represented by multiple PDBs, choose the
# lexicographically last PDB instead.
alternative = (
    df.sort_values(["vhh_key", "pdb_id"])
      .drop_duplicates("vhh_key", keep="last")
      .copy()
)

def correlation(data, metric):
    x = data[metric].to_numpy(dtype=float)
    y = data["pKD"].to_numpy(dtype=float)

    if len(x) < 4 or len(np.unique(x)) < 2:
        return np.nan, np.nan

    result = spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)

rows = []
for analysis_name, data in [
    ("primary_first_pdb", primary),
    ("median_per_vhh", median_data),
    ("alternative_last_pdb", alternative),
]:
    for metric in metrics:
        rho, p = correlation(data, metric)

        loo = []
        for i in range(len(data)):
            part = data.drop(data.index[i])
            r, _ = correlation(part, metric)
            if np.isfinite(r):
                loo.append(r)

        rows.append({
            "analysis": analysis_name,
            "metric": metric,
            "n_unique_vhh": len(data),
            "spearman_rho": rho,
            "p_value_exploratory": p,
            "loo_rho_min": min(loo) if loo else np.nan,
            "loo_rho_max": max(loo) if loo else np.nan,
        })

results = pd.DataFrame(rows)
results.to_csv(
    OUT / "lysozyme_spearman_sensitivity.csv",
    index=False
)

primary.to_csv(
    OUT / "lysozyme_primary_10_vhh.csv",
    index=False
)

median_data.to_csv(
    OUT / "lysozyme_median_per_vhh.csv",
    index=False
)

print("\n=== PRIMARY: 10 DISTINCT VHH ===")
print(primary[
    ["pdb_id", "vhh_key", "pKD"] + metrics
].sort_values("pKD", ascending=False).to_string(index=False))

print("\n=== SPEARMAN AND SENSITIVITY RESULTS ===")
print(results.to_string(index=False, float_format="%.4f"))

# Visualize the primary analysis.
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, metric in zip(axes.flat, metrics):
    rho, p = correlation(primary, metric)
    ax.scatter(primary[metric], primary["pKD"], s=50)
    ax.set_xlabel(metric.replace("_", " "))
    ax.set_ylabel("Reported pKD")
    ax.set_title(f"Spearman rho={rho:.3f}, p={p:.3g}")
    ax.grid(alpha=0.2)

fig.suptitle(
    "Lysozyme: experimental-complex geometry vs. reported affinity\n"
    "Exploratory; one PDB per unique VHH (n=10)",
    fontsize=13
)
fig.tight_layout()

figure = OUT / "lysozyme_affinity_correlations.png"
fig.savefig(figure, dpi=180)
plt.close(fig)

summary = {
    "primary_unique_vhh": len(primary),
    "source_structure_records": len(df),
    "metrics": metrics,
    "primary_selection": (
        "Lexicographically first PDB per unique VHH; "
        "does not optimize the measured correlation"
    ),
    "sensitivity": [
        "Median per VHH",
        "Lexicographically last PDB per VHH",
        "Leave-one-VHH-out"
    ],
    "interpretation_limitations": [
        "Exploratory, small cohort",
        "Assay comparability not established",
        "Experimental structures, not sequence-only prediction",
        "Unadjusted p-values across multiple metrics",
        "PDB selection may affect associations",
        "Geometry alone cannot establish binding free energy"
    ]
}

(OUT / "lysozyme_analysis_summary.json").write_text(
    json.dumps(summary, indent=2)
)

print("\nSaved:")
print(OUT / "lysozyme_spearman_sensitivity.csv")
print(OUT / "lysozyme_primary_10_vhh.csv")
print(OUT / "lysozyme_median_per_vhh.csv")
print(figure)
