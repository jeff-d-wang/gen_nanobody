from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT = Path("rq2_affinity/results")
df = pd.read_csv(OUT / "rq2_affinity_adjudication.csv")

metrics = [
    "residue_contacts",
    "atom_contacts",
    "vhh_residues",
    "ag_residues",
]

# Freeze the existing analysis definition.
# Exclude 4PGJ in the strict lysozyme cohort.
df = df[
    ~(
        (df.antigen_cohort == "lysozyme") &
        (df.pdb_id == "4PGJ")
    )
].copy()

df["paper_kd_nM"] = pd.to_numeric(
    df["paper_kd_nM"], errors="coerce"
)

def evaluate(data, label, mode, antigen):
    rows = []

    for metric in metrics:
        sub = data.dropna(subset=[metric, label])

        if (
            len(sub) < 4 or
            sub[metric].nunique() < 2 or
            sub[label].nunique() < 2
        ):
            rho = np.nan
            p = np.nan
        else:
            result = spearmanr(sub[metric], sub[label])
            rho = result.statistic
            p = result.pvalue

        rows.append({
            "antigen": antigen,
            "scenario": mode,
            "metric": metric,
            "n": len(sub),
            "rho": rho,
            "p_exploratory": p
        })

    return rows

results = []

for antigen, group in df.groupby("antigen_cohort"):
    # One observation per VHH, same first-PDB convention.
    primary = (
        group.sort_values(["vhh_key", "pdb_id"])
             .drop_duplicates("vhh_key")
             .copy()
    )

    primary["label_original"] = primary["pKD"]

    # Hypothetical substitution for documented conflicts.
    # Not a validated corrected dataset.
    primary["label_paper_candidate"] = np.where(
        primary["paper_kd_nM"].notna() &
        (primary["paper_kd_nM"] > 0),
        9 - np.log10(primary["paper_kd_nM"]),
        primary["pKD"]
    )

    # Alternative: remove records with known conflicts.
    clean_subset = primary[
        primary.review_status !=
        "NUMERIC_CONFLICT_NEEDS_ADJUDICATION"
    ].copy()

    results += evaluate(
        primary, "label_original",
        "ANDD_original", antigen
    )
    results += evaluate(
        primary, "label_paper_candidate",
        "paper_value_scenario_UNVERIFIED", antigen
    )
    results += evaluate(
        clean_subset, "label_original",
        "exclude_known_conflicts", antigen
    )

result = pd.DataFrame(results)
path = OUT / "rq2_label_sensitivity.csv"
result.to_csv(path, index=False)

print("\n=== LABEL SENSITIVITY ===")
print(result.to_string(
    index=False,
    float_format=lambda x: f"{x:.4f}"
))

print("\nSaved:", path)
print(
    "\nWARNING: Paper-value substitution is hypothetical; "
    "none of these scenarios is a fully verified benchmark."
)
