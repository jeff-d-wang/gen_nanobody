from pathlib import Path
import pandas as pd
import numpy as np
import json

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"

cohorts = ["lysozyme", "rnase", "spike"]
frames = []

for cohort in cohorts:
    df = pd.read_csv(
        OUT / f"{cohort}_interface_audit.csv"
    )
    df = df[df["status"] == "MAPPED_CONTACT"].copy()
    df["antigen_cohort"] = cohort
    frames.append(df)

data = pd.concat(frames, ignore_index=True)
data["pdb_id"] = data["pdb_id"].str.upper()

# Source-backed values from the primary publications.
# Not automatic replacements: these require construct/assay
# reconciliation before use in a final statistical analysis.
paper_values = {
    "6JB2": {
        "construct": "D3-L11 Y102A",
        "paper_kd_nM": 295.0,
        "method": "SPR kinetics",
        "doi": "10.1038/s41598-019-50722-y",
        "evidence": "Table 2"
    },
    "6JB8": {
        "construct": "D3-L11 WT",
        "paper_kd_nM": 0.129,
        "method": "SPR kinetics",
        "doi": "10.1038/s41598-019-50722-y",
        "evidence": "Table 2"
    },
    "7F5G": {
        "construct": "DL4",
        "paper_kd_nM": 0.245,
        "method": "BLI",
        "doi": "10.1016/j.ijbiomac.2022.04.096",
        "evidence": "DL4-RBD binding kinetics"
    },
    "7F5H": {
        "construct": "DL28",
        "paper_kd_nM": 1.56,
        "method": "BLI",
        "doi": "10.3389/fmicb.2022.875840",
        "evidence": "Figure 1C"
    }
}

data["andd_kd_nM"] = data["kd_M"] * 1e9
data["paper_kd_nM"] = np.nan
data["paper_construct"] = ""
data["paper_method"] = ""
data["paper_doi"] = ""
data["evidence_location"] = ""
data["review_status"] = "PENDING_SOURCE_REVIEW"

for pdb, evidence in paper_values.items():
    mask = data["pdb_id"] == pdb
    data.loc[mask, "paper_kd_nM"] = evidence["paper_kd_nM"]
    data.loc[mask, "paper_construct"] = evidence["construct"]
    data.loc[mask, "paper_method"] = evidence["method"]
    data.loc[mask, "paper_doi"] = evidence["doi"]
    data.loc[mask, "evidence_location"] = evidence["evidence"]
    data.loc[mask, "review_status"] = (
        "NUMERIC_CONFLICT_NEEDS_ADJUDICATION"
    )

data["paper_to_andd_ratio"] = (
    data["paper_kd_nM"] / data["andd_kd_nM"]
)

data["corrected_pKD_candidate"] = np.where(
    data["paper_kd_nM"] > 0,
    9 - np.log10(data["paper_kd_nM"]),
    np.nan
)

# This is NOT a cleaned final analysis dataset.
data.to_csv(
    OUT / "rq2_affinity_adjudication.csv",
    index=False
)

conflicts = data[
    data["review_status"] ==
    "NUMERIC_CONFLICT_NEEDS_ADJUDICATION"
].copy()

conflicts.to_csv(
    OUT / "rq2_affinity_numeric_conflicts.csv",
    index=False
)

summary = {
    "total_structural_observations": len(data),
    "unique_pdbs": int(data.pdb_id.nunique()),
    "numeric_conflict_pdbs": len(paper_values),
    "verified_final_statistical_cohort": False,
    "previous_correlations": (
        "Exploratory ANDD-label associations only; "
        "not validated affinity benchmark results"
    ),
    "next_actions": [
        "Reconcile exact constructs and assay definitions",
        "Trace remaining KD labels to original sources",
        "Freeze independently verified measurements",
        "Recalculate correlations without tuning descriptors"
    ]
}

(OUT / "rq2_provenance_checkpoint.json").write_text(
    json.dumps(summary, indent=2)
)

print("\n=== DOCUMENTED KD DISCREPANCIES ===")
print(conflicts[
    ["antigen_cohort", "pdb_id", "paper_construct",
     "andd_kd_nM", "paper_kd_nM",
     "paper_to_andd_ratio", "paper_method"]
].to_string(index=False, float_format="%.5g"))

print("\n=== REVIEW STATUS ===")
print(data.review_status.value_counts().to_string())

print("\nSaved:")
print(OUT / "rq2_affinity_adjudication.csv")
print(OUT / "rq2_affinity_numeric_conflicts.csv")
print(OUT / "rq2_provenance_checkpoint.json")
print("\nExisting geometry and Spearman outputs unchanged.")
