from pathlib import Path
import pandas as pd
import numpy as np

OUT = Path("rq2_affinity/results")

data = pd.read_csv(
    OUT / "rq2_affinity_adjudication.csv"
)

refs = pd.read_csv(
    OUT / "rq2_pdb_publication_audit.csv",
    dtype={"pubmed_id": "string"}
)

refs = refs[
    ["pdb_id", "pubmed_id", "citation_title"]
].drop_duplicates("pdb_id")

data = data.merge(
    refs, on="pdb_id",
    how="left", validate="many_to_one"
)

# One study can contain multiple PDB structures for
# the same VHH. Keep them visible, not independent.
summary = (
    data.groupby(
        ["antigen_cohort", "pubmed_id"],
        dropna=False
    )
    .agg(
        pdbs=("pdb_id", "nunique"),
        records=("pdb_id", "size"),
        unique_vhh=("vhh_key", "nunique"),
        kd_min_nM=("andd_kd_nM", "min"),
        kd_max_nM=("andd_kd_nM", "max"),
        conflicts=(
            "review_status",
            lambda x: int(
                (x == "NUMERIC_CONFLICT_NEEDS_ADJUDICATION").sum()
            )
        )
    )
    .reset_index()
)

summary["kd_fold_range"] = (
    summary["kd_max_nM"] / summary["kd_min_nM"]
)

summary["ranking_candidate_ge3"] = (
    (summary["unique_vhh"] >= 3)
    & (summary["kd_fold_range"] > 1)
)

summary.to_csv(
    OUT / "rq2_study_level_coverage.csv",
    index=False
)

details = data[[
    "antigen_cohort", "pubmed_id", "pdb_id",
    "vhh_key", "andd_kd_nM", "paper_kd_nM",
    "paper_method", "review_status", "citation_title"
]].sort_values(
    ["antigen_cohort", "pubmed_id", "pdb_id"]
)

details.to_csv(
    OUT / "rq2_study_level_records.csv",
    index=False
)

print("\n=== STUDY-LEVEL COVERAGE ===")
print(summary.to_string(
    index=False, float_format="%.4g"
))

print("\n=== STUDIES WITH >=3 DISTINCT VHH ===")
eligible = summary[summary["ranking_candidate_ge3"]]
print(eligible.to_string(
    index=False, float_format="%.4g"
))

print("\nSaved:")
print(OUT / "rq2_study_level_coverage.csv")
print(OUT / "rq2_study_level_records.csv")

print("\nWARNING: KD values still require source verification.")
