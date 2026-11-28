from pathlib import Path
import pandas as pd
import numpy as np
import json
import hashlib
import re

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw" / "andd_v2"
PROCESSED = ROOT / "processed"
RESULTS = ROOT / "results"

PROCESSED.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

files = list(RAW.glob("*.xlsx"))
if not files:
    raise SystemExit(f"No ANDD Excel found in {RAW}")

print(f"Reading: {files[0]}", flush=True)
df = pd.read_excel(files[0], sheet_name="Sheet3")
print(f"Loaded {len(df):,} rows and {len(df.columns)} columns", flush=True)

def column(prefix):
    matches = [c for c in df.columns if str(c).startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one column starting {prefix!r}; found {matches}")
    return matches[0]

kd_col = column("Affinity_Kd(M)")
seq_col = "Ab/Nano H_Chain AA"

required = [
    "Ab_or_Nano", seq_col, "Ag_Name", "Ag_Seq",
    "Affinity_Method", "Predicted_or_Not",
    "Provenance", "PDB_ID", "Source", "Reason_Code",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

def clean(v):
    if pd.isna(v):
        return ""
    return str(v).strip()

def norm(v):
    return clean(v).lower()

def sha(v):
    return hashlib.sha256(clean(v).encode()).hexdigest()[:16]

print("\n=== DATA FIELD DISTRIBUTIONS ===")
for col in ["Ab_or_Nano", "Predicted_or_Not",
            "Affinity_Method", "Reason_Code"]:
    print(f"\n{col}")
    print(df[col].fillna("<MISSING>").astype(str)
          .value_counts().head(25).to_string())

print("\n=== PROVENANCE EXAMPLES ===")
print(df["Provenance"].fillna("<MISSING>").astype(str)
      .value_counts().head(20).to_string())

# Retain the entire original table through the candidate-selection stage.
kind = df["Ab_or_Nano"].map(norm)
nano_mask = kind.str.contains("nano|vhh|single.domain|sdab", regex=True)

nano = df.loc[nano_mask].copy()

nano["kd_M"] = pd.to_numeric(nano[kd_col], errors="coerce")
nano["valid_kd"] = np.isfinite(nano["kd_M"]) & (nano["kd_M"] > 0)
nano["pKD"] = np.where(nano["valid_kd"], -np.log10(nano["kd_M"]), np.nan)

nano["prediction_flag_raw"] = nano["Predicted_or_Not"].map(norm)
nano["affinity_method_raw"] = nano["Affinity_Method"].map(clean)
nano["provenance_raw"] = nano["Provenance"].map(clean)

# Conservative triage. These are candidates, NOT confirmed experimental labels.
negative_flags = {
    "no", "false", "0", "real", "not predicted", "experimental",
    "experiment", "measured"
}
positive_flags = {
    "yes", "true", "1", "predicted", "prediction",
    "computational", "in silico", "in-silico"
}
missing_methods = {
    "", "nan", "none", "unknown", "not_reported",
    "not reported", "-", "\\"
}

def classify(row):
    if not row["valid_kd"]:
        return "invalid_or_missing_kd"

    flag = row["prediction_flag_raw"]
    method = norm(row["affinity_method_raw"])
    provenance = norm(row["provenance_raw"])

    if flag in positive_flags:
        return "predicted_exclude"

    if flag not in negative_flags:
        return "unknown_prediction_status"

    if any(token in method for token in
           ("predicted", "computational", "in silico", "antipasti")):
        return "predicted_exclude"

    if "antipasti" in provenance:
        return "provenance_requires_review"

    if method in missing_methods:
        return "missing_affinity_method"

    return "experimental_candidate_needs_review"

nano["curation_status"] = nano.apply(classify, axis=1)

nano["antigen_key"] = nano["Ag_Seq"].map(sha)
nano["vhh_key"] = nano[seq_col].map(sha)
nano["method_group"] = nano["Affinity_Method"].map(
    lambda v: norm(v) or "<missing>"
)

# Full nanobody subset, with original ANDD fields preserved.
nano.to_csv(PROCESSED / "affinity_candidates.csv", index=False)

experimental = nano[
    nano["curation_status"] == "experimental_candidate_needs_review"
].copy()

experimental.to_csv(
    PROCESSED / "affinity_experimental_candidates_NEEDS_REVIEW.csv",
    index=False
)

# Group only where antigen sequence and reported method agree.
if len(experimental):
    coverage = (
        experimental.groupby(
            ["antigen_key", "method_group"], dropna=False
        )
        .agg(
            rows=("kd_M", "size"),
            unique_vhhs=("vhh_key", "nunique"),
            antigen_name=("Ag_Name", "first"),
            kd_min_M=("kd_M", "min"),
            kd_max_M=("kd_M", "max"),
            pKD_min=("pKD", "min"),
            pKD_max=("pKD", "max"),
            unique_pdbs=("PDB_ID", lambda x: x.dropna().nunique()),
            sources=("Source", lambda x: "; ".join(
                sorted(set(x.dropna().astype(str)))[:8]
            )),
        )
        .reset_index()
    )
    coverage["pKD_range"] = coverage["pKD_max"] - coverage["pKD_min"]
    coverage["candidate_ranking_group"] = (
        (coverage["unique_vhhs"] >= 5)
        & (coverage["pKD_range"] > 0)
    )
    coverage = coverage.sort_values(
        ["candidate_ranking_group", "unique_vhhs", "pKD_range"],
        ascending=[False, False, False]
    )
else:
    coverage = pd.DataFrame(columns=[
        "antigen_key", "method_group", "rows", "unique_vhhs",
        "candidate_ranking_group"
    ])

coverage.to_csv(RESULTS / "affinity_coverage.csv", index=False)

status_counts = {
    str(k): int(v)
    for k, v in nano["curation_status"].value_counts().items()
}
audit = {
    "source_file": files[0].name,
    "total_andd_rows": int(len(df)),
    "nanobody_rows": int(len(nano)),
    "valid_nanobody_kd_rows": int(nano["valid_kd"].sum()),
    "experimental_candidates_needing_review": int(len(experimental)),
    "status_counts": status_counts,
    "candidate_ranking_groups": int(
        coverage["candidate_ranking_group"].sum()
    ) if len(coverage) else 0,
    "note": (
        "Experimental candidates are NOT confirmed ground truth. "
        "Verify source publications, measurement method, exact constructs, "
        "assay conditions and duplicate records before RQ2 inference."
    ),
}

(RESULTS / "curation_audit.json").write_text(
    json.dumps(audit, indent=2, ensure_ascii=False)
)

print("\n" + "=" * 55)
print("RQ2 — ANDD INITIAL CURATION RESULTS")
print("=" * 55)
print(json.dumps(audit, indent=2, ensure_ascii=False))

if len(coverage):
    print("\n=== TOP 15 CANDIDATE AFFINITY GROUPS ===")
    cols = [
        "antigen_name", "method_group", "rows",
        "unique_vhhs", "pKD_range", "unique_pdbs",
        "candidate_ranking_group"
    ]
    print(coverage[cols].head(15).to_string(index=False))

print("\nOutputs:")
for p in [
    PROCESSED / "affinity_candidates.csv",
    PROCESSED / "affinity_experimental_candidates_NEEDS_REVIEW.csv",
    RESULTS / "affinity_coverage.csv",
    RESULTS / "curation_audit.json"
]:
    print(" ", p)
