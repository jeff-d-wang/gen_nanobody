from pathlib import Path
import pandas as pd
import numpy as np
import hashlib

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

df = pd.read_csv(ROOT / "processed/affinity_candidates.csv",
                 low_memory=False)

def txt(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s.lower() in {"nan", "none", "\\", "unknown"} else s

def seq(x):
    return "".join(txt(x).upper().split())

def key(x):
    return hashlib.sha256(x.encode()).hexdigest()[:16] if x else ""

df["flag"] = df["Predicted_or_Not"].map(lambda x: txt(x).lower())
df["source_name"] = df["Source"].map(txt)
df["provenance_name"] = df["Provenance"].map(txt)
df["antigen_seq"] = df["Ag_Seq"].map(seq)
df["vhh_seq"] = df["Ab/Nano H_Chain AA"].map(seq)
df["ag_key"] = df["antigen_seq"].map(key)
df["vhh_key"] = df["vhh_seq"].map(key)
df["pdb_clean"] = df["PDB_ID"].map(txt)

df["kd_M"] = pd.to_numeric(df["kd_M"], errors="coerce")

real = df[
    (df["flag"] == "real") &
    np.isfinite(df["kd_M"]) &
    (df["kd_M"] > 0)
].copy()

real["pKD"] = -np.log10(real["kd_M"])
real["pair_key"] = real["ag_key"] + ":" + real["vhh_key"]

real["complete_pair"] = (
    real["antigen_seq"].str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+") &
    real["vhh_seq"].str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+")
).fillna(False)

print("\n=== REAL-LABELED NANOBODY KD AUDIT ===")
print("Rows:", len(real))
print("Complete VHH-antigen sequences:",
      int(real["complete_pair"].sum()))
print("With PDB:", int((real["pdb_clean"] != "").sum()))
print("Unique nonempty PDB:",
      real.loc[real["pdb_clean"] != "", "pdb_clean"].nunique())
print("Unique complete sequence pairs:",
      real.loc[real["complete_pair"], "pair_key"].nunique())

print("\n=== KD SUMMARY ===")
print(real["kd_M"].describe().to_string())

source = (
    real.groupby(["source_name", "provenance_name"], dropna=False)
    .agg(
        rows=("kd_M", "size"),
        unique_vhh=("vhh_key", lambda x: x[x != ""].nunique()),
        antigen_names=("Ag_Name", "nunique"),
        with_pdb=("pdb_clean", lambda x: int((x != "").sum())),
        complete_pairs=("complete_pair", "sum"),
        kd_min=("kd_M", "min"),
        kd_max=("kd_M", "max")
    )
    .reset_index()
    .sort_values("rows", ascending=False)
)
source.to_csv(OUT / "affinity_source_audit.csv", index=False)

print("\n=== TOP SOURCES ===")
print(source.head(20).to_string(index=False, max_colwidth=48))

complete = real[real["complete_pair"]].copy()

# Remove identical exact records for coverage counts, but retain
# differing KD values and source provenance for later adjudication.
dedup = complete.drop_duplicates(
    subset=["ag_key", "vhh_key", "kd_M"]
).copy()

groups = (
    dedup.groupby("ag_key")
    .agg(
        rows=("kd_M", "size"),
        unique_vhh=("vhh_key", "nunique"),
        antigen_name=("Ag_Name", "first"),
        kd_min=("kd_M", "min"),
        kd_max=("kd_M", "max"),
        pkd_min=("pKD", "min"),
        pkd_max=("pKD", "max"),
        pdb_count=("pdb_clean", lambda x: x[x != ""].nunique()),
        source_count=("source_name", "nunique")
    )
    .reset_index()
)

groups["pKD_span"] = groups["pkd_max"] - groups["pkd_min"]
groups["candidate_ge5"] = (
    (groups["unique_vhh"] >= 5) &
    (groups["pKD_span"] > 0)
)

groups = groups.sort_values(
    ["candidate_ge5", "unique_vhh", "pKD_span"],
    ascending=[False, False, False]
)

groups.to_csv(OUT / "affinity_antigen_groups.csv", index=False)

print("\n=== TOP 20 ANTIGEN GROUPS ===")
print(groups.head(20).to_string(
    index=False, max_colwidth=45
))

print("\n=== GROUP COUNTS ===")
print("Groups:", len(groups))
print("Groups with >= 2 VHH:",
      int((groups["unique_vhh"] >= 2).sum()))
print("Groups with >= 5 VHH:",
      int((groups["unique_vhh"] >= 5).sum()))
print("Candidates >= 5 VHH and variable KD:",
      int(groups["candidate_ge5"].sum()))

print("\n=== MISSING ANTIGEN SEQUENCES ===")
print(int((real["antigen_seq"] == "").sum()))

print("\n=== NONEMPTY AFFINITY METHODS ===")
print(real["Affinity_Method"]
      .fillna("<missing>").value_counts().head(20).to_string())

# Preserve every real-labeled row for literature-level audit.
real.to_csv(
    OUT / "affinity_real_records_NEEDS_SOURCE_VALIDATION.csv",
    index=False
)

print("\nSaved reports in:", OUT)
print("NOTE: real label is not independent experimental validation.")
