from pathlib import Path
import subprocess
import sys
import re
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"

groups = pd.read_csv(OUT / "affinity_antigen_groups.csv")
records = pd.read_csv(
    OUT / "affinity_real_records_NEEDS_SOURCE_VALIDATION.csv",
    low_memory=False
)

AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
PDB = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

def valid_seq(x):
    return isinstance(x, str) and bool(AA.fullmatch(x.strip().upper()))

candidates = []

for _, group in groups.iterrows():
    if str(group["antigen_name"]).strip().lower() != "spike glycoprotein":
        continue

    sub = records[
        records["ag_key"].astype(str) == str(group["ag_key"])
    ].copy()

    sub["PDB_ID"] = sub["PDB_ID"].astype(str).str.strip().str.upper()
    sub["kd_M"] = pd.to_numeric(sub["kd_M"], errors="coerce")

    sub = sub[
        sub["Ab/Nano H_Chain AA"].map(valid_seq)
        & sub["Ag_Seq"].map(valid_seq)
        & sub["PDB_ID"].map(lambda x: bool(PDB.fullmatch(x)))
        & np.isfinite(sub["kd_M"])
        & (sub["kd_M"] > 0)
    ].copy()

    sub = sub.drop_duplicates(["PDB_ID", "vhh_key", "kd_M"])

    candidates.append({
        "ag_key": group["ag_key"],
        "unique_vhh": sub["vhh_key"].nunique(),
        "pdb_count": sub["PDB_ID"].nunique(),
        "records": len(sub),
        "data": sub,
    })

if not candidates:
    raise SystemExit("No exact-name spike groups found.")

print("\n=== SPIKE ANTIGEN SEQUENCE GROUPS ===")
for c in sorted(candidates, key=lambda x: -x["unique_vhh"]):
    print(
        c["ag_key"],
        "VHH:", c["unique_vhh"],
        "PDB:", c["pdb_count"],
        "records:", c["records"]
    )

eligible = [c for c in candidates if c["unique_vhh"] >= 5]

if not eligible:
    raise SystemExit(
        "\nNo spike group has >=5 unique VHH with PDB. "
        "This is a coverage limitation, not a mapping failure."
    )

eligible.sort(
    key=lambda c: (c["unique_vhh"], c["pdb_count"]),
    reverse=True
)

best = eligible[0]
ties = [
    c for c in eligible
    if (c["unique_vhh"], c["pdb_count"]) ==
       (best["unique_vhh"], best["pdb_count"])
]

if len(ties) > 1:
    raise SystemExit("Ambiguous best groups; manual review required.")

cohort = best["data"].copy()
cohort["pKD"] = -np.log10(cohort["kd_M"])

cohort.to_csv(OUT / "spike_detailed_audit.csv", index=False)

print("\n=== SELECTED COHORT ===")
print("Antigen sequence key:", best["ag_key"])
print("Unique VHH:", best["unique_vhh"])
print("PDB count:", best["pdb_count"])
print("Source distribution:")
print(cohort["Source"].value_counts().to_string())

# Reuse the original frozen geometry implementation.
for script_type in ["map", "analyze"]:
    source = ROOT / f"{script_type}_lysozyme.py"
    target = ROOT / f"{script_type}_spike.py"

    code = source.read_text()
    code = code.replace("lysozyme", "spike")
    code = code.replace("Lysozyme", "Spike")
    code = code.replace("10 DISTINCT VHH", "UNIQUE VHH")

    target.write_text(code)

    print("\nRUNNING", target.name, flush=True)
    subprocess.run(
        [sys.executable, str(target)],
        cwd=ROOT.parent,
        check=True
    )

print("\n=== SPIKE ANALYSIS COMPLETE ===")
