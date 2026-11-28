from pathlib import Path
import subprocess
import sys
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
PYTHON = sys.executable

COHORTS = [
    ("rnase", "ribonuclease pancreatic"),
    ("spike", "spike glycoprotein"),
]

groups = pd.read_csv(RESULTS / "affinity_antigen_groups.csv")
records = pd.read_csv(
    RESULTS / "affinity_real_records_NEEDS_SOURCE_VALIDATION.csv",
    low_memory=False
)

map_template = (ROOT / "map_lysozyme.py").read_text()
analysis_template = (ROOT / "analyze_lysozyme.py").read_text()

AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
PDB = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

def valid_sequence(value):
    return (
        isinstance(value, str)
        and bool(AA.fullmatch(value.strip().upper()))
    )

def run(path):
    print(f"\nRUNNING {path.name}", flush=True)
    subprocess.run(
        [PYTHON, str(path)],
        cwd=ROOT.parent,
        check=True
    )

summaries = []

for slug, antigen in COHORTS:
    print(f"\n{'=' * 60}")
    print(f"COHORT: {antigen}")
    print(f"{'=' * 60}", flush=True)

    matches = groups[
        groups["antigen_name"].astype(str).str.lower().eq(antigen)
    ]

    if len(matches) != 1:
        print(
            "SKIP: expected exactly one antigen sequence group;",
            "found", len(matches)
        )
        continue

    ag_key = matches.iloc[0]["ag_key"]

    cohort = records[
        records["ag_key"].astype(str) == str(ag_key)
    ].copy()

    cohort = cohort[
        cohort["Ab/Nano H_Chain AA"].map(valid_sequence)
        & cohort["Ag_Seq"].map(valid_sequence)
    ].copy()

    cohort["PDB_ID"] = (
        cohort["PDB_ID"].astype(str).str.strip().str.upper()
    )
    cohort = cohort[
        cohort["PDB_ID"].map(
            lambda x: bool(PDB.fullmatch(x))
        )
    ].copy()

    cohort["kd_M"] = pd.to_numeric(
        cohort["kd_M"], errors="coerce"
    )
    cohort = cohort[
        np.isfinite(cohort["kd_M"])
        & (cohort["kd_M"] > 0)
    ].copy()

    cohort["pKD"] = -np.log10(cohort["kd_M"])

    cohort = cohort.drop_duplicates([
        "PDB_ID", "vhh_key", "kd_M"
    ])

    data_path = RESULTS / f"{slug}_detailed_audit.csv"
    cohort.to_csv(data_path, index=False)

    n_vhh = cohort["vhh_key"].nunique()
    n_pdb = cohort["PDB_ID"].nunique()

    print("Rows:", len(cohort))
    print("Unique VHH:", n_vhh)
    print("Unique PDB:", n_pdb)

    if n_vhh < 5:
        print("SKIP: fewer than 5 VHH with PDB records")
        continue

    # Reuse the frozen Lysozyme structural analysis.
    # Only cohort paths and display names change.
    mapper = map_template.replace(
        "lysozyme", slug
    ).replace(
        "Lysozyme", antigen
    )

    analyzer = analysis_template.replace(
        "lysozyme", slug
    ).replace(
        "Lysozyme", antigen
    )

    map_path = ROOT / f"map_{slug}.py"
    analyze_path = ROOT / f"analyze_{slug}.py"

    map_path.write_text(mapper)
    analyze_path.write_text(analyzer)

    run(map_path)

    mapped = pd.read_csv(
        RESULTS / f"{slug}_interface_audit.csv"
    )

    successes = mapped[
        mapped["status"] == "MAPPED_CONTACT"
    ]

    mapped_vhh = successes["vhh_key"].nunique()

    print(f"\nMapped unique VHH: {mapped_vhh}/{n_vhh}")

    if mapped_vhh < 5:
        print("SKIP correlation: insufficient mapped VHH")
        continue

    run(analyze_path)

    correlations = pd.read_csv(
        RESULTS / f"{slug}_spearman_sensitivity.csv"
    )

    primary = correlations[
        correlations["analysis"] == "primary_first_pdb"
    ].copy()

    print("\n=== PRIMARY SPEARMAN ===")
    print(primary[
        ["metric", "n_unique_vhh",
         "spearman_rho", "p_value_exploratory",
         "loo_rho_min", "loo_rho_max"]
    ].to_string(index=False))

    primary.insert(0, "antigen", antigen)
    summaries.append(primary)

if summaries:
    combined = pd.concat(summaries, ignore_index=True)
    combined.to_csv(
        RESULTS / "rq2_cross_antigen_exploratory.csv",
        index=False
    )

    print("\n=== CROSS-ANTIGEN SUMMARY ===")
    print(combined[
        ["antigen", "metric", "n_unique_vhh",
         "spearman_rho", "p_value_exploratory"]
    ].to_string(index=False))
else:
    print("\nNo additional cohort passed the mapping gate.")

print("\nDONE. Original Lysozyme files were preserved.")
