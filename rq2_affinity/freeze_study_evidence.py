from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

OUT = Path("rq2_affinity/results")
queue = pd.read_csv(
    OUT / "rq2_three_study_verification_queue.csv",
    dtype={"pubmed_id": str}
)

source_values = {
    "2P42": ("native cAb-RN05 / SE3", 23.0, "Yeast-display binding titration", "18445622"),
    "2P43": ("native cAb-RN05 / SE3", 23.0, "Yeast-display binding titration", "18445622"),
    "2P44": ("clone 7 / SE5a", 20.0, "Yeast-display binding titration", "18445622"),
    "2P45": ("clone 22 / SE5b", 116.0, "Yeast-display binding titration", "18445622"),
    "2P46": ("clone 22 / SE5b", 116.0, "Yeast-display binding titration", "18445622"),
    "2P47": ("clone 22 / SE5b", 116.0, "Yeast-display binding titration", "18445622"),
    "2P48": ("clone 22 / SE5b", 116.0, "Yeast-display binding titration", "18445622"),
    "7C8W": ("MR17", 83.7, "BLI", "34330908"),
}

records = []
for _, row in queue.iterrows():
    rec = row.to_dict()
    pdb = rec["pdb_id"]
    rec["source_construct"] = ""
    rec["source_kd_nM"] = np.nan
    rec["source_method"] = ""
    rec["verification_status"] = "PENDING_EXACT_VALUE"

    if pdb in source_values:
        construct, kd, method, pmid = source_values[pdb]
        if str(rec["pubmed_id"]) != pmid:
            raise ValueError("PMID mismatch: " + pdb)
        rec["source_construct"] = construct
        rec["source_kd_nM"] = kd
        rec["source_method"] = method
        rec["verification_status"] = (
            "SOURCE_CONSISTENT"
            if np.isclose(float(rec["andd_kd_nM"]), kd, rtol=0.015)
            else "NUMERIC_CONFLICT"
        )
    records.append(rec)

audit = pd.DataFrame(records)
audit.to_csv(
    OUT / "rq2_source_checked_study_records.csv",
    index=False
)

rnase = audit[
    (audit["antigen"] == "RNase")
    & (audit["verification_status"] == "SOURCE_CONSISTENT")
].copy()

rnase = (
    rnase.sort_values(["vhh_key", "pdb_id"])
    .drop_duplicates("vhh_key")
)

geometry = pd.read_csv(OUT / "rnase_interface_audit.csv")

rnase = rnase.merge(
    geometry[
        ["pdb_id", "vhh_key", "residue_contacts",
         "atom_contacts", "vhh_residues", "ag_residues"]
    ].drop_duplicates(["pdb_id", "vhh_key"]),
    on=["pdb_id", "vhh_key"],
    how="left",
    validate="one_to_one"
)

rnase["source_pKD"] = 9 - np.log10(rnase["source_kd_nM"])
rnase.to_csv(
    OUT / "rnase_source_supported_three_vhh.csv",
    index=False
)

print("\n=== SOURCE AUDIT ===")
print(audit[
    ["antigen", "pdb_id", "source_construct",
     "andd_kd_nM", "source_kd_nM", "verification_status"]
].to_string(index=False))

print("\n=== RNASE: THREE SOURCE-SUPPORTED VHH ===")
print(rnase[
    ["pdb_id", "source_construct", "source_kd_nM",
     "vhh_residues", "ag_residues"]
].to_string(index=False))

print("\n=== DESCRIPTIVE CORRELATIONS (n=3) ===")
for metric in [
    "residue_contacts", "atom_contacts",
    "vhh_residues", "ag_residues"
]:
    rho = spearmanr(
        rnase[metric], rnase["source_pKD"]
    ).statistic
    print(f"{metric:20s} rho={rho:+.3f}")

print("\nWARNING: n=3; descriptive only.")
print("No generalization or significance claim.")
