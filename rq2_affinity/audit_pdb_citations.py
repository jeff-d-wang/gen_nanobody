from pathlib import Path
import json
import time
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"

cohorts = {
    "lysozyme": "lysozyme_interface_audit.csv",
    "rnase": "rnase_interface_audit.csv",
    "spike": "spike_interface_audit.csv",
}

session = requests.Session()
session.headers.update({
    "User-Agent": "RQ2-affinity-provenance-audit/1.0"
})

rows = []
cache = {}

for antigen, filename in cohorts.items():
    df = pd.read_csv(OUT / filename)

    for pdb in sorted(df["pdb_id"].dropna().unique()):
        pdb = str(pdb).upper()

        if pdb not in cache:
            url = (
                "https://data.rcsb.org/rest/v1/"
                f"core/entry/{pdb}"
            )
            try:
                response = session.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()

                citation = data.get(
                    "rcsb_primary_citation", {}
                ) or {}

                cache[pdb] = {
                    "citation_title": citation.get("title"),
                    "pubmed_id": citation.get("pdbx_database_id_PubMed"),
                    "publication_doi": None,
                    "year": citation.get("year"),
                    "journal": citation.get("rcsb_journal_abbrev"),
                    "status": "retrieved",
                }
                print("FETCHED", pdb, flush=True)
            except Exception as exc:
                cache[pdb] = {
                    "status": "error",
                    "error": str(exc)
                }
                print("ERROR", pdb, exc, flush=True)

            time.sleep(0.2)

        rows.append({
            "antigen": antigen,
            "pdb_id": pdb,
            **cache[pdb]
        })

result = pd.DataFrame(rows)

# The PDB primary citation is a STRUCTURE reference.
# It does not automatically validate the reported KD.
result.to_csv(
    OUT / "rq2_pdb_publication_audit.csv",
    index=False
)

print("\n=== PDB PUBLICATION AUDIT ===")
print(result.to_string(
    index=False, max_colwidth=65
))

print("\nSaved:",
      OUT / "rq2_pdb_publication_audit.csv")
