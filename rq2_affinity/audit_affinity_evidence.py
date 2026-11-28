from pathlib import Path
from urllib.parse import quote
import requests
import pandas as pd
import xml.etree.ElementTree as ET
import re
import time

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

sources = {
    "lysozyme": "lysozyme_interface_audit.csv",
    "rnase": "rnase_interface_audit.csv",
    "spike": "spike_interface_audit.csv",
}

session = requests.Session()
session.headers["User-Agent"] = "RQ2-research-audit/1.0"

def fetch_json(url, params=None):
    r = session.get(url, params=params, timeout=35)
    r.raise_for_status()
    return r.json()

# 1. Load structural observations and previously retrieved PMIDs.
observations = []
for antigen, filename in sources.items():
    df = pd.read_csv(OUT / filename)
    df = df[df.status == "MAPPED_CONTACT"].copy()
    df["antigen"] = antigen
    observations.append(df)

obs = pd.concat(observations, ignore_index=True)

refs = pd.read_csv(
    OUT / "rq2_pdb_publication_audit.csv",
    dtype={"pubmed_id": "string"}
)

refs = refs[["pdb_id", "pubmed_id", "citation_title"]]
obs = obs.merge(refs, on="pdb_id", how="left", validate="many_to_one")

# 2. Retrieve exact publication DOI from RCSB.
doi_by_pdb = {}
for pdb in sorted(obs.pdb_id.unique()):
    try:
        d = fetch_json(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb}"
        )
        cite = d.get("rcsb_primary_citation") or {}
        doi_by_pdb[pdb] = cite.get("pdbx_database_id_DOI")
        print("RCSB", pdb, doi_by_pdb[pdb], flush=True)
    except Exception as e:
        print("RCSB ERROR", pdb, e, flush=True)
        doi_by_pdb[pdb] = None
    time.sleep(0.15)

obs["publication_doi"] = obs.pdb_id.map(doi_by_pdb)

# 3. Query publication metadata, and accessible full texts.
pmids = sorted(
    {str(x) for x in obs.pubmed_id.dropna()
     if str(x).isdigit()}
)

KEYWORDS = re.compile(
    r"(?i)\b(?:"
    r"K\s*[_\-]?\s*[dD]|dissociation|"
    r"affinity|binding kinetics|surface plasmon|"
    r"\bSPR\b|isothermal titration|"
    r"\bITC\b|biolayer interferometry|"
    r"\bBLI\b|equilibrium constant"
    r")\b"
)

paper_rows = []
evidence_rows = []

for pmid in pmids:
    record = {
        "pubmed_id": pmid,
        "metadata_status": "not_retrieved",
        "fulltext_status": "not_checked",
    }

    try:
        data = fetch_json(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            {"query": f"EXT_ID:{pmid} AND SRC:MED",
             "format": "json", "pageSize": 1}
        )
        found = data.get("resultList", {}).get("result", [])
        if not found:
            raise ValueError("No Europe PMC record")

        p = found[0]
        pmcid = p.get("pmcid")
        record.update({
            "metadata_status": "retrieved",
            "title": p.get("title"),
            "doi": p.get("doi"),
            "pmcid": pmcid,
            "year": p.get("pubYear"),
            "is_open_access": p.get("isOpenAccess"),
        })

        # Abstract is a clue, never proof of exact KD provenance.
        abstract = p.get("abstractText") or ""
        record["abstract_affinity_terms"] = bool(
            KEYWORDS.search(abstract)
        )

        if not pmcid:
            record["fulltext_status"] = "no_pmcid"
        else:
            url = (
                "https://www.ebi.ac.uk/europepmc/webservices/"
                f"rest/{quote(pmcid)}/fullTextXML"
            )
            response = session.get(url, timeout=45)
            if response.status_code == 404:
                record["fulltext_status"] = "not_available"
            else:
                response.raise_for_status()
                xml = ET.fromstring(response.content)

                # Extract paragraphs and table cells separately,
                # preserving local context around affinity terms.
                passages = []
                for elem in xml.iter():
                    if elem.tag.rsplit("}", 1)[-1] in {
                        "p", "td", "th", "caption"
                    }:
                        t = " ".join("".join(
                            elem.itertext()).split()
                        )
                        if t:
                            passages.append(t)

                matches = []
                for passage in passages:
                    hit = KEYWORDS.search(passage)
                    if hit:
                        start = max(0, hit.start() - 65)
                        end = min(len(passage), hit.end() + 115)
                        snippet = passage[start:end]
                        if snippet not in matches:
                            matches.append(snippet)

                record["fulltext_status"] = "retrieved"
                record["affinity_passage_count"] = len(matches)

                for j, snippet in enumerate(matches[:15], 1):
                    evidence_rows.append({
                        "pubmed_id": pmid,
                        "passage_number": j,
                        "snippet": snippet,
                        "evidence_type": "keyword_context_only",
                        "verified_kd": False,
                    })

        print(
            "PAPER", pmid, record["fulltext_status"],
            flush=True
        )
    except Exception as e:
        record["error"] = str(e)
        print("PAPER ERROR", pmid, e, flush=True)

    paper_rows.append(record)
    time.sleep(0.2)

papers = pd.DataFrame(paper_rows)
evidence = pd.DataFrame(
    evidence_rows,
    columns=[
        "pubmed_id", "passage_number", "snippet",
        "evidence_type", "verified_kd"
    ]
)

# 4. Attach paper evidence to each KD observation.
obs["pubmed_id"] = obs["pubmed_id"].astype("string")
papers["pubmed_id"] = papers["pubmed_id"].astype("string")

merged = obs.merge(
    papers,
    on="pubmed_id",
    how="left",
    suffixes=("_pdb", "_paper"),
    validate="many_to_one"
)

merged["affinity_validation"] = "requires_manual_review"

# Expose repeated VHH measurements and cross-paper provenance.
merged["same_vhh_observations"] = (
    merged.groupby(["antigen", "vhh_key"])["pdb_id"]
    .transform("size")
)

cols = [
    "antigen", "pdb_id", "vhh_key", "kd_M", "pKD",
    "pubmed_id", "publication_doi", "title",
    "doi", "year", "fulltext_status",
    "affinity_passage_count",
    "same_vhh_observations", "affinity_validation"
]
cols = [c for c in cols if c in merged.columns]

merged[cols].to_csv(
    OUT / "rq2_affinity_evidence_audit.csv",
    index=False
)

papers.to_csv(
    OUT / "rq2_publication_metadata.csv", index=False
)
evidence.to_csv(
    OUT / "rq2_affinity_text_evidence.csv", index=False
)

print("\n=== PUBLICATION COVERAGE ===")
print("Unique PDB:", obs.pdb_id.nunique())
print("Unique PubMed IDs:", len(pmids))
print(papers.fulltext_status.value_counts().to_string())

print("\n=== PUBLICATION METADATA ===")
print(papers[
    ["pubmed_id", "year", "doi", "fulltext_status",
     "affinity_passage_count"]
    if "affinity_passage_count" in papers.columns
    else ["pubmed_id", "year", "doi", "fulltext_status"]
].to_string(index=False))

print("\n=== OUTPUTS ===")
print("rq2_affinity_evidence_audit.csv")
print("rq2_publication_metadata.csv")
print("rq2_affinity_text_evidence.csv")
print("\nIMPORTANT: No reported KD has been automatically validated.")
