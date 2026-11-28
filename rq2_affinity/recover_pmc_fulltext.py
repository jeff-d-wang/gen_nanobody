from pathlib import Path
import re
import time
import xml.etree.ElementTree as ET

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
CACHE = ROOT / "raw" / "pmc_bioc"
CACHE.mkdir(parents=True, exist_ok=True)

papers_path = OUT / "rq2_publication_metadata.csv"
evidence_path = OUT / "rq2_affinity_text_evidence.csv"

papers = pd.read_csv(papers_path, dtype=str).fillna("")
previous = pd.read_csv(evidence_path, dtype=str).fillna("")

session = requests.Session()
session.headers.update({
    "User-Agent": "NanobodyRQ2AcademicAudit/1.0"
})

retry = Retry(
    total=2,
    backoff_factor=1.0,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)

session.mount(
    "https://",
    HTTPAdapter(max_retries=retry)
)

AFFINITY = re.compile(
    r"(?i)(?:\baffinit(?:y|ies)\b|"
    r"\bdissociation\b|\bbinding kinetics\b|"
    r"\bSPR\b|\bITC\b|\bBLI\b|"
    r"\bsurface plasmon\b|\bkinexa\b|"
    r"\bK\s*[_-]?\s*[dD]\b)"
)

UNIT = re.compile(
    r"(?i)(?:\b\d+(?:\.\d+)?\s*"
    r"(?:pM|nM|uM|µM|μM|mM)\b|"
    r"\bK\s*[_-]?\s*[dD]\s*[=:≈])"
)

new_evidence = []
status_updates = []

for _, paper in papers.iterrows():
    pmid = paper["pubmed_id"]
    pmcid = paper.get("pmcid", "")
    old_status = paper.get("fulltext_status", "")

    if old_status == "retrieved":
        continue

    if not pmcid.startswith("PMC"):
        status_updates.append({
            "pubmed_id": pmid,
            "pmcid": pmcid,
            "recovery_status": "no_pmcid",
            "passages": 0,
            "error": "",
        })
        continue

    url = (
        "https://www.ncbi.nlm.nih.gov/research/bionlp/"
        f"RESTful/pmcoa.cgi/BioC_xml/{pmcid}/unicode"
    )

    print("FETCH", pmid, pmcid, flush=True)

    try:
        cached = CACHE / f"{pmcid}.xml"

        if cached.exists():
            content = cached.read_bytes()
        else:
            response = session.get(url, timeout=(10, 30))
            response.raise_for_status()
            content = response.content

        root = ET.fromstring(content)

        passages = []
        for element in root.iter():
            if element.tag.rsplit("}", 1)[-1] != "passage":
                continue

            parts = [
                x.text or ""
                for x in element.iter()
                if x.tag.rsplit("}", 1)[-1] == "text"
            ]
            passage = " ".join(" ".join(parts).split())

            if passage:
                passages.append(passage)

        if not passages:
            raise ValueError(
                "BioC returned no readable full-text passages"
            )

        if not cached.exists():
            cached.write_bytes(content)

        seen = set()
        hits = []

        for passage in passages:
            match = AFFINITY.search(passage)
            if not match:
                continue

            start = max(0, match.start() - 100)
            end = min(len(passage), match.end() + 250)
            snippet = passage[start:end]

            if snippet not in seen:
                seen.add(snippet)
                hits.append({
                    "pubmed_id": pmid,
                    "passage_number": len(hits) + 1,
                    "snippet": snippet,
                    "evidence_type": "bioc_keyword_context_only",
                    "contains_numeric_affinity_clue": bool(
                        UNIT.search(passage)
                    ),
                    "verified_kd": False,
                })

        new_evidence.extend(hits)

        status_updates.append({
            "pubmed_id": pmid,
            "pmcid": pmcid,
            "recovery_status": "bioc_retrieved",
            "passages": len(hits),
            "error": "",
        })

        print(
            "SUCCESS", pmid,
            "affinity passages:", len(hits),
            flush=True
        )

    except Exception as exc:
        status_updates.append({
            "pubmed_id": pmid,
            "pmcid": pmcid,
            "recovery_status": "bioc_unavailable",
            "passages": 0,
            "error": str(exc)[:300],
        })
        print("UNAVAILABLE", pmid, str(exc)[:160])

    time.sleep(0.5)

recovery = pd.DataFrame(status_updates)
recovery.to_csv(
    OUT / "rq2_fulltext_recovery_status.csv",
    index=False
)

new = pd.DataFrame(new_evidence)

if not new.empty:
    combined = pd.concat(
        [previous, new],
        ignore_index=True
    )
    combined = combined.drop_duplicates(
        subset=["pubmed_id", "snippet"]
    )
else:
    combined = previous.copy()

combined.to_csv(
    OUT / "rq2_affinity_text_evidence_extended.csv",
    index=False
)

print("\n=== FULLTEXT RECOVERY SUMMARY ===")
if len(recovery):
    print(
        recovery["recovery_status"]
        .value_counts().to_string()
    )
    print("\n=== DETAILS ===")
    print(
        recovery[
            ["pubmed_id", "pmcid",
             "recovery_status", "passages"]
        ].to_string(index=False)
    )

print("\nNew keyword passages:", len(new))
print("Total saved passages:", len(combined))
print("\nOutputs:")
print("rq2_fulltext_recovery_status.csv")
print("rq2_affinity_text_evidence_extended.csv")
print("\nNo KD values automatically verified.")
