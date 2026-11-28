from pathlib import Path
import pandas as pd
import re

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"

evidence = pd.read_csv(
    OUT / "rq2_affinity_text_evidence_extended.csv",
    dtype=str
).fillna("")

observations = pd.read_csv(
    OUT / "rq2_affinity_evidence_audit.csv",
    dtype={"pubmed_id": str}
)

papers = pd.read_csv(
    OUT / "rq2_publication_metadata.csv",
    dtype={"pubmed_id": str}
)

# Numeric affinity candidates. These are search clues only.
numeric_pattern = re.compile(
    r"(?i)(?:"
    r"\b\d+(?:\.\d+)?(?:\s*[×x]\s*10\s*[-−]?\d+)?"
    r"\s*(?:fM|pM|nM|µM|μM|uM|mM)\b"
    r"|"
    r"\bK\s*[_-]?\s*d\b.{0,35}\d"
    r")"
)

method_pattern = re.compile(
    r"(?i)\b(?:SPR|ITC|BLI|KinExA|"
    r"surface plasmon resonance|"
    r"biolayer interferometry|"
    r"isothermal titration calorimetry)\b"
)

evidence["numeric_clue"] = evidence["snippet"].map(
    lambda x: bool(numeric_pattern.search(x))
)
evidence["method_clue"] = evidence["snippet"].map(
    lambda x: bool(method_pattern.search(x))
)

evidence["priority"] = (
    evidence["numeric_clue"].astype(int) * 2
    + evidence["method_clue"].astype(int)
)

evidence = evidence.sort_values(
    ["pubmed_id", "priority"],
    ascending=[True, False]
)

review_rows = []
report = [
    "# RQ2 Affinity Evidence Review",
    "",
    "Automatic text matching; no KD value verified.",
    ""
]

for pmid, group in observations.groupby("pubmed_id"):
    pmid = str(pmid)
    unique = group.drop_duplicates(
        ["antigen", "pdb_id", "vhh_key", "kd_M"]
    )

    paper = papers[
        papers["pubmed_id"].astype(str) == pmid
    ]
    title = (
        str(paper.iloc[0].get("title", ""))
        if len(paper) else ""
    )

    matches = evidence[evidence["pubmed_id"] == pmid]

    report += [
        f"## PMID {pmid}",
        "",
        f"Title: {title}",
        "",
        "### ANDD observations",
        ""
    ]

    for _, row in unique.iterrows():
        kd_nM = float(row["kd_M"]) * 1e9
        report.append(
            f"- {row['antigen']} | {row['pdb_id']} | "
            f"VHH {str(row['vhh_key'])[:10]} | "
            f"KD = {kd_nM:.6g} nM"
        )

    report += [
        "",
        f"### Retrieved text clues ({len(matches)})",
        ""
    ]

    for _, clue in matches.head(15).iterrows():
        snippet = " ".join(str(clue["snippet"]).split())
        report.append(
            f"- [priority {clue['priority']}] {snippet}"
        )

        review_rows.append({
            "pubmed_id": pmid,
            "title": title,
            "snippet": snippet,
            "numeric_clue": clue["numeric_clue"],
            "method_clue": clue["method_clue"],
            "review_status": "UNVERIFIED",
        })

    report += [
        "",
        "Verification: PENDING — exact construct, "
        "reported KD, method and conditions",
        "",
    ]

(OUT / "rq2_kd_evidence_review.md").write_text(
    "\n".join(report),
    encoding="utf-8"
)

pd.DataFrame(review_rows, columns=[
    "pubmed_id", "title", "snippet",
    "numeric_clue", "method_clue", "review_status"
]).to_csv(
    OUT / "rq2_kd_evidence_prioritized.csv",
    index=False
)

print("\n=== TEXT EVIDENCE SCREENING ===")
print("Total passages:", len(evidence))
print("Numeric affinity clues:",
      int(evidence["numeric_clue"].sum()))
print("Method clues:",
      int(evidence["method_clue"].sum()))

print("\n=== BY PUBLICATION ===")
print(
    evidence.groupby("pubmed_id")
    .agg(
        passages=("snippet", "size"),
        numeric=("numeric_clue", "sum"),
        method=("method_clue", "sum")
    )
    .to_string()
)

print("\nSaved:")
print(OUT / "rq2_kd_evidence_review.md")
print(OUT / "rq2_kd_evidence_prioritized.csv")
