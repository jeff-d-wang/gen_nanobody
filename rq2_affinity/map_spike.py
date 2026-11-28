from pathlib import Path
from Bio.PDB import MMCIFParser, is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
from scipy.spatial import cKDTree
import pandas as pd
import numpy as np
import requests
import json
import re

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
STRUCT = ROOT / "structures/spike"
STRUCT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(OUT / "spike_detailed_audit.csv")
df = df.drop_duplicates(["PDB_ID", "vhh_key", "kd_M"])

aligner = PairwiseAligner(mode="local")
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -6
aligner.extend_gap_score = -0.5

AA = set("ACDEFGHIKLMNPQRSTVWY")

def clean(x):
    if pd.isna(x):
        return ""
    s = re.sub(r"\s+", "", str(x)).upper()
    return s if len(s) > 30 and set(s) <= AA else ""

def match(ref, observed):
    if not ref or not observed:
        return 0., 0.
    a = aligner.align(ref, observed)[0]
    aligned = 0
    identical = 0
    for (i, j), (k, l) in zip(*a.aligned):
        n = min(j-i, l-k)
        aligned += n
        identical += sum(
            ref[i+t] == observed[k+t] for t in range(n)
        )
    return (
        identical / aligned if aligned else 0.,
        aligned / len(ref)
    )

def chain_info(chain):
    residues = []
    letters = []
    for r in chain:
        if is_aa(r, standard=False) and "CA" in r:
            residues.append(r)
            letters.append(
                seq1(r.resname, custom_map={"MSE": "M"},
                     undef_code="X")
            )
    return "".join(letters), residues

def coords(residues):
    xyz, rid = [], []
    for i, r in enumerate(residues):
        for atom in r:
            if atom.element in ("C", "N", "O", "S"):
                xyz.append(atom.coord)
                rid.append(i)
    return np.asarray(xyz).reshape(-1, 3), np.asarray(rid)

def interface(vhh, ag):
    a, ai = coords(vhh)
    b, bi = coords(ag)
    if not len(a) or not len(b):
        return dict(atom_contacts=0, residue_contacts=0,
                    clashes=0, vhh_residues=0, ag_residues=0)

    tree = cKDTree(b)
    pairs = [
        (i, j)
        for i, js in enumerate(tree.query_ball_point(a, 5.0))
        for j in js
    ]
    rp = {(int(ai[i]), int(bi[j])) for i, j in pairs}
    return dict(
        atom_contacts=len(pairs),
        residue_contacts=len(rp),
        vhh_residues=len({x for x, y in rp}),
        ag_residues=len({y for x, y in rp}),
        clashes=sum(
            len(x) for x in tree.query_ball_point(a, 2.0)
        )
    )

parser = MMCIFParser(QUIET=True)
cache = {}
results = []
chain_checks = []

for _, row in df.iterrows():
    pdb = str(row["PDB_ID"]).upper()
    ref_v = clean(row["Ab/Nano H_Chain AA"])
    ref_a = clean(row["Ag_Seq"])
    base = dict(
        pdb_id=pdb,
        vhh_key=row["vhh_key"],
        kd_M=float(row["kd_M"]),
        pKD=float(row["pKD"])
    )

    try:
        if pdb not in cache:
            path = STRUCT / f"{pdb.lower()}.cif"
            if not path.exists():
                response = requests.get(
                    f"https://files.rcsb.org/download/{pdb}.cif",
                    timeout=60
                )
                response.raise_for_status()
                if len(response.content) < 1000:
                    raise ValueError("Empty mmCIF response")
                path.write_bytes(response.content)

            model = next(
                parser.get_structure(pdb, str(path)).get_models()
            )
            cache[pdb] = {
                c.id: chain_info(c)
                for c in model
            }

        chains = cache[pdb]
        vhhs, antigens = [], []

        for cid, (sequence, residues) in chains.items():
            vi, vc = match(ref_v, sequence)
            ai, ac = match(ref_a, sequence)

            chain_checks.append(dict(
                pdb_id=pdb, chain=cid,
                vhh_identity=vi, vhh_coverage=vc,
                antigen_identity=ai, antigen_coverage=ac
            ))

            if vi >= .90 and vc >= .80:
                vhhs.append((cid, residues, vi, vc))
            if ai >= .90 and ac >= .80:
                antigens.append((cid, residues, ai, ac))

        candidates = []
        for vid, vr, vi, vc in vhhs:
            for aid, ar, ai, ac in antigens:
                if vid == aid:
                    continue
                candidates.append(dict(
                    **base,
                    vhh_chain=vid,
                    antigen_chain=aid,
                    vhh_identity=vi,
                    vhh_coverage=vc,
                    antigen_identity=ai,
                    antigen_coverage=ac,
                    **interface(vr, ar)
                ))

        if not candidates:
            results.append(dict(
                **base, status="NO_VALID_CHAIN_PAIR",
                vhh_candidates=len(vhhs),
                antigen_candidates=len(antigens)
            ))
            print(pdb, "NO_VALID_CHAIN_PAIR")
            continue

        # Select an interacting copy for exploratory geometry.
        # This selection is NOT independent affinity prediction.
        best = max(
            candidates,
            key=lambda x: (
                x["residue_contacts"], x["atom_contacts"]
            )
        )
        best["candidate_chain_pairs"] = len(candidates)
        best["status"] = (
            "MAPPED_CONTACT"
            if best["residue_contacts"] > 0
            else "MAPPED_NO_CONTACT"
        )
        results.append(best)
        print(
            pdb, best["status"],
            f'{best["vhh_chain"]}/{best["antigen_chain"]}',
            "residue contacts:", best["residue_contacts"],
            "clashes:", best["clashes"]
        )

    except Exception as e:
        results.append(dict(
            **base, status="ERROR", error=str(e)
        ))
        print(pdb, "ERROR", str(e))

out = pd.DataFrame(results)
out.to_csv(OUT / "spike_interface_audit.csv", index=False)
pd.DataFrame(chain_checks).to_csv(
    OUT / "spike_chain_mapping.csv", index=False
)

summary = {
    "analyzed_pdb_vhh_kd_records": int(len(out)),
    "status_counts": out["status"].value_counts().to_dict(),
    "mapped_unique_vhhs": int(
        out.loc[
            out["status"] == "MAPPED_CONTACT", "vhh_key"
        ].nunique()
    ),
    "warning": (
        "Exploratory structural descriptors only. "
        "Biological assembly, assay provenance and repeated KD "
        "measurements remain to be verified."
    )
}

(OUT / "spike_mapping_summary.json").write_text(
    json.dumps(summary, indent=2)
)
print("\n=== FINAL SUMMARY ===")
print(json.dumps(summary, indent=2))
