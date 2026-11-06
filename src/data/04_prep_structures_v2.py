from pathlib import Path
from collections import defaultdict
import csv
import traceback

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


ROOT = Path.cwd()
TARGETS = ROOT / "results" / "targets.csv"
OUTDIR = ROOT / "results" / "prepped_v2"

EPITOPE_CUTOFF = 5.0
CROP_CUTOFF = 15.0

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "MSE", "HYP", "CSO",
    "CSD", "CME", "MLY", "M3L", "KCX",
    "LLP", "FME", "SEP", "TPO", "PTR",
    "HIC", "PCA", "CYX", "ORN", "SAR",
    "NLE", "ABA",
}

BACKBONE = {"N", "CA", "C", "O"}
HEAVY_ELEMENTS = {"C", "N", "O", "S", "P", "SE"}


def split_chains(value):
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split("/") if x.strip()]


def is_blank_altloc(value):
    return value in ("", ".", "?", "\x00", " ")


def residue_key(chain_id, residue):
    return (
        chain_id,
        residue.seqid.num,
        str(residue.seqid.icode).replace("\x00", "").strip(),
        residue.name,
    )


def choose_residue_atoms(residue):
    """
    Keep blank-altloc atoms plus one coherent nonblank conformer.
    Select the conformer with the highest summed occupancy.
    For each atom name, prefer the selected altloc over blank
    if both are present.
    """
    atoms = list(residue)
    altloc_scores = defaultdict(float)

    for atom in atoms:
        alt = str(atom.altloc)
        if not is_blank_altloc(alt):
            altloc_scores[alt] += float(atom.occ)

    selected_alt = None
    if altloc_scores:
        selected_alt = sorted(
            altloc_scores,
            key=lambda x: (-altloc_scores[x], x)
        )[0]

    chosen = {}
    for atom in atoms:
        alt = str(atom.altloc)

        if is_blank_altloc(alt):
            candidate = (0, atom)
        elif alt == selected_alt:
            candidate = (1, atom)
        else:
            continue

        name = atom.name.strip()
        if name not in chosen or candidate[0] > chosen[name][0]:
            chosen[name] = candidate

    return [entry[1] for entry in chosen.values()], selected_alt


def heavy_atom(atom):
    return atom.element.name.upper() in HEAVY_ELEMENTS


def get_residue_records(chain, chain_id):
    records = []

    for residue in chain:
        resname = residue.name.upper()
        if resname not in PROTEIN_RESIDUES:
            continue

        atoms, selected_alt = choose_residue_atoms(residue)
        atoms = [a for a in atoms if heavy_atom(a)]

        if not atoms:
            continue

        coords = np.array(
            [[a.pos.x, a.pos.y, a.pos.z] for a in atoms],
            dtype=float
        )

        names = [a.name.strip() for a in atoms]

        records.append({
            "key": residue_key(chain_id, residue),
            "chain": chain_id,
            "resid": residue.seqid.num,
            "icode": str(residue.seqid.icode).replace("\x00", "").strip(),
            "resname": residue.name,
            "atoms": atoms,
            "coords": coords,
            "atom_names": names,
            "selected_altloc": selected_alt,
        })

    return records


def min_distance(coords_a, coords_b):
    if len(coords_a) == 0 or len(coords_b) == 0:
        return float("inf")

    best = float("inf")
    for start in range(0, len(coords_a), 2000):
        d = cdist(coords_a[start:start + 2000], coords_b)
        best = min(best, float(d.min()))
    return best


def make_atom(source):
    atom = gemmi.Atom()
    atom.name = source.name
    atom.element = source.element
    atom.pos = source.pos
    atom.occ = source.occ
    atom.b_iso = source.b_iso
    atom.altloc = "\x00"
    return atom


def make_residue(record, output_number, insertion_code=" "):
    residue = gemmi.Residue()
    residue.name = record["resname"]
    residue.seqid = gemmi.SeqId(output_number, insertion_code)

    for atom in record["atoms"]:
        residue.add_atom(make_atom(atom))

    return residue


def write_structure(path, h_records, antigen_records, antigen_map):
    structure = gemmi.Structure()
    structure.name = path.stem

    model = gemmi.Model("1")

    h_chain = gemmi.Chain("H")
    for record in h_records:
        residue = make_residue(
            record,
            record["resid"],
            record["icode"] if record["icode"] else " "
        )
        h_chain.add_residue(residue)

    t_chain = gemmi.Chain("T")
    for record in antigen_records:
        new_number = antigen_map[record["key"]]
        residue = make_residue(record, new_number)
        t_chain.add_residue(residue)

    model.add_chain(h_chain)
    model.add_chain(t_chain)
    structure.add_model(model)

    structure.setup_entities()
    structure.write_pdb(str(path))


def prepare_one(row):
    pdb_id = str(row["PDB_ID"])
    h_chain_id = str(row["Hchain"])
    antigen_chain_ids = split_chains(row["agchains"])
    cif_path = Path(row["cif"])

    if not cif_path.is_file():
        raise FileNotFoundError(cif_path)

    structure = gemmi.read_structure(str(cif_path))
    model = structure[0]

    available = {chain.name for chain in model}

    if h_chain_id not in available:
        raise ValueError(
            f"Nanobody chain {h_chain_id} not found; "
            f"available={sorted(available)}"
        )

    missing = [c for c in antigen_chain_ids if c not in available]
    if missing:
        raise ValueError(f"Missing antigen chains: {missing}")

    h_records = get_residue_records(model[h_chain_id], h_chain_id)
    antigen_records = []

    for chain_id in antigen_chain_ids:
        antigen_records.extend(
            get_residue_records(model[chain_id], chain_id)
        )

    if not h_records:
        raise ValueError("No nanobody protein residues")

    if not antigen_records:
        raise ValueError("No antigen protein residues")

    h_coords = np.concatenate([r["coords"] for r in h_records])

    epitope_records = []
    for record in antigen_records:
        distance = min_distance(record["coords"], h_coords)
        record["distance_to_h"] = distance
        if distance <= EPITOPE_CUTOFF:
            epitope_records.append(record)

    if not epitope_records:
        raise ValueError("No antigen epitope residues found")

    epitope_coords = np.concatenate(
        [r["coords"] for r in epitope_records]
    )

    crop_records = []
    for record in antigen_records:
        distance = min_distance(record["coords"], epitope_coords)
        record["distance_to_epitope"] = distance
        if distance <= CROP_CUTOFF:
            crop_records.append(record)

    # Preserve original chain/residue identity in a separate mapping.
    antigen_map = {
        record["key"]: i + 1
        for i, record in enumerate(crop_records)
    }

    output_path = OUTDIR / f"{pdb_id}_H_T.pdb"
    write_structure(
        output_path,
        h_records,
        crop_records,
        antigen_map
    )

    epitope_rows = []
    for record in epitope_records:
        epitope_rows.append({
            "PDB_ID": pdb_id,
            "original_chain": record["chain"],
            "original_res_id": record["resid"],
            "original_insertion_code": record["icode"],
            "original_res_name": record["resname"],
            "min_distance_to_nanobody_A": record["distance_to_h"],
            "epitope_cutoff_A": EPITOPE_CUTOFF,
        })

    mapping_rows = []
    for record in crop_records:
        mapping_rows.append({
            "PDB_ID": pdb_id,
            "original_chain": record["chain"],
            "original_res_id": record["resid"],
            "original_insertion_code": record["icode"],
            "original_res_name": record["resname"],
            "output_chain": "T",
            "output_res_id": antigen_map[record["key"]],
        })

    # Report excluded non-protein antigen components.
    nonprotein_rows = []
    for chain_id in antigen_chain_ids:
        for residue in model[chain_id]:
            if residue.name.upper() not in PROTEIN_RESIDUES:
                nonprotein_rows.append({
                    "PDB_ID": pdb_id,
                    "chain": chain_id,
                    "res_id": residue.seqid.num,
                    "insertion_code": str(
                        residue.seqid.icode
                    ).replace("\x00", "").strip(),
                    "res_name": residue.name,
                })

    # Altloc selection audit.
    altloc_rows = []
    for record in h_records + crop_records:
        if record["selected_altloc"] is not None:
            altloc_rows.append({
                "PDB_ID": pdb_id,
                "chain": record["chain"],
                "res_id": record["resid"],
                "insertion_code": record["icode"],
                "res_name": record["resname"],
                "selected_altloc": record["selected_altloc"],
            })

    summary = {
        "PDB_ID": pdb_id,
        "Hchain_original": h_chain_id,
        "antigen_chains_original": "/".join(antigen_chain_ids),
        "nanobody_residues": len(h_records),
        "epitope_residues": len(epitope_records),
        "antigen_residues_after_crop": len(crop_records),
        "nanobody_atoms": sum(len(r["atoms"]) for r in h_records),
        "antigen_atoms_after_crop": sum(
            len(r["atoms"]) for r in crop_records
        ),
        "output": str(output_path),
        "status": "OK",
    }

    return summary, epitope_rows, mapping_rows, nonprotein_rows, altloc_rows


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    targets = pd.read_csv(TARGETS)
    required = ["PDB_ID", "Hchain", "agchains", "cif"]
    missing = [c for c in required if c not in targets.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    summaries = []
    epitopes = []
    mappings = []
    nonproteins = []
    altlocs = []
    failures = []

    for _, row in targets.iterrows():
        pdb_id = str(row["PDB_ID"])
        print(f"\nProcessing {pdb_id}", flush=True)

        try:
            result = prepare_one(row)
            summary, ep, mp, np_rows, al = result

            summaries.append(summary)
            epitopes.extend(ep)
            mappings.extend(mp)
            nonproteins.extend(np_rows)
            altlocs.extend(al)

            print(
                f"OK | H={summary['nanobody_residues']} residues "
                f"| epitope={summary['epitope_residues']} "
                f"| crop={summary['antigen_residues_after_crop']}"
            )

        except Exception as exc:
            failures.append({
                "PDB_ID": pdb_id,
                "error": repr(exc),
                "traceback": traceback.format_exc(),
            })
            print(f"FAILED: {exc}")

    pd.DataFrame(summaries).to_csv(
        ROOT / "results" / "prepped_v2.csv", index=False
    )
    pd.DataFrame(epitopes).to_csv(
        ROOT / "results" / "epitope_residues_v2.csv", index=False
    )
    pd.DataFrame(mappings).to_csv(
        ROOT / "results" / "antigen_residue_mapping_v2.csv",
        index=False
    )
    pd.DataFrame(nonproteins).to_csv(
        ROOT / "results" / "nonprotein_components_v2.csv",
        index=False
    )
    pd.DataFrame(altlocs).to_csv(
        ROOT / "results" / "altloc_choices_v2.csv",
        index=False
    )
    pd.DataFrame(failures).to_csv(
        ROOT / "results" / "prep_failures_v2.csv",
        index=False
    )

    print("\n" + "=" * 55)
    print(f"Total: {len(targets)}")
    print(f"Successful: {len(summaries)}")
    print(f"Failed: {len(failures)}")
    print(f"Output: {OUTDIR}")
    print("=" * 55)


if __name__ == "__main__":
    main()
