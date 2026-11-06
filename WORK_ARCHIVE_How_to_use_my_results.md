# Work Archive: How to Use My Results

**Contributor:** Charlie (C)  
**Project:** Nanobody Design  
**Branch:** `charlie/benchmark-midterm`  
**Scope:** Structural benchmark construction, preprocessing, QC, benchmark analysis, and midterm figures.

All paths below are relative to the repository root. The data products described here are intended to be directly usable from the repository.

---

## 1. Overview

This contribution provides a ready-to-use structural benchmark of 21 nanobody–antigen complexes.

Included deliverables:

- 21 standardized complex structures.
- 21 ProteinMPNN-compatible structure files.
- Target metadata and reference nanobody sequences.
- Epitope annotations and antigen residue mappings.
- Structural quality-control records.
- Benchmark analysis tables.
- Midterm figures and summary statistics.
- Structural preprocessing code.

No access to the contributor's personal working files is required to use the committed data products.

## 2. Benchmark Construction

### 2.1 Source and filtering

The benchmark was constructed from SAbDAb2 ML v0.2.0 metadata (`abag_split_sd.csv`).

Filtering criteria:

- `type == "SD-H"`
- `holo == True`
- Protein antigen
- Nonempty antigen-chain annotation
- `ab_ag_split == "test"`
- PDB deposition date later than 2023-01-01
- Best-resolution structure selected per `ab_ag_cluster`

The final benchmark contains **21 complexes**.

PDB IDs:

`8qf4, 8yj8, 9ty2, 9g2g, 8idi, 9h77, 9min, 9my8, 9dc3, 9rxi, 9b71, 9mnz, 9q6n, 9yc3, 9npt, 9rye, 8rwf, 9g3e, 9ljj, 9rmk, 9e7l`

### 2.2 Target metadata

**File:** `data/processed/targets.csv`

Contains:

- `ab_ag_cluster`
- `PDB_ID`
- `Hchain`
- `agchains`
- `agtypes`
- `resolution`
- `PDBdepo`
- `CDRH1`, `CDRH2`, `CDRH3`
- `Hseq`
- `cif`
- `cif_exists`

Use this file as the target-level metadata reference for target identity, chain selection, reference sequence, and CDR annotations.

The `cif` field records source metadata. Downstream workflows should use the committed standardized structures rather than depending on original source paths.

## 3. Standardized Structures

### 3.1 Canonical v2 structures

**Directory:** `data/processed/prepped_v2/`  
**Index:** `data/processed/prepped_v2.csv`

These are the canonical standardized nanobody–antigen complex structures.

Preprocessing steps:

1. Read mmCIF structures using Gemmi.
2. Retain protein residues and heavy atoms.
3. Resolve alternate locations using summed occupancy.
4. Preserve nanobody H-chain residue numbering and insertion codes.
5. Merge antigen chains into chain T and renumber the merged antigen sequentially.
6. Define epitope residues as antigen residues within 5 Å of the nanobody.
7. Crop antigen residues to within 15 Å of the epitope.

All 21 selected complexes completed preprocessing.

Use `prepped_v2/` for standard structural inspection and structural evaluation.

### 3.2 Epitope annotations and antigen mapping

| File | Description |
|---|---|
| `data/processed/epitope_residues_v2.csv` | Antigen residues within 5 Å of the nanobody |
| `data/processed/antigen_residue_mapping_v2.csv` | Correspondence for cropped antigen residues |

The epitope is operationally defined using a 5 Å distance cutoff.

Use the antigen mapping to translate cropped antigen residue identifiers to their original identifiers.

Do not assume sequential chain-T numbering is identical to original antigen-chain numbering.

## 4. ProteinMPNN-Compatible Structures

**Directory:** `data/processed/prepped_mpnn/`  
**Mapping:** `data/processed/mpnn_residue_map.csv`

A separate structure derivative was generated to address an insertion-code parsing issue observed in `8qf4`.

### 4.1 Parser issue

The standard v2 structures preserve original H-chain numbering and insertion codes.

The official ProteinMPNN parser was observed to mishandle insertion-code ordering in `8qf4`. The reference CDRH3 sequence `ISKG` was consequently parsed as `GKSI`.

The ProteinMPNN-compatible derivative:

- Iterates residues in structural order using Gemmi.
- Renumbers residues sequentially.
- Removes insertion codes.
- Records correspondence to original residue identifiers in `mpnn_residue_map.csv`.

The parsed H1, H2, and H3 sequences were checked against target metadata and matched for all 21 structures.

### 4.2 Recommended input

| Task | Input |
|---|---|
| Standard structural analysis | `data/processed/prepped_v2/` |
| ProteinMPNN parsing and design | `data/processed/prepped_mpnn/` |
| Translate ProteinMPNN indices to original residue identifiers | `data/processed/mpnn_residue_map.csv` |

**Important:** ProteinMPNN sequential residue indices must not be interpreted as original PDB residue numbers. Use the mapping file for residue-level correspondence.

The ProteinMPNN-compatible derivative does not replace the canonical v2 structures.

## 5. Structural Quality Control

QC files:

| File | Description |
|---|---|
| `data/processed/structural_qc_final.csv` | Final structural QC summary |
| `data/processed/prep_failures_v2.csv` | Preprocessing failure records |
| `data/processed/altloc_choices_v2.csv` | Alternate-location decisions |
| `data/processed/nonprotein_components_v2.csv` | Non-protein component records |

### 5.1 QC summary

- 20 structures: PASS
- 1 structure: CHECK (`8yj8`)
- H-chain sequence status: 15 MATCH, 6 TERMINAL_UNRESOLVED
- H-chain duplicate residue checks: 0
- T-chain duplicate residue checks: 0
- H-chain C–N geometry anomalies: 1
- T-chain C–N geometry anomalies: 0

### 5.2 Important case: 8yj8

The selected alternate locations for ASP113-B and PHE114-C produce a C–N distance of approximately 2.194 Å.

The structure was retained and flagged for inspection, not silently repaired.

The six `TERMINAL_UNRESOLVED` sequence statuses concern unresolved terminal residues or tags. They should not automatically be interpreted as internal sequence mismatches.

Inspect QC records before making geometry-sensitive conclusions.

## 6. Benchmark Analysis

### 6.1 Analysis tables

| File | Description |
|---|---|
| `results/benchmark_analysis/cropping_comparison.csv` | Full-antigen versus cropped-antigen residue counts and cropping statistics |
| `results/benchmark_analysis/cdr_contact_analysis.csv` | CDR/framework contact analysis |

### 6.2 Benchmark summary

| Metric | Result |
|---|---:|
| Number of complexes | 21 |
| Resolution mean | 2.501 Å |
| Resolution median | 2.70 Å |
| Resolution range | 1.02–3.33 Å |
| Mean nanobody length | 122.286 residues |
| Mean full antigen length | 370.333 residues |
| Mean cropped antigen length | 145.476 residues |
| Mean epitope size | 21.238 residues |
| Median epitope size | 21 residues |
| Mean cropped/full antigen fraction | 0.584837 |
| Mean antigen size reduction | 0.415163 |

### 6.3 CDR contact fractions

At a 5 Å cutoff:

- H1: 0.373
- H2: 0.504
- H3: 0.540
- Framework: 0.048

These values represent the fraction of residues within each region that have at least one contact. They are not fractions of the total interface contributed by each region.

## 7. Midterm Figures and Tables

**Directory:** `results/midterm_figures/`

| File | Description |
|---|---|
| `Figure1_benchmark.pdf` | Benchmark composition and structural/contact summaries |
| `Figure1_benchmark.png` | PNG version of Figure 1 |
| `Figure2_cropping.pdf` | Antigen cropping analysis |
| `Figure2_cropping.png` | PNG version of Figure 2 |
| `Table1_target_statistics.csv` | Per-target statistics |
| `Table1_summary.csv` | Aggregate benchmark statistics |

Interpretation notes:

- Figure 1C represents within-region contact fraction.
- Figure 2B represents cropped/full antigen fraction, not epitope retention.
- The epitope-to-crop ratio is not a strict epitope-retention metric.

## 8. Downstream Usage

### 8.1 Sequence design and ProteinMPNN

Recommended starting files:

1. `data/processed/targets.csv`
2. `data/processed/prepped_mpnn/`
3. `data/processed/mpnn_residue_map.csv`
4. `data/processed/epitope_residues_v2.csv`
5. `data/processed/structural_qc_final.csv`

Suggested workflow:

1. Read target metadata to identify the target, reference sequence, H chain, antigen chains, and CDR annotations.
2. Use `prepped_mpnn/` as input to the ProteinMPNN parser workflow.
3. Verify that parsed H1, H2, and H3 sequences agree with the reference annotations.
4. Preserve `mpnn_residue_map.csv` for interpreting ProteinMPNN residue indices.
5. Consult epitope and antigen mapping tables for interface interpretation.
6. Review structural QC before geometry-sensitive evaluation.

The structures and metadata needed to begin this workflow are included in this repository.

### 8.2 Structural evaluation

For structural evaluation:

- Use `prepped_v2/` as the canonical structural reference.
- Use `epitope_residues_v2.csv` for the 5 Å epitope definition.
- Use `antigen_residue_mapping_v2.csv` to interpret cropped antigen identifiers.
- Use `structural_qc_final.csv` to identify structures requiring inspection.

## 9. Limitations and Cautions

1. The benchmark contains 21 complexes and is not a comprehensive sample of all nanobody complexes.
2. One representative structure was selected per cluster.
3. The epitope is operationally defined using a 5 Å cutoff.
4. The antigen crop retains residues within 15 Å of the defined epitope.
5. `8yj8` has a retained C–N geometry anomaly and should be inspected.
6. Terminal unresolved sequence statuses require terminal-context interpretation.
7. The ProteinMPNN-compatible derivative addresses the parser behavior observed in `8qf4`.
8. No claim is made that every v2 coordinate is universally more accurate than every v1 coordinate.
9. Available cropping summaries do not establish strict epitope retention through a complete original-to-cropped residue audit.

## 10. Completion Status

Completed:

- Benchmark selection and construction.
- Standardized v2 structural preprocessing for 21 complexes.
- ProteinMPNN-compatible structure preparation.
- Residue mapping.
- Structural QC.
- Epitope annotation.
- Cropping and CDR-contact analysis.
- Midterm figures and summary tables.
- Documentation of downstream usage and known caveats.

The committed data products are intended to provide a ready-to-use starting point for downstream sequence design and structural evaluation.
