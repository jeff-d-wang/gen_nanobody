# RQ2 Training-Free Feature Rank Benchmark — Predeclared v1

## Purpose

Test whether **existing** experimental-complex geometry features preserve *within-antigen* ordering of source-reported ANDD KD values. **No model training**, no absolute KD prediction, no feature-weight tuning, no result-driven selection of antigens or structures. This is exploratory until label provenance / assay comparability is validated, and is NOT a substitute for evaluating independently predicted complexes from generated sequences.

## Data source

Read-only reuse of `~/Research/nanobody_rq2_clean/tables/{structure_scores.csv,pair_replicate_audit.csv}` and optional `verification/publication_review.csv` from the preceding exhaustive run. Do not re-run RCSB downloads or geometry extraction. Never overwrite original outputs.

## Frozen scores

For each exact VHH–antigen pair in a group, evaluate `residue_contacts`, `atom_contacts`, `vhh_residues`, `ag_residues` (hypothetically higher indicates stronger binding), and `clashes` (hypothetically lower indicates stronger binding). The contact radius is inherited from original 5 Å geometry; clash cutoff 2 Å. Zero/constant scores are retained and shown as ties, **never dropped for performance**. Not evaluated: BSA, ipTM, PAE, pLDDT and Rosetta energy (not measured on the cached structures by this pipeline).

## Label integrity

- ANDD `real` is NOT publication-verified. Main outputs must say **exploratory/unverified**.
- Main exploratory tiers quarantine preexisting known-discrepancy watchlist, extreme-KD flags, and globally conflicting repeated KD for the identical sequence pair (max/min >1.10). Exclusions are logged.
- KD for consistent duplicate records: geometric median (median on log KD). Geometry for duplicate structures of same exact pair: median across **distinct PDBs**, independent of their KD. No pseudoreplication of PDB copies.
- If a strict group also has internal KD conflict >10%, exclude the exact pair from the group and record reason.
- A separate verified tier only accepts full reviewed annotations with verified experimental KD, assay group, affinity-paper DOI, construct identity and reviewer. Conflicted source labels may enter only via explicit human evidence, never by automatic correction.

## Grouping and primary metrics

1. `structure_paper_proxy_EXPLORATORY`: exact antigen sequence × PDB primary citation PMID or DOI (a *proxy*, never assumed to prove matching affinity assay). Includes all groups of n>=2 for pairwise comparisons. Spearman only shown for n>=5; no p-values.
2. `antigen_only_mixed_assays_SENSITIVITY`: exact antigen sequence, even across structural studies. A **sensitivity upper bound only**; measurements may not be comparable.
3. `source_verified_assay_only`: exact antigen sequence × verified affinity-paper DOI × verified assay group. No verified status is inferred from a PDB citation.

### Pre-fixed pairwise rule

Within a group, for distinct pairs i,j, with KD fold separation >= **3** (`0.4771 pKD`), compare sign of `(S_i-S_j)(pKD_i-pKD_j)`. When feature score is tied, give 0.5 concordance credit (neutral abstention); correct=1, reversed=0. Define per-group concordance as mean over informative pairs; **macro aggregate gives each group equal weight**. Micro pooled comparison score is reported *only as descriptive* because comparisons within groups overlap and are dependent. Report number of groups, informative comparisons, ties, and largest-group share. Pre-fixed sensitivity KD separations = 1.5x, 3x, 10x. Spearman for groups n>=5, no p-values. A group of n=2 may provide one pairwise comparison; it is NOT treated as independent statistical evidence.

### Formal go/no-go gate

Keep prior study's conservative gate: >=8 separate source-verified antigen×study×assay groups with n>=5, >=5 distinct antigens, >=50 eligible verified pairs, **before** considering group-level inference. Even passing this rough gate is not proof of statistical power or validity. If unavailable, label `NO_GO_UNVERIFIED` / `DESCRIPTIVE_ONLY`. No estimated KD, feature calibration, or transfer claim for the 84 RFantibody generated sequences.

## Mandatory outcomes

Produce **all** five feature results at 1.5x,3x,10x in both main/sensitivity tiers; full group census even if no informative pairs; full per-comparison audit; verified tier if available; report; JSON decision. Large-group dominance, low group count, many score ties, label/provenance gaps are reasons to lower confidence, not to revise the protocol.
