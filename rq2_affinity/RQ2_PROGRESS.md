# RQ2 Affinity Analysis — Progress Checkpoint

## Completed

- ANDD v2 downloaded and audited
- 48,800 ANDD records; 30,119 labeled Nanobody/VHH
- 1,645 real-labeled nanobody rows with positive Kd
- 1,304 complete VHH-antigen sequence records
- 928 unique complete sequence pairs
- 13 candidate antigen groups with >=5 VHH sequences
- Lysozyme: 11 structures mapped, 10 unique VHH
- 4PGJ excluded in strict camelid-VHH sensitivity cohort
- Lysozyme strict n=9, first PDB:
  - VHH residues rho=+0.7236, unadjusted p=0.0276
  - Antigen residues rho=-0.6441, unadjusted p=0.0612

## Frozen exploratory protocol

- Contact cutoff: 5 Angstrom
- Clash cutoff: 2 Angstrom
- Metrics: residue_contacts, atom_contacts,
  vhh_residues, ag_residues
- Within-antigen Spearman correlation with reported pKD
- Unique VHH as the analysis unit
- First/last PDB sensitivity analysis
- No metric or threshold tuning on validation cohorts

## Important limitations

- Missing affinity assay details in ANDD
- Original affinity publications not yet audited
- Repeated VHH measurements may differ
- Experimental crystal structures, not Boltz predictions
- Small cohort and multiple exploratory comparisons
- HER2 AbSci 422-record cohort excluded pending
  antibody-format verification
- Results are hypothesis-generating, not validated
  affinity prediction performance

## Next

1. Replicate fixed geometry pipeline for RNase and Spike
2. Check cross-antigen direction and sensitivity
3. Validate original affinity evidence and assay comparability
4. Write RQ2 results with explicit limitations
