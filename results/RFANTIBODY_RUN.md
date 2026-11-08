# RFantibody benchmark execution

- Run: `rfab_benchmark_1ce3344e4822`
- Complete conditions: 42/42 (21 targets, two conditions).
- Candidates passing lineage checks: 84/84.
- Recorded model-process time, including failed attempts: 4.34 hours.

Settings: two deterministic backbones per target/condition; CDR lengths 7/6/8; 50 diffusion steps; one ProteinMPNN sequence per backbone at temperature 0.1; RF2 with one recycle and seed 37.

Read `status.csv` for progress and failures, `candidate_metrics.csv` for RF2 scores, `designed_sequences.fasta` for sequences, and `audit.json` for validation details. Each target/condition folder contains per-stage attempt logs, structures, and `state.json`. `config.json`, `inputs/`, and `framework.pdb` identify the inputs.

Both conditions use the same reference-derived antigen crop. Only the hotspot condition receives explicit hotspot labels. This is a workflow feasibility run; no binding claims or threshold calibration are supported. RF2 is a pipeline predictor, not independent final scoring.

Raw confidence matrices are not exported by this runner; the saved RF2 PDB score lines and tables contain aggregate confidence metrics.
