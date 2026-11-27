# A1 cost accounting — implementation and results

Implemented November 26, 2026 against local branch `RFantibody`, source commit
`44d4b710c6ec9530113ce106d41f4ac0ba4baefa`. No remote branch or PR status is asserted.
This implements A1 in `Brief Status check & Next steps.md`; A2 is unchanged.

## Verified results

| Measure | All | Hotspot | Unguided |
|---|---:|---:|---:|
| Completed target/condition runs | 42 | 21 | 21 |
| Accepted generated backbones | 84 | 42 | 42 |
| Designed sequence records | 84 | 42 | 42 |
| Exact unique heavy-chain sequences | 83 | 42 | 42 |
| Re-audited RF2 predictions | 84 | 42 | 42 |
| Recorded model-process wall time (seconds) | 15,634.59 | 7,760.41 | 7,874.18 |
| Recorded model-process wall time (hours) | 4.34294 | 2.15567 | 2.18727 |
| All recorded stage attempts | 127 | 63 | 64 |
| Unaccepted attempt time (seconds, included above) | 177.09 | 0 | 177.09 |

The duplicate exact heavy-chain sequence belongs to `pdb_00008qf4`, backbone 0,
in both arms. Therefore the global unique count is not the sum of the arm counts.
Deduplication is by full heavy-chain residue sequence, not sequence similarity,
CDR clustering, predicted structural diversity, or uniqueness of antigen pairs.

Stage totals, including unaccepted attempts:

| Stage | Recorded seconds |
|---|---:|
| RFdiffusion | 12,003.82 |
| ProteinMPNN | 973.44 |
| RF2 | 2,657.33 |

The one unaccepted attempt is in `pdb_00009mnz_unguided`'s RFdiffusion stage.
It contributes cost but contributes no accepted candidate yield.

**Computational-filter pass counts are missing, not zero.** No scientific filter
has been fixed. File/structure/lineage audit success is not a binding or interface
quality filter. The results establish execution yield and recorded cost only.

## Files added and changes made

- `src/evaluate/rfab_compute.py`: standard-library-only, reproducible analysis.
  Reads accepted stage directories from state files, relocates archived WSL paths
  by attempt name, validates structures using the existing generation validators,
  checks sequence lineage, hashes inputs, deduplicates sequences, and reports costs.
- `src/evaluate/test_rfab_compute.py`: four regression/integration tests.
- `results/evaluation/rfab_compute_summary.csv`: three rows (all, hotspot,
  unguided), integer output/attempt counts, exact unique-sequence counts, seconds,
  hours, per-stage seconds, unaccepted attempt costs, and explicit missing-filter status.
- `results/evaluation/rfab_compute_conditions.csv`: 42 target/arm rows with
  accepted counts, within-condition unique sequences, and total/unaccepted seconds.
- `results/evaluation/rfab_compute_attempts.csv`: 127 attempt rows with stage,
  accepted flag, return code when recorded, UTC timestamps, recorded seconds,
  audited-candidate increments, cumulative seconds, and cumulative candidate counts.
- `results/evaluation/rfab_compute_candidates.csv`: 84 candidate identities,
  full-heavy-sequence hashes, repository-relative accepted structure paths, and
  lineage-audit status. Sequence hashes encode ordered three-letter residue names
  joined with hyphens; the independent FASTA check compares one-letter sequences.
- `results/evaluation/rfab_compute_inputs.csv`: SHA-256 manifest of the configuration,
  summary/status inputs, all 42 states, 21 target inputs, and 252 accepted stage PDBs.
- `results/evaluation/rfab_compute_provenance.json`: source commit, script hash,
  input-manifest hash, archive location, and interpretation definitions.
- `results/evaluation/rfab_compute_yield.svg`: standalone cumulative cost/yield plot.
- `A1_COST_ACCOUNTING_REPORT.md`: this change report, methods, results, and validation.

Existing scripts, tracked results, the raw archive, and the supplied brief were
not modified. No generation, prediction, scientific-filter selection, or A2 work
was performed. No commit was created.

## Input/output contract and methods

The archive is `runs/rfab_benchmark_1ce3344e4822/`.
`config.json` supplies 21 unique target IDs and target-file SHA-256 values;
each target has exactly `hotspot` and `unguided` states. Each `state.json` must
identify its target/arm, report completion, and contain exactly the RFdiffusion,
ProteinMPNN, and RF2 stages. Accepted directories and output basenames must agree
with the state. The authoritative candidate key is
`(pdb_id, arm, backbone_index)`; sequence sample 0 is the sole MPNN sample in this
frozen run. Each accepted stage must contain two validated outputs.

Required attempt fields are stage, log path, nonnegative finite `seconds`, and
timezone-aware start/finish timestamps. Accepted attempts additionally require
return code zero. Missing timing or accepted-attempt records fail the analysis;
they are not silently assigned zero. An unaccepted attempt may lack a return code,
which remains blank. All recorded attempts contribute cost even when their outputs
are abandoned. Each accepted RF2 attempt contributes two audited candidates.

Chain H is the generated heavy chain and chain T the cropped antigen. Existing
validators check finite coordinates, H/T chains, CDR remarks, and applicable RF2
confidence scores; light-chain-only NaNs remain inapplicable. Lineage checks verify
antigen sequence preservation, framework sequence preservation during MPNN, and
heavy-chain sequence preservation during RF2. No fixed RF2 residue offset is used;
sequence comparisons use residue order and the existing CDR-label convention.
No native-epitope or experimental-label join is involved in A1.

Recorded seconds are reconciled per condition with
`results/generation_status.csv` and globally with the archive's `summary.json`.
The script sums the recorded duration field rather than subtracting timestamps:
these are distinct measurements in the archive. The graph orders attempts by
recorded completion time and plots cumulative **summed model-process wall time**
against audited candidate count. It is not an elapsed calendar-time timeline.
All data validation completes before output files are written.

The summaries describe the frozen two-design-per-arm run, not a sampling-depth
experiment. No significance tests, uncertainty intervals, or multiple-testing
adjustments apply to these bookkeeping counts. Costs exclude preflight work and
gaps between attempts; they do not measure utilization, allocation, electricity,
monetary cost, or GPU-hours. The plot does not demonstrate biological success or
superiority of either condition.

## Reproduction and validation

From the repository root in WSL/Linux (Python 3, no additional packages needed):

```bash
python3 src/evaluate/rfab_compute.py
python3 -m unittest discover -s src/evaluate -p test_rfab_compute.py -v
```

The full archive is required; Git-tracked summary CSVs alone are insufficient.
Output paths can be overridden with `--out`; `--run` accepts an archive within the
checkout. Team-facing paths remain repository-relative. Rerunning replaces only
the named A1 outputs, leaving archive records untouched.

Validation completed with WSL Python 3.12.3:

- Re-audited all 42 conditions and all 84 candidates across 252 accepted PDBs.
- Checked all 21 target input hashes and unique candidate identities.
- Matched accepted output basenames and exactly one accepted attempt per stage.
- Reconciled every condition's time/counts with the committed generation-status
  table and the total with the archived summary.
- Independently summed the generation-status CSV in PowerShell: 15,634.59 seconds.
- Four tests passed: retry costs without duplicate yield; missing/nonfinite/negative
  timing rejection; missing accepted-attempt rejection; and an independent match of
  all 84 structure-derived sequences to archived FASTA records, including 83 unique
  sequences. The synthetic retry test explicitly selects attempt 2.

No required A1 input remains blocked. A scientific-filter yield remains unavailable
until a filter is independently specified and evaluated; the output records this
limitation explicitly. Recipient-side archive delivery/audit (A0) is not established
by this local analysis.
