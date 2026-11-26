# RFantibody installation pilot

For the completed 21-target run, output locations, metric definitions, and next
evaluation steps, read [RFantibody handoff](../../WORK_ARCHIVE_RFantibody.md).

## 21-target benchmark

Use `benchmark_rfantibody.py` for the benchmark; the scripts under `rfantibody_scripts/`
are retained as reference material and are not the active runner.

```bash
python3 src/generate/benchmark_rfantibody.py --prepare-only
python3 src/generate/benchmark_rfantibody.py --pilot-only
python3 src/generate/benchmark_rfantibody.py
```

The first command validates all 21 epitope mappings and writes target-only
PDBs to `data/processed/rfab_inputs/`. The second runs the smallest target's
hotspot condition through all three models. The third runs both conditions
for all targets, reusing validated completed stages from the pilot.

The frozen small-sweep settings are two deterministic backbones (design seeds
0 and 1) per target/condition, H1/H2/H3 lengths 7/6/8, 50 diffusion steps, one
MPNN sequence per backbone at temperature 0.1, and one RF2 recycle with seed
37. These deliberately small settings test the workflow, not scientific
significance. Both conditions share the same antigen crop and framework.

The runner uses the sibling installation's environment and CLI, including its
local weights. A configuration hash identifies the run directory. It locks
against concurrent invocation, validates completed stages before reusing
them, and retries failed stages only with `--retry-failed`. Each retry uses
a fresh attempt directory and retains earlier outputs and timings. Model
commands have a one-hour per-stage limit; timeouts terminate their process
group. The first target must pass before a new sweep proceeds.

`runs/rfab_benchmark_<hash>/` contains inputs, the framework, scripts, weight
hashes, upstream revision, per-attempt logs, per-condition state files,
structures, designed sequences in FASTA, confidence metrics, and status
summaries. `results/generation_status.csv` and
`results/rfab_candidate_metrics.csv` provide small report-facing tables.
Time totals include failed attempts and are process wall time, not measured
GPU utilization. Partial files remain in their attempt directories; they are
not automatically passed downstream. Raw outputs remain gitignored.

Input preparation uses four closest-contact hotspots with an 8 Å CA spacing
preference, relaxing that spacing to fill four slots when necessary. The
manifest records this relaxation. The MSE at 8rwf T85 is converted to MET for
model parsing, preserving coordinates and numbering (SE becomes SD/S).
Canonical reference files are unchanged. All epitope residues must map, all
chosen hotspots must have CA coordinates, and each antigen's residue numbers
must match the mapping table. Known structural QC flags are retained.

These are held-out benchmark inputs. No thresholds are calibrated here and
RF2 confidence is not an independent final evaluation or evidence of binding.

CPU checks, after input preparation:

```bash
python3 -m unittest discover -s src/generate -p test_benchmark_rfantibody.py -v
```

## Standalone installation check

Run from WSL, at the `gen_nanobody` repository root:

```bash
python3 src/generate/pilot_rfantibody.py
```

The script uses the existing `../RFantibody/.venv` explicitly, including its
Python interpreter for subprocesses. It does not install packages. Override
the installation with `--rfantibody-root /path/to/RFantibody` if needed.

Each invocation creates a new timestamped directory under `runs/`. To choose
a location, use `--output-dir runs/my_pilot`; the directory must be empty.

Outputs:

- `summary.json`: overall pass/fail, errors, commands, stage timings, upstream
  commit, input hashes, residue counts, and RF2 confidence scores.
- `pilot.log`: full combined stdout/stderr from all stages.
- `inputs/`: lysozyme antigen and independent NbBCII10 framework.
- `backbones/`: one RFdiffusion backbone plus its upstream metadata.
- `sequences/`: one ProteinMPNN-designed structure.
- `predictions/`: RF2 structures, including the best prediction with scores.

The target is chain T from the committed 1MEL test fixture (hen lysozyme).
The reference nanobody is not supplied to generation. This unguided test uses
one backbone, fixed H1/H2/H3 lengths of 7/6/8, 50 diffusion steps, one sequence
at temperature 0.1, and one RF2 recycle. Generation and MPNN use upstream
deterministic mode; RF2 receives seed 37. This is a small installation test,
not an experimental setting or a benchmark result.

Success requires all commands to exit successfully, exactly one backbone and
sequence structure, exactly one best RF2 prediction, H/T chains and all three
CDR labels, finite coordinates, and finite RF2 `interaction_pae` and
`pred_lddt` scores. Low confidence does not fail the installation check.
Light-chain loop RMSDs are not applicable to nanobodies: upstream RF2 emits
`nan` for these three scores, which the summary records as JSON `null`.
All other nonfinite scores fail validation.
The process exits nonzero on failure and preserves partial outputs and logs.
An interrupted or failed run should be retried in a fresh output directory.
To recheck saved outputs without rerunning the models, use
`python3 src/generate/pilot_rfantibody.py --validate-only --output-dir runs/my_pilot`.
This retains the previous validation status in the summary and appends to the log.

This does not validate binding, assess benchmark performance, or replace
independent final scoring. Raw run outputs are gitignored by repository policy.
