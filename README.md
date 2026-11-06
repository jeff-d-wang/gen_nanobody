# gen_nanobody

Targeted de novo nanobody design through backbone generation and sequence optimization.

Pipeline: target structure + epitope hotspots → RFantibody (CDR backbones, framework
fixed) → ProteinMPNN (CDR sequences) → re-prediction without the designed backbone as
template → confidence + geometry + Rosetta interface filters → ranked shortlist.

Proposal and report live in [report/](report/) (`main.tex`, sections in [report/sec/](report/sec/)).

## Layout

| Path | What goes there | Owner |
|---|---|---|
| [src/data/](src/data/) | AVIDa labels, $K_D$ curation, SAbDab2 target set, splits, date filtering | C, Charlie |
| [src/generate/](src/generate/) | RFantibody and RFdiffusion wrappers, hotspot specs, sampling sweeps | A, Aidan |
| [src/evaluate/](src/evaluate/) | ProteinMPNN CDR redesign, structure prediction, Rosetta scoring, metrics | B, Jeffery |
| [third_party/](third_party/) | cloned upstream repos, gitignored: `git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git third_party/ProteinMPNN` | shared |
| [configs/](configs/) | one YAML per run; its hash is the run id | shared |
| [data/raw/](data/raw/) | untouched downloads (gitignored) | shared |
| [data/processed/](data/processed/) | splits, cluster representatives, label and $K_D$ tables (gitignored; commit the split CSVs if small) | shared |
| [runs/](runs/) | per-run outputs and manifest, named by config hash (gitignored) | shared |
| [results/](results/) | the few tables and figures that go into the report (tracked) | shared |

`runs/` is regenerable and large; `results/` is small and cited by the report. That is the
only reason they are separate.

## Rules that aren't negotiable

Three come straight from the proposal, and breaking them invalidates the experiments:

1. The final-scoring predictor is never used during generation or candidate selection.
2. Bootstrap over nanobody and antigen clusters, not individual sequences.
3. Held-out test data is touched only after thresholds are frozen on the development set.

## Status

Midterm due Friday Nov 6. PSC access ends Monday Nov 9, so anything GPU-bound has to run
before then. See [PLAN.md](PLAN.md) for the schedule and owners. Final report Nov 29.
