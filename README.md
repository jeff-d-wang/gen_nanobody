# gen_nanobody

Targeted de novo nanobody design through backbone generation and sequence optimization.

Pipeline: target structure + epitope hotspots → RFantibody (CDR backbones, framework fixed)
→ ProteinMPNN (CDR sequences) → re-prediction **without** the designed backbone as template
→ confidence + geometry + Rosetta interface filters → ranked shortlist.

Proposal and report live in [report/](report/) (`main.tex`, sections in [report/sec/](report/sec/)).

## Layout

| Path | What goes there | Owner |
|---|---|---|
| [src/data/](src/data/) | AVIDa labels, $K_D$ curation, SAbDab2 target set, splits, date filtering | Data |
| [src/generate/](src/generate/) | RFantibody / ProteinMPNN wrappers, hotspot specs, sampling sweeps | Generation |
| [src/evaluate/](src/evaluate/) | structure prediction, Rosetta scoring, metrics, threshold calibration | Evaluation |
| [configs/](configs/) | one YAML per run; its hash is the run id | — |
| [data/raw/](data/raw/) | untouched downloads (gitignored) | — |
| [data/processed/](data/processed/) | splits, cluster representatives, label/$K_D$ tables (gitignored; commit the split CSVs if small) | — |
| [runs/](runs/) | per-run outputs + manifest, named by config hash (gitignored) | — |
| [results/](results/) | the few tables/figures that go into the report (tracked) | — |

`runs/` is regenerable and large; `results/` is small and cited by the report. That's the only
reason they're separate.

## Rules that aren't negotiable

Three come straight from the proposal, and breaking them invalidates the experiments:

1. The final-scoring predictor is **never** used during generation or candidate selection.
2. Bootstrap over nanobody/antigen **clusters**, not individual sequences.
3. Held-out test data is touched only after thresholds are frozen on the development set.

## Status

Midterm shipped Nov 1. Current week: full generation sweep across held-out targets and seeds.
Final report Nov 29.
