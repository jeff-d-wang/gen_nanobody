# Revised plan: Nov 4 to Nov 9

Supersedes the schedule in `Nanobody_Team_Execution_Plan.docx`. Roles and the
`C1 → C2 → B1 → B2` critical path are unchanged. What changes is the ordering of the next
five days, plus three deliverables.

A (Aidan): generation and infra. B (Jeffery): sequence design and evaluation.
C (Charlie): data and benchmark.

## Why this changes

1. PSC access ends Monday Nov 9. The old plan scheduled nearly all GPU work *after* the
   midterm: scaling RFantibody across held-out targets, hotspot vs unguided, seed and
   sampling-depth sweeps, structure prediction. None of that can run after Nov 9.
2. We were spending the expiring resource on the wrong tasks. ProteinMPNN runs fine on a
   laptop CPU next week. RFantibody and folding are the only things that need PSC, and the
   old plan marked both "additional if time permits."
3. CDR sequence recovery doesn't test any of our three research questions. It is a useful
   sanity check that ProteinMPNN is wired up correctly and it stays in the report, but it
   can't be the headline preliminary result.

New rule for this week: if it needs a GPU, it happens before Nov 9. Everything else waits.

## One artifact, not three (added Nov 4)

RFantibody bundles base ProteinMPNN with its own CDR-only wrapper (`proteinmpnn -l "H1,H2,H3"`), and its HLT format carries CDR definitions as `REMARK PDBinfo-LABEL` lines derived from Chothia numbering (H1 26-32, H2 52-56, H3 95-102).

So the CDR definition is produced once, by `scripts/util/chothia2HLT.py`, and the
resulting HLT files feed both A's generation and B's recovery scoring. That script is
standalone: argparse, numpy and biotite, with no container, no GPU and no hmmer. The
RFantibody README misnames it `chothia_to_HLT.py`. No CDR columns in `prepped.csv`, and
no ANARCI.

| Overlap | Resolution |
|---|---|
| CDR-only MPNN wrapper | B keeps a standalone one only until A's container is proven, then switches to `proteinmpnn -l H1,H2,H3` |
| CDR definitions | one source: HLT remarks from `chothia2HLT.py` |
| MPNN on *generated* backbones (post-midterm) | A only. It is a stage inside RFantibody, and RQ3 needs both arms from the same stage |
| RF2 vs Boltz | not duplicate. RF2 is inside the pipeline, so it cannot be the excluded final scorer. That must be Boltz or AF3 |

## Decide today (Nov 4), before anything else

| Question | Who | Why it blocks |
|---|---|---|
| Boltz or AF3 for folding? | A | Determines whether we can honor "final-scoring predictor excluded from generation and selection." With only one predictor we cannot, and then it goes in the report as a stated limitation rather than being dropped silently. |
| RFantibody fallback if install fails | A | It is the single point of failure for RQ3. BoltzGen has nanobody support and public weights, and is already in our related work. Name the backup now instead of improvising on Saturday. |
| Date cutoff: is Jan 1 2023 actually enough? | C | The proposal says *after the training cutoffs of the models we use*. AF3 is Sept 2021, and RFdiffusion's weights are separate. One-line check, then move on. |

## Schedule

### Wed Nov 4, today

| Owner | Task | Done when |
|---|---|---|
| C | SAbDab2 download, SD-H + protein-antigen + bound + test-split filter, date filter, one representative per cluster by resolution | `targets.csv` exists and is pushed |
| C | Download targets from SAbDab in Chothia numbering, then `chothia2HLT.py` to get HLT files (epitope definition and 15 Å crop on top) | `hlt/` pushed. This is the blocker for both A and B, so ship it today even if imperfect |
| A | PSC login, `sinfo -s`, remaining GPU-hours and storage, RFantibody install (container route first) | Install works, or fallback decision made |
| B | ProteinMPNN install, PyTorch/CUDA check, CDR fixed-position mask code ready against a dummy input | Can run ProteinMPNN on any PDB end to end |

A and B: don't idle waiting on `hlt/`. Both setup tasks are independent of it.

The install budget for RFantibody is half a day rather than two hours, because it gates a
third of the project. If it is still dead by end of Wednesday, switch to the fallback and
stop.

### Thu Nov 5, GPU day

| Owner | Task | Done when |
|---|---|---|
| A | RFantibody pilot on the 3 smallest targets; verify output format and GPU memory | Backbones exist on disk |
| A | Submit the big batch jobs: held-out targets, hotspot vs unguided, fixed seed set | Jobs queued as *batch*, not interactive |
| B | ProteinMPNN CDR redesign, 8 to 16 seqs/target; CDR recovery vs native | `cdr_recovery.csv` and Figure 3 |
| B | One end-to-end target: prep, ProteinMPNN, fold, then pLDDT / interchain PAE / hotspot-contact recovery | One complete row of real metrics |
| C | Figure 1 (benchmark construction), Table 1 (dataset stats), Figure 2 (cropping) | Pushed to `results/` |

The one end-to-end target is now required rather than optional. It proves the pipeline
exists and forces the metric plumbing we need for every experiment after Nov 9. One
complete case beats 21 targets of recovery statistics.

Queue depth matters more than polish today. A messy backbone set we can analyze in two
weeks is worth more than a clean one we never generated.

### Fri Nov 6, midterm due in the evening

| Owner | Task |
|---|---|
| A | Keep GPU jobs flowing; report allocation used and remaining |
| B | Folding on whatever designed sequences exist; recovery stats and preliminary analysis |
| C | Intro / Related Work / Methodology drafts, limitations, remaining-work timeline |
| All | Report integration |

The limitations section must say that compute access ends Nov 9 and what that costs us. A
stated constraint reads better than a Nov 22 ablation plan that can't run.

### Sat Nov 7 to Sun Nov 8, last GPU weekend

| Owner | Task |
|---|---|
| A | Largest remaining generation sweep; record GPU-hours per run and candidate yield |
| B | Fold every candidate worth folding. This is the last chance |
| A | Archive off PSC by Sunday, verified (checklist below) |
| C | AVIDa-hIL6 curation (CPU, no PSC needed) |

### Mon Nov 9, PSC expires

Buffer only. Confirm the archive restores on someone else's machine. Do not plan work here.

## Midterm deliverables (revised)

| # | Deliverable | Owner | Change |
|---|---|---|---|
| 1 | Figure 1, benchmark construction | C | unchanged |
| 2 | Table 1, dataset statistics | C | unchanged |
| 3 | Figure 2, antigen cropping | C | unchanged |
| 4 | Figure 3, CDR sequence recovery | B | unchanged, but framed as a pipeline sanity check |
| 5 | One end-to-end target with real metrics | B | now required, was "if available" |
| 6 | Generation status and GPU-hours used/remaining | A | new |
| 7 | Written sections and remaining-work timeline including the Nov 9 cliff | All, C coordinates | expanded |

## After Nov 9, CPU only

| Owner | Work |
|---|---|
| A | AVIDa binding classification (RQ1), moved from C. Confidence-only vs geometry-only vs energy-only vs combined, bootstrapped over nanobody clusters. Plus run manifests and cost accounting from the archived runs. |
| B | Rosetta interface scoring and geometry checks on archived structures, candidate ranking, threshold calibration on dev data, failure-mode analysis (framework distortion, CDR refolding, wrong pose, clashes). |
| C | Affinity curation (ANDD, SAbDab-nano, sdAb-DB) and affinity-ranking experiments (RQ2), evaluation harness, leakage and data-separation checks. |
| All | Ablations on cached outputs, error analysis, statistics, final report, shortlist. |

RQ3 (generation yield vs GPU-hours) can only be analyzed from what we archive before Nov 9.
Whatever isn't on disk by Sunday is not in the final report.

## Archive checklist, A owns, due Sun Nov 8

Copy off PSC, then verify the copy opens on a different machine:

- [ ] Generated backbones
- [ ] Designed sequences
- [ ] Folded structures
- [ ] Confidence matrices and summaries (pLDDT, PAE, ipTM)
- [ ] Benchmark and evaluation tables
- [ ] Scripts, job files, model versions, run parameters
- [ ] `sacct` output or equivalent GPU-hour accounting (needed for RQ3)

An unverified transfer discovered on Nov 10 is unrecoverable.

## Non-negotiables (from the proposal)

1. The final-scoring predictor is never used during generation or candidate selection, or
   we state in the report that we couldn't honor this.
2. Bootstrap over nanobody and antigen clusters, not individual sequences.
3. Held-out test data is touched only after thresholds are frozen on dev data.
4. In-silico success is not binding. No claim that a generated nanobody binds.
