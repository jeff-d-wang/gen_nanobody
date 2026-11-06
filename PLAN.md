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
   laptop CPU next week: measured Nov 6, all 21 targets at 16 sequences each took about
   12 seconds per target, roughly 4 minutes in total, and it never touches the GPU
   (ProteinMPNN only tests `torch.cuda`, so Apple MPS goes unused anyway). RFantibody and
   folding are the only things that need PSC, and the old plan marked both "additional if
   time permits."
3. CDR sequence recovery doesn't test any of our three research questions. It is a useful
   sanity check that ProteinMPNN is wired up correctly and it stays in the report, but it
   can't be the headline preliminary result.

New rule for this week: if it needs a GPU, it happens before Nov 9. Everything else waits.

## One artifact, not three (added Nov 4)

RFantibody bundles base ProteinMPNN with its own CDR-only wrapper (`proteinmpnn -l "H1,H2,H3"`), and its HLT format carries CDR definitions as `REMARK PDBinfo-LABEL` lines derived from Chothia numbering (H1 26-32, H2 52-56, H3 95-102).

So the CDR definition is produced once and the resulting HLT files feed both A's
generation and B's recovery scoring. One artifact still holds, but the producer turned
out to be ours rather than RFantibody's. See the next section. (The RFantibody README
also misnames the script `chothia_to_HLT.py`; the real path is
`scripts/util/chothia2HLT.py`.)

| Overlap | Resolution |
|---|---|
| CDR-only MPNN wrapper | B keeps a standalone one only until A's container is proven, then switches to `proteinmpnn -l H1,H2,H3` |
| CDR definitions | one source: `data/processed/hlt/`, written by `src/evaluate/make_hlt.py`. Not `chothia2HLT.py`, because C's structures are IMGT-numbered |
| MPNN on *generated* backbones (post-midterm) | A only. It is a stage inside RFantibody, and RQ3 needs both arms from the same stage |
| RF2 vs Boltz | not duplicate. RF2 is inside the pipeline, so it cannot be the excluded final scorer. That must be Boltz or AF3 |

## A: do not run chothia2HLT.py on C's structures (added Nov 6)

C's `prepped_v2/` and `prepped_mpnn/` are IMGT-numbered, not Chothia. All 21 structures
put the conserved anchors at IMGT 23, 41 and 104, the heavy chain runs to 128, and 8qf4
carries CDR3 insertions at 111A to 112C. Chothia stops at 113 and inserts at 100A to
100K instead. The numbering came in with the deposited mmCIFs; C's prep script uses
gemmi and preserves whatever it is given, so this is nobody's bug.

`chothia2HLT.py` assumes Chothia and never checks the numbering of its input. Pointed at
these files it applies H1 26-32, H2 52-56 and H3 95-102 to IMGT positions and writes
wrong CDR labels with no error, and nothing downstream would catch it. The two
definitions really do disagree: 26.9 designed positions per target under Chothia against
31.9 under IMGT.

Use `src/evaluate/make_hlt.py` instead. ANARCI renumbers each heavy chain to Chothia and
the same `REMARK PDBinfo-LABEL` lines are appended, in a format byte-compatible with
`chothia2HLT.py`. Coordinates are copied unchanged, because C's `prepped_mpnn/` files
already fix a real ProteinMPNN bug: the parser misorders insertion codes and read 8qf4's
CDRH3 as `GKSI` instead of `ISKG`. Use `prepped_mpnn/`, never `prepped_v2/`, for anything
going through the MPNN parser.

```
conda create -y -n nbnum -c bioconda -c conda-forge python=3.10 anarci hmmer biopython
conda run -n nbnum python src/evaluate/make_hlt.py \
    --pdbs data/processed/prepped_mpnn --out data/processed/hlt \
    --targets data/processed/targets.csv
```

ANARCI needs that separate env. Installing it into a conda base env aborts on conda's
own dependencies and changes nothing. No other stage needs it.

Checked on all 21 targets: CDR3 is flanked by the conserved Cys92 and Trp103 on every
one, Chothia H3 falls entirely inside C's independently derived IMGT CDRH3 on every one,
and each output file differs from its input only by the added remarks. 8yj8, 9h77 and
9q6n carry Arg at Chothia 103 rather than the near-invariant Trp, which is real
sequence variation and not a numbering error.

`data/processed/hlt/` is gitignored on purpose. Regenerate it rather than committing it,
so the CDR definition cannot fork into a stale copy.

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
