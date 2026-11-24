# Revised plan: Nov 4 to Nov 9

> **Start at the Nov 24 addendum at the end of this file.** It is the current state:
> the EC2 box, the runbook, what is done and what is next. The Nov 16 addendum before
> it still defines the *work* (B0–B5) but its compute section is superseded. PSC access
> is gone; GPU work moved to AWS and B owns it. Everything in this first section is kept
> only as the record of what was decided and why.

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

---

# Nov 16 addendum: no PSC, AWS instead. B owns the compute.

Supersedes every compute assumption above. The Nov 4–9 schedule is history; keep it
for the record. PSC is gone and the archive checklist was never verified against a
second machine, so that is now a prerequisite, not a formality.

## What changed since Nov 9

| Fact | Consequence |
|---|---|
| PSC access ended | `src/evaluate/psc_boltz.sbatch` is dead as a submission path. Keep it as the recorded run config; the Boltz flags, GPU preflight and native-first ordering all carry over |
| A finished a 21-target RFantibody run on a local RTX 3060 Ti (`origin/RFantibody`, Nov 8) | 84 candidates exist. Generation is no longer the bottleneck; evaluation is |
| A's 84 structures are gitignored | Scripts and 2 summary CSVs are on the branch. The structures must arrive by file transfer. Nothing B owns can start on real candidates until they do |
| RF2 ran *inside* A's pipeline | RF2 can never be the independent final scorer. Non-negotiable #1 is now entirely B's: fold from sequence only, no template, no designed pose |
| B has $250 AWS credit, no PSC | Boltz moves to one EC2 GPU instance. See the compute section |

## B's work, in dependency order

### B0. Unblock: get A's run directory (do this first, it gates B2–B5)

Ask A for `runs/rfab_benchmark_1ce3344e4822/` as an archive, through shared storage,
not Git. Then, per A's handoff:

```bash
python3 src/generate/audit_rfantibody.py runs/rfab_benchmark_1ce3344e4822
# expect passed: true, target_conditions: 42, validated_candidates: 84, errors: []
```

That audit doubles as the Nov 8 archive verification the old checklist deferred.
A's CSVs carry absolute WSL paths (`/mnt/c/Users/adnja/...`); rewrite the repo prefix,
keep the rest of the layout. Candidate identity is `(pdb_id, arm, backbone_index)` and
must not depend on any absolute path.

### B1. Finish midterm deliverable 5: one end-to-end target, real metrics

Only item still outstanding from the midterm list. Inputs are already built:
`runs/boltz_in/pdb_00008qf4_H_T/` holds `native.yaml` plus 4 designs (samples 5, 16, 2, 9
by `global_score`). Needs ~5 folds, then:

```bash
python3 src/evaluate/fold_metrics.py validate --target pdb_00008qf4   # run this first
python3 src/evaluate/fold_metrics.py report  --target pdb_00008qf4 \
    --boltz runs/boltz_out --inputs runs/boltz_in
```

`validate` reproduces C's 5 Å epitope set from the HLT file before any design number is
believed. Both scripts are still untracked — commit them with the first real results.

### B2. Independent final scoring of A's 84 candidates (the non-negotiable)

Fold each designed heavy chain plus the cropped antigen **from sequence alone**. Boltz
never sees the generated backbone, the RF2 prediction, or the docked pose as a template.
Sequences come from `RUN/designed_sequences.fasta`; antigen sequences from the same
`parsed.jsonl` route `make_fold_inputs.py` already uses, so the antigen is byte-identical
across every fold and contact recovery stays a set comparison on chain-T indices.

**Inputs built Nov 16: 109 jobs staged in `runs/boltz_in/`** — 21 native controls, 84
RFantibody candidates, and 8qf4's 4 MPNN redesigns, merged into one manifest per target
so a single `fold_metrics.py report` covers all three kinds under the same scorer. Run
`src/evaluate/run_boltz.sh` from the input root. Nothing is left to prepare but the GPU.

Checked while building, all 84 candidates: every `designed_sequences.fasta` record matches
the chain-H sequence of its accepted ProteinMPNN output (so candidate identity is sound,
and the accepted attempt is read from `state.json`, never guessed), and stripping each
candidate's own CDR positions leaves one shared 95-aa framework string — the fixed-framework
invariant, confirmed without any alignment. 95 + 21 designed positions = 116 aa.

B's own MPNN redesigns for the other 20 targets are deliberately not staged. Add them
when a question needs them, not before.

Verified Nov 16, all 21 targets: the chain-T sequence in A's `RUN/inputs/targets/*_T.pdb`
is byte-identical to `seq_chain_T` in B's `runs/mpnn_seed37/parsed.jsonl`, numbered 1..N
contiguously. So `fold_metrics.py`'s core assumption survives onto A's candidates
unchanged — antigen residue *i* is antigen residue *i*, contact recovery stays a set
comparison on T indices, no structural alignment anywhere — and C's epitope sets and A's
hotspot lists in T numbering apply directly. This is the invariant that would otherwise
have forced a rewrite of the metrics code. Re-run the check if either side regenerates.

Largest antigen is 300 aa (9dc3), so the biggest complex is ~416 residues. That is what
sets the GPU memory floor below, not the median target.

MSA cost note: the 21 antigens are identical across arms and designs, so fetch 21 antigen
MSAs once and reuse them. Do not send 105 queries to the public MSA server. Nanobody
chains go single-sequence — the server's hits for a designed VHH are framework noise.

### B3. Interface and geometry evaluation on both structure sets

Per candidate, on the generated complex, A's RF2 prediction, and B's independent Boltz
prediction: hotspot contacts, full-epitope contacts, retained-contact fraction (state the
cutoff and the denominator, and keep the three counts separate), steric clashes, buried
area, interface geometry.

A's own numbers already flag where the interesting failures are: that first row has
`pred_lddt` 0.87 with `target_aligned_antibody_rmsd` 36.3 Å — a confident fold that left
the epitope entirely. Confidence-vs-placement disagreement is the headline failure mode,
not a footnote.

Rosetta `InterfaceAnalyzer` is CPU-only and stays on the laptop. If the academic license
or the build costs more than an afternoon, drop it and report contacts, clashes and buried
area instead; it is a secondary scorer, not a dependency.

### B4. Candidate ranking and threshold calibration

Freeze thresholds on dev data before touching these 21 held-out targets. Until a dev set
exists (C's curation), report descriptive distributions only and say so. Ranking must not
use `pred_lddt` alone.

### B5. Hotspot vs unguided, and failure-mode analysis

Summarize within each target/condition first, then compare as 21 **paired target-level**
outcomes. Never 84 independent observations. Bootstrap over nanobody/antigen clusters.
Count sequence duplicates among the 84 and state the clustering rule.

> **Duplicate count, answered Nov 24: 83 distinct sequences of 84.** The single
> duplicate is on 8qf4, where one hotspot candidate's designed sequence is identical
> to one unguided candidate's. 8qf4 is an included target, so state that its arm
> comparison is partly a self-comparison. Nothing else is duplicated, in either
> direction, on any other target.
>
> This also rules out duplicates as the explanation for the three exact arm ties on
> the included targets (9q6n 0.086 vs 0.086, 9rye and 9ty2 both 0.000 vs 0.000).
> Those are genuine, and two of them are both arms failing to contact the epitope
> at all. The clustering rule for the bootstrap is still owed.

Failure taxonomy: framework distortion, CDR refolding, wrong pose, clashes, and the
confident-but-displaced class above.

## Compute: one EC2 GPU instance, B owns it

> Superseded by the Nov 24 addendum, which records the instance as actually built and
> the measured costs. The quota reasoning below is still the explanation of why; the
> approved quota turned out to be 8 vCPU on-demand and 4 spot, not the 16 requested.

Feasible, with room to spare. The only GPU work left in B's queue is Boltz folding of
roughly 105–200 complexes of about 190 residues each.

| Item | Estimate |
|---|---|
| Per fold, 1 diffusion sample, 3 recycles, ~190 residues | 2–5 min |
| 105 folds + retries | ~6–12 GPU-hours |
| Instance | `g5.xlarge`, A10G 24 GB, ≈$1.0/hr on-demand in us-east-1 (confirm in console) |
| Folding cost | ~$10–15. Under $5 on spot |
| Whole-project ceiling if B2+B3 are rerun twice | well under $50 of $250 |

Credit is not the constraint. Three things actually are:

1. **GPU quota.** A new AWS account has an on-demand G-instance vCPU quota of 0.
   In Service Quotas → Amazon EC2, raise **Running On-Demand G and VT instances**
   (code `L-DB2E81BA`), measured in **vCPUs**, not instances. `g5.xlarge` is 4 vCPUs,
   so 4 permits exactly one at a time; **ask for 16** to allow a parallel batch or a
   bigger instance later, since a second request costs another day and a quota is a
   cap, not a reservation. Spot is a separate quota (**All G and VT Spot Instance
   Requests**, `L-3819A6DF`) — request it in the same sitting. Both are **per region**,
   so file them in the region you will actually launch in. Do this **today**: approval
   takes hours to a couple of days and nothing else can start until it lands.
2. **Credit terms.** Check Billing → Credits for expiry and service restrictions
   before planning around the $250.
3. **A stopped instance still bills EBS.** Keep the root volume at 100 GB (~$8/mo),
   `stop` the instance between sessions, `terminate` and delete the volume when B2 is
   done. Set an AWS Budget alert at $50. A GPU instance left running over a weekend is
   the only realistic way to lose this credit.

Setup, in order: quota request → Deep Learning AMI (Ubuntu 22.04, PyTorch) so the CUDA
driver is already correct → `pip install boltz` → keep `BOLTZ_CACHE` on the root volume →
run the preflight GPU check lifted from `psc_boltz.sbatch` before spending GPU time →
batch folds under `nohup`, native first.

Claude can drive the instance over SSH from this repo (jobs backgrounded, polled), so the
workflow is the same as the sbatch one minus the scheduler: security group open on port 22
to B's IP only, or SSM Session Manager if no open port is preferred. Keys and credentials
stay on B's machine.

Done Nov 16: `src/evaluate/run_boltz.sh` is the plain-bash successor to
`psc_boltz.sbatch`. The Slurm directives were the only part that died. The GPU preflight,
`--output_format pdb`, native-first ordering and the resume check all carried over, and
the resume check matters more on EC2 than it did on Slurm — it is what makes a dropped
ssh session, an overnight interruption or a spot reclaim survivable. The sbatch file is
kept, marked superseded, as the record of the PSC config.

Open: populating the cached antigen MSAs. `--msa cached` and the runner's missing-a3m
abort are in place, but Boltz's on-disk a3m layout is not guessed in code — fold the
first native with `MSA=server`, look at what Boltz wrote, then wire the 21 reused files.
Until then `MSA=server` works and costs 105 server queries instead of 21.

## Still true, and still binding

The four non-negotiables are unchanged. B2 is what makes #1 survivable, and B4 is what
keeps #3 honest. In-silico success is still not binding.

---

# Nov 24 addendum: the EC2 box is live. Runbook, state, and one real problem.

Supersedes the compute section of the Nov 16 addendum. The *work* items (B2–B5) are
unchanged; what changes is that the instance exists, the MSA open item is closed, and
the native controls have raised a question that has to be answered before any design
number is interpreted.

## State of play

| Item | State |
|---|---|
| B0 unblock | Done. `runs/rfab_benchmark_1ce3344e4822/` is local, `audit.json` reads `passed: true, target_conditions: 42, validated_candidates: 84, errors: []`. That is also the Nov 8 archive verification the old checklist deferred |
| 109 fold inputs | Staged and synced to the instance |
| MSA open item | Closed. Boltz emits a merged CSV, not a reusable a3m. See below |
| GPU environment | Verified working. Native 8qf4 folded: iptm 0.940, complex_plddt 0.968 |
| 21 native controls | Folded (one per target), MSAs harvested |
| 88 design folds | Not started. This is the next action |
| B3–B5 | Not started, CPU-only, no AWS needed |

## The instance

| | |
|---|---|
| Region | us-east-2 (Ohio). Quota is per region, so it must be approved here, not us-east-1 |
| Host | `ec2-18-116-80-200.us-east-2.compute.amazonaws.com` |
| Connect | `ssh -i ~/.ssh/aws.pem ubuntu@<host>` (key is `0400`) |
| Type | `g5.xlarge`, A10G 24 GB, `sm_86`, 4 vCPU |
| Quota approved | 8 vCPU on-demand G/VT (`L-DB2E81BA`) = 2× g5.xlarge; 4 vCPU spot (`L-3819A6DF`) = 1× g5.xlarge |
| AMI | Ubuntu 24.04.5 LTS (Canonical, owner `099720109477`), not the DLAMI |
| Python / torch / boltz | 3.12.3 / 2.15.0+cu132 / 2.2.1 with the `[cuda]` extra |
| Driver | 595.91.07, `nvidia-utils-595-server` |
| Root volume | 43 GB usable. Tight. Boltz cache alone is 7.6 GB |
| Venv | `~/bz`, `BOLTZ_CACHE=~/.boltz` |
| Remote layout | Mirrors the repo: `~/gen_nanobody/runs/boltz_in`, `~/gen_nanobody/src/evaluate` |

Boltz stays the excluded final scorer: sequence only, no template, no designed pose.
Nothing about moving to EC2 changes non-negotiable #1.

## Setup, as actually performed (replaces the Nov 16 "Setup, in order" line)

```bash
# 1. instance: Ubuntu 24.04 LTS from Quick Start (NOT the search box, which only
#    offers a paid SQL Server bundle), g5.xlarge, 45 GB gp3, SSH from your IP only.
# 2. driver + headers
sudo apt-get update
sudo apt-get install -y build-essential python3-venv python3-dev ubuntu-drivers-common
sudo ubuntu-drivers install --gpgpu
sudo reboot
# 3. after reboot
sudo apt-get install -y nvidia-utils-595-server     # for nvidia-smi; match the driver
nvidia-smi                                          # expect A10G, 23028MiB
# 4. boltz
python3 -m venv ~/bz && source ~/bz/bin/activate
pip install --upgrade pip && pip install boltz[cuda] -U
echo 'export BOLTZ_CACHE=$HOME/.boltz' >> ~/.bashrc
# 5. ship inputs (from the repo root, on the laptop)
rsync -az -e "ssh -i ~/.ssh/aws.pem" runs/boltz_in ubuntu@<host>:~/gen_nanobody/runs/
rsync -az -e "ssh -i ~/.ssh/aws.pem" src/evaluate ubuntu@<host>:~/gen_nanobody/src/
```

`python3-dev` is not optional. Boltz uses Triton, which JIT-compiles a CUDA shim
with gcc on first GPU use and needs `Python.h`. Without it every fold dies after the
MSA fetch with a `CalledProcessError` from gcc, which reads like a CUDA problem and is
not one. `build-essential` and `python3-venv` do not pull the headers in. This is the
EC2 analogue of the PSC `sm_70` trap: a late failure once per fold.

## MSA: the open item is closed, and the answer is a CSV

With `--use_msa_server` Boltz writes the merged MSA to

```
<out_dir>/boltz_results_<job>/msa/<job>_<chain_idx>.csv     # header: key,sequence
<out_dir>/boltz_results_<job>/msa/<job>_unpaired_tmp_env/   # raw per-database a3m
```

`<job>_1.csv` is chain index 1, the antigen (YAMLs always list nanobody A then antigen B).
The a3m files in `_unpaired_tmp_env/` are the separate MMseqs2 databases (uniref,
bfd/mgnify/metaeuk/smag), unmerged, and concatenating them is not the same input, so do
not try to rebuild an a3m. `boltz/main.py` dispatches on file suffix and parses `.a3m` or
`.csv` equally, so the harvested CSV is reused directly.

Consequences, already applied:

- `make_fold_inputs.py:msa_path()` emits `<pdb_id>.csv`, not `.a3m`. The `--demo`
  assertion was updated with it and passes.
- All 109 staged YAMLs were switched, each by exactly one line, verified per file
  (line counts equal, one differing line, differing only in the extension). 21 unique
  MSA references, all `.csv`.
- The harvest recipe is `/tmp/harvest.sh` on the instance: for each target, strip the
  cached-MSA line from `native.yaml`, fold with `--use_msa_server`, copy
  `boltz_results_native/msa/native_1.csv` to `runs/boltz_in/msa/<pdb_id>.csv`.
  21 server queries instead of 109, and the native controls come out of the same
  folds, so no GPU time is spent only on MSAs.

The harvested `msa/` directory is worth keeping: it makes a re-run reproducible without
re-querying the public server. Pull it back alongside the outputs.

## PROBLEM: most native controls do not reproduce their known interface

This is the finding of the session and it outranks everything else here.

PLAN.md has said since the PSC days that natives run first because "if the predictor
cannot reproduce a known interface the design numbers are not interpretable anyway."
That gate now fails or is marginal for 12 of the 21 targets (9 clear 0.7, 4 sit
between 0.5 and 0.7, 8 fall below 0.5). Native iptm, 1 diffusion sample, 3 recycling
steps, antigen MSA from the server, nanobody single-sequence:

| Target | native iptm | complex_plddt | MSA rows |
|---|---|---|---|
| pdb_00008qf4 | 0.940 | 0.968 | 4903 |
| pdb_00009rye | 0.936 | 0.904 | 15851 |
| pdb_00009e7l | 0.918 | 0.594 | 1717 |
| pdb_00009g3e | 0.915 | 0.727 | 2589 |
| pdb_00009q6n | 0.902 | 0.818 | 1284 |
| pdb_00009min | 0.876 | 0.947 | 5111 |
| pdb_00009ty2 | 0.873 | 0.932 | 2312 |
| pdb_00009b71 | 0.765 | 0.927 | 7565 |
| pdb_00009h77 | 0.762 | 0.893 | 14196 |
| pdb_00009mnz | 0.666 | 0.750 | 1207 |
| pdb_00008idi | 0.646 | 0.834 | 204 |
| pdb_00009yc3 | 0.581 | 0.847 | 42 |
| pdb_00009g2g | 0.523 | 0.799 | 10597 |
| pdb_00009rxi | 0.446 | 0.687 | 12060 |
| pdb_00009my8 | 0.370 | 0.839 | 106 |
| pdb_00008rwf | 0.325 | 0.895 | 165 |
| pdb_00009ljj | 0.294 | 0.663 | 1749 |
| pdb_00008yj8 | 0.273 | 0.765 | 10804 |
| pdb_00009rmk | 0.244 | 0.656 | 732 |
| pdb_00009dc3 | 0.230 | 0.488 | 615 |
| pdb_00009npt | 0.174 | 0.724 | 66 |

Complete, all 21 targets. Harvest: 20 folded, 0 failed, 1 already done. 9 clear iptm
0.7; 8 fall below 0.5.

> **This 1-sample table is superseded for all decision-making** by the 5-sample table
> further down. It is kept because the 1-vs-5 comparison is itself a result about
> single-sample confidence, and because `runs/boltz_out/` still holds these folds.

MSA depth does not explain it. The shallowest MSA in the set, 9yc3 at 42 rows,
reaches 0.581 — better than 8yj8 at 10804 rows (0.273) and 9rxi at 12060 rows (0.446).
9dc3, our largest antigen at 300 aa, reaches only 0.230. Whatever drives this, it is a
property of the complex, not of how much coevolutionary signal we fed it. Do not "fix"
the weak targets by chasing deeper MSAs.

complex_plddt is a bad proxy for interface quality, and the natives prove it.
8rwf pairs iptm 0.325 with plddt 0.895; 9npt pairs iptm 0.174 with plddt 0.724; and in
the other direction 9e7l pairs iptm 0.918 with plddt 0.594. Monomer confidence is high
while interface confidence is low, and vice versa. This is the same
confidence-vs-placement disagreement A already saw in RF2 (`pred_lddt` 0.87 with
`target_aligned_antibody_rmsd` 36.3 Å), now reproduced by an independent predictor on
*native* complexes where the right answer is known. It is direct empirical support for
the existing rule that ranking must never use `pred_lddt` alone — B4's constraint is no
longer just a precaution.
What must be decided before looking at a single design number. Native control
quality is a property of the *target and the predictor*, not of any candidate, so
stratifying on it leaks nothing about the designs. But it is still a threshold, and
non-negotiable #3 says thresholds are frozen before held-out data is touched. So:

1. Pick and write down a native-iptm floor now, before the 88 design folds are read.
2. Report the designs split into interpretable and uninterpretable target sets, never
   pooled. A design "beating" its native on a target whose native sits at 0.17 is noise.
3. Do not quietly drop the failing targets. Six of 21 clearing a reasonable bar is a
   result about Boltz-as-independent-scorer on nanobody complexes, and it belongs in the
   report either way.

### ANSWERED Nov 24: one sample was unlucky. Do not set a floor on 1-sample data.

Re-folded the bottom 5 natives at `--diffusion_samples 5`, recycling held at 3 so only
one variable moved, plus 8qf4 as a positive control (max-of-5 beats a single draw even
from an unchanged distribution, so a control is needed to size that artifact).

| Target | 1 sample | model_0 of 5 | min | max | spread |
|---|---|---|---|---|---|
| pdb_00009npt | 0.174 | 0.377 | 0.226 | 0.508 | 0.282 |
| pdb_00009dc3 | 0.230 | 0.458 | 0.183 | 0.458 | 0.275 |
| pdb_00009rmk | 0.244 | 0.454 | 0.150 | 0.454 | 0.304 |
| pdb_00008yj8 | 0.273 | 0.592 | 0.341 | 0.592 | 0.251 |
| pdb_00009ljj | 0.294 | 0.830 | 0.364 | 0.830 | 0.467 |
| pdb_00008qf4 (control) | 0.940 | 0.950 | 0.926 | 0.950 | 0.023 |

The control's 5 samples span 0.023 and model_0 gains +0.010 over the single draw, so
sample-selection is worth about one hundredth of an iptm point on a target Boltz handles.
The weak targets span 0.251 to 0.467, 11× to 20× wider, and gain +0.20 to +0.54.
That is far too large to be the selection artifact the control measures.

So the weak natives are not "Boltz cannot do this complex." They are multi-modal: one
diffusion sample is noise. 9ljj moves from the worst band to clearing 0.7 outright.
Three others roughly double and still land below 0.5, so the stratification does not
vanish, but it cannot be read off 1-sample numbers.

Five samples are nearly free. 55–92 s per target for 5, against ~75 s for 1: the
trunk pass (MSA processing, pair representation) dominates and runs once regardless.
There is no reason to run 1 sample anywhere in this project.

Consequences, binding:

1. The 21-target 1-sample table above is **not a basis for any threshold**. Re-run all
   21 natives at 5 samples before the floor is set (~30 min, well under $1).
2. Run the 88 design folds at 5 samples too. At 1 sample every design iptm carries a
   ±0.3 uncertainty, which makes B5's paired target-level comparisons meaningless.
3. Use `model_0` as the per-job value. That is Boltz's own composite ranking and is
   usually, but not always, the iptm argmax — on 9npt, model_0 is 0.377 while model_2
   reaches 0.508. Taking the iptm argmax across samples would mean selecting on the very
   statistic being reported; `model_0` is pre-specified and avoids that. Record n=5 for
   every job so the selection pressure is identical everywhere.
4. `fold_metrics.py` needs no change: its `confidence_*model_0.json` and
   `*model_0.pdb` globs each match exactly one file in a 5-sample directory, so the
   `one()` helper is satisfied and it reads Boltz's top-ranked sample by construction.

### The authoritative native table: all 21 at 5 samples

`runs/boltz_out_native5/`, 5 diffusion samples, recycling 3, cached MSAs, zero server
queries. 15 folded Nov 24 + the 6 from the experiment above reused at identical settings.
0 failures. Sorted by `min`, which is the criterion recommended below.

| Target | model_0 | min | median | max | spread |
|---|---|---|---|---|---|
| pdb_00008qf4 | 0.950 | 0.926 | 0.939 | 0.950 | 0.023 |
| pdb_00009g3e | 0.904 | 0.904 | 0.941 | 0.949 | 0.046 |
| pdb_00009e7l | 0.947 | 0.903 | 0.927 | 0.963 | 0.060 |
| pdb_00009rye | 0.920 | 0.894 | 0.920 | 0.930 | 0.035 |
| pdb_00009min | 0.892 | 0.883 | 0.884 | 0.893 | 0.010 |
| pdb_00009ty2 | 0.886 | 0.869 | 0.880 | 0.886 | 0.017 |
| pdb_00009mnz | 0.958 | 0.867 | 0.954 | 0.962 | 0.095 |
| pdb_00009q6n | 0.922 | 0.838 | 0.869 | 0.922 | 0.084 |
| pdb_00009b71 | 0.858 | 0.767 | 0.788 | 0.858 | 0.092 |
| pdb_00009g2g | 0.566 | 0.481 | 0.530 | 0.566 | 0.086 |
| pdb_00009ljj | 0.830 | 0.364 | 0.649 | 0.830 | 0.467 |
| pdb_00008yj8 | 0.592 | 0.341 | 0.472 | 0.592 | 0.251 |
| pdb_00009h77 | 0.726 | 0.333 | 0.608 | 0.726 | 0.394 |
| pdb_00009rxi | 0.616 | 0.296 | 0.600 | 0.658 | 0.362 |
| pdb_00009my8 | 0.655 | 0.260 | 0.483 | 0.655 | 0.395 |
| pdb_00009yc3 | 0.613 | 0.255 | 0.420 | 0.613 | 0.357 |
| pdb_00008idi | 0.656 | 0.255 | 0.348 | 0.656 | 0.402 |
| pdb_00008rwf | 0.722 | 0.236 | 0.403 | 0.722 | 0.486 |
| pdb_00009npt | 0.377 | 0.226 | 0.300 | 0.508 | 0.282 |
| pdb_00009dc3 | 0.458 | 0.183 | 0.276 | 0.458 | 0.275 |
| pdb_00009rmk | 0.454 | 0.150 | 0.312 | 0.454 | 0.304 |

**The spread is strictly bimodal.** 10 targets at spread ≤ 0.10, nothing between 0.10
and 0.25, 11 targets at ≥ 0.25. Boltz is either consistent on a complex or it is a
coin flip, with nothing in between.

### Recommended floor: min across 5 samples ≥ 0.7 → n = 9

> **Superseded.** This condition was adopted, then *extended* with
> `native epitope recovery >= 0.5` once the report revealed that 9g3e folds
> confidently at the wrong interface. The rule in force is in "The inclusion rule,
> as adopted" below, and it gives 8 targets, not 9. The reasoning here still
> stands as the justification for the iptm half.

Five independently reasonable criteria pick out exactly the same 9 targets:

| Criterion | n | Targets |
|---|---|---|
| `min ≥ 0.7` ("every sample reproduces it") | 9 | 8qf4 9g3e 9e7l 9rye 9min 9ty2 9mnz 9q6n 9b71 |
| `median ≥ 0.7` | 9 | identical |
| `min ≥ 0.6` | 9 | identical |
| `model_0 ≥ 0.7` and `spread ≤ 0.10` | 9 | identical |
| `median ≥ 0.7` and `spread ≤ 0.10` | 9 | identical |
| `model_0 ≥ 0.7` alone | 12 | the 9 plus 9ljj, 9h77, 8rwf |

Because the criteria agree, the stratification is a property of the data rather than of
our threshold choice, and the floor cannot reasonably be accused of being tuned. Say so
in the report: it is the strongest available defence of the inclusion rule.

Use `min ≥ 0.7`. It is the most interpretable form ("the predictor reproduces this
interface on every one of 5 samples"), needs no second clause, and is insensitive to the
threshold: 0.6 and 0.7 give the same answer, so the number is not load-bearing.

**`model_0 ≥ 0.7` alone is the one criterion to avoid.** It admits three targets whose
point estimate hides a coin flip:

| Target | model_0 | min | spread |
|---|---|---|---|
| pdb_00009ljj | 0.830 | 0.364 | 0.467 |
| pdb_00009h77 | 0.726 | 0.333 | 0.394 |
| pdb_00008rwf | 0.722 | 0.236 | 0.486 |

8rwf passes on model_0 at 0.722 while its worst sample is 0.236. A design evaluated on
that target is being compared against a native the predictor only sometimes gets right.

### What the small included set costs, and what to do about it

> Written when the rule gave 9 targets. With the recovery condition added it gives 8,
> and the numbers below are the 8-target ones. The argument is unchanged in kind.

This is a real loss and must not be glossed. B5 was written for "21 paired target-level
outcomes"; the included set gives 8 pairs, and 8 targets × 4 candidates = 32 of the 84
candidates. Bootstrap over 8 clusters has little power to separate hotspot from unguided.

So report both, and say which is primary:

- Primary: the 8 included targets. Thresholds, ranking, and the hotspot-vs-unguided
  comparison live here.
- Secondary: all 21, with each target's native quality reported alongside every
  number, so a reader can see which rows rest on an unreliable control.

Do not drop the 13 excluded targets from the report. "Boltz reproduces 8 of 21 known
nanobody-antigen interfaces well enough to score designs against, and is unreliable or
misplaced on the rest" is a finding about the independent scorer, and it is arguably more
defensible than anything the 84 designs will support.

This also sharpens B3. A's RF2 row with `pred_lddt` 0.87 and
`target_aligned_antibody_rmsd` 36.3 Å was already the headline failure mode; now there
is a second, independent axis of confidence-vs-reality disagreement to compare it to.

## Next actions, in order

1. ~~Finish the harvest and complete the iptm table.~~ Done Nov 24: 21/21 MSAs in
   `runs/boltz_in/msa/`, 21 native controls folded, 0 failures.
2. ~~Re-run all 21 natives at 5 samples.~~ Done Nov 24: `runs/boltz_out_native5/`,
   15 folded + 6 reused, 0 failures, ~16 min. Table above.
3. ~~Freeze the native floor.~~ Done Nov 24, then amended: the rule in force is
   `min of 5 iptm >= 0.7` AND `native epitope recovery >= 0.5`, giving n = 8.
   See "The inclusion rule, as adopted". Original note: `min of 5 ≥ 0.7`, n = 9. Five
   criteria agree on the same 9 targets, so this is close to forced by the data — but it
   is still B's call to sign off, and it must be signed off before any design number is
   *read*. Design folds may be *computed* first: folding leaks nothing into a floor
   defined purely on natives.
4. ~~Run the 88 design folds at 5 samples.~~ Done Nov 24: 88 folded, 0 failed, 21
   natives skipped by the resume check, 109/109 jobs at n=5 in 94 min. The tree is
   `runs/boltz_out5/`, seeded by copying the 5-sample natives in so one
   `fold_metrics.py report` covers natives, RFantibody candidates and MPNN redesigns
   under the same scorer. Recipe, for a re-run:
   ```bash
   cd ~/gen_nanobody/runs/boltz_in
   PATH=$HOME/bz/bin:$PATH OUT=../boltz_out5 SAMPLES=5 \
     nohup bash ../../src/evaluate/run_boltz.sh > ../fold.log 2>&1 &
   ```
   `SAMPLES=5` matters: `run_boltz.sh` defaults to 1 and 1 sample is unusable. Its resume
   check skips any job with a `*model_0.pdb`, so a stale 1-sample output is silently kept
   so use a fresh `OUT` or clear it.
5. ~~Pull the results back.~~ Done Nov 24, verified locally before the instance was
   stopped: 109 model_0 + 109 model_4 + 109 each of confidence/pae/plddt, 21 1-sample
   natives, 21 cached MSAs; every confidence JSON parses with iptm in range and no
   structure has <100 ATOM records. `runs/boltz_out5` 375M, `runs/boltz_out` 61M,
   `runs/boltz_in/msa` 15M, plus the run logs.
6. Instance stopped Nov 24 08:21Z via `sudo shutdown -h now` (no IAM role on the box
   and no local AWS session, so the API path was unavailable). Unreachable on 6 probes.
   **Still open, and never verified all session:**
   - Confirm the console reads Stopped, not Terminated. A shutdown from inside only
     stops if `InstanceInitiatedShutdownBehavior` is `stop`, which is the default but
     cannot be read from inside the instance. Nothing is lost either way: results are
     verified locally and `msa/` removes the need to re-query, but know which it is.
   - Billing → Credits: expiry and service restrictions on the $250.
   - Whether the $50 Budget alert actually exists.
   - `terminate` and delete the 45 GB volume once B3 is done. A stopped instance
     still bills EBS at ~$3.60/mo, and termination does not always take the root volume.

   Instance for reference: `i-0f7ecdbdc88435d6e`, g5.xlarge, us-east-2a, on-demand.
   Total compute spend across the whole session was roughly $4–5.
7. ~~`fold_metrics.py validate` then `report`.~~ Done Nov 24: validate passes on
   21/21, report merged to `results/end_to_end_metrics.csv`. See the Phase 7 findings
   section. Note `cmd_report` opens `--out` with mode "w", so it overwrites, so loop
   with a per-target `--out` and merge, do not reuse one path. Original note: `validate`
   reproduces C's
   5 Å epitope set from the HLT file before any design number is believed; run it first,
   every time. Then commit `fold_metrics.py`, `make_fold_inputs.py` and `run_boltz.sh`
   with the first real results; all three are still untracked.
8. ~~B3.~~ Done Nov 24: `results/interface_metrics.csv`, 252 rows. See the B3 section.
9. ~~B5.~~ Done Nov 24: `results/b5_paired.csv` and `results/b5_failures.csv`. Clustering
   rule, duplicate count, cluster bootstrap and failure taxonomy all discharged. See the
   B5 section.
10. **B4 is the only B item left, and it is blocked on C.** Threshold calibration needs a
   dev set to freeze on, and none exists (no split file in `data/processed/`). Until C
   delivers it, report descriptive distributions only and say so. Ranking must not use
   `pred_lddt` alone, which the native controls now demonstrate rather than assume.

# Nov 24, later: B2, B3 and B5 are complete. Results and findings.

Floor signed off by B: `min iptm over 5 samples >= 0.7`, frozen before any design
number was read. Encoded in `src/evaluate/native_stability.py` (self-check included) and
emitted to `results/native_stability.csv`. 9 of 21 targets stable.

`fold_metrics.py validate` passes on all 21 targets: our contact code reproduces C's
independently derived 5 Å epitope exactly, every target. The metric plumbing is sound.

`fold_metrics.py report` run per target against `runs/boltz_out5`, merged with the
stability flag into `results/end_to_end_metrics.csv`: 109 rows, 24 columns, 21 native +
84 rfab + 4 MPNN redesigns. Per-target files kept in `results/per_target/`.

## The inclusion rule, as adopted

**`min iptm over 5 samples >= 0.7` AND `native epitope recovery >= 0.5`.**
8 of 21 targets included: 8qf4, 9e7l, 9rye, 9min, 9ty2, 9mnz, 9q6n, 9b71.
Encoded with a 4-corner self-check in `src/evaluate/native_stability.py`,
emitted to `results/native_stability.csv`.

Ordering, for the record. The iptm condition was frozen before any design number
was read, as non-negotiable #3 requires. The recovery condition was added *after* the
first report run, on native-only evidence. It is therefore **not blind**, and the report
must say so rather than present 8 targets as the original plan. What makes it defensible
is that it rests entirely on whether a native reproduces its own annotated epitope, which
is a property of the target and the predictor and cannot encode anything about a design.

### Why both conditions are needed

iptm and placement are decoupled in both directions:

| Target | min iptm | native epitope recovery | excluded by |
|---|---|---|---|
| pdb_00009g3e | 0.904 | 0.056 | recovery only — confident at the wrong interface |
| pdb_00009g2g | 0.481 | 0.917 | iptm only — right place, unreliably |
| pdb_00008rwf | 0.236 | 0.733 | iptm only |
| pdb_00009rmk | 0.150 | 0.607 | iptm only |

9g3e contacts 1 of 18 annotated epitope residues at iptm 0.904 with spread 0.045.
Eight natives recover <= 0.333 of their own epitope, five of them exactly 0.000.

Six targets are excluded by exactly one condition (five by iptm only, one by recovery
only), so neither condition is redundant. This is the `pred_lddt` lesson one level up:
confidence of any kind is not placement. Only epitope recovery measures placement.

## Finding 1: the rfab candidates do not engage the intended epitope

Epitope recovery, Boltz from sequence alone, 8 included targets:

| kind | n | median | mean | max |
|---|---|---|---|---|
| native | 8 | 0.887 | 0.820 | 0.955 |
| rfab | 32 | 0.067 | 0.154 | 0.640 |
| MPNN redesign (8qf4 only) | 4 | 0.960 | 0.950 | 0.960 |

Best rfab candidate per target against its native:

| Target | native | best rfab | arm |
|---|---|---|---|
| pdb_00009min | 0.955 | 0.091 | hotspot |
| pdb_00009rye | 0.952 | 0.000 | unguided |
| pdb_00009b71 | 0.941 | 0.412 | hotspot |
| pdb_00009ty2 | 0.895 | 0.000 | unguided |
| pdb_00008qf4 | 0.880 | 0.640 | hotspot |
| pdb_00009q6n | 0.862 | 0.172 | unguided |
| pdb_00009mnz | 0.545 | 0.545 | hotspot |
| pdb_00009e7l | 0.533 | 0.067 | unguided |

**Best rfab candidate reaches native epitope recovery on 1 of 8 targets** (9mnz, a tie),
and on two targets no candidate contacts a single epitope residue. This is the headline
B2 result and it is negative: scored by a predictor that never saw the designed backbone,
RFantibody candidates largely fail to reproduce the intended interface.

The amendment strengthened this result by removing a false positive. Under the
iptm-only rule the count was 2 of 9, and the extra "win" was 9g3e, where the candidate
scored 0.167 against a native of 0.056 — a design beating a control that was itself
broken. Excluding targets whose native misses its own epitope removes exactly that
artifact.

The MPNN redesigns scoring 0.960, above their native, is expected and is a useful
positive control for the metric: they keep the native backbone and pose and change only
CDR sequence, so the interface survives by construction. They are not a de novo arm and
must never be pooled with the rfab candidates.

## Finding 2: hotspot conditioning helps placement, and only placement

Median epitope recovery per arm, paired within each included target:

| Target | hotspot | unguided | |
|---|---|---|---|
| pdb_00008qf4 | 0.580 | 0.540 | H |
| pdb_00009b71 | 0.324 | 0.088 | H |
| pdb_00009mnz | 0.364 | 0.228 | H |
| pdb_00009e7l | 0.067 | 0.034 | H |
| pdb_00009min | 0.045 | 0.022 | H |
| pdb_00009q6n | 0.086 | 0.086 | = |
| pdb_00009rye | 0.000 | 0.000 | = |
| pdb_00009ty2 | 0.000 | 0.000 | = |

**hotspot better on 5, unguided better on 0, tied on 3.** Hotspot is never worse on any
target, and where it wins the margin is large — 3.7x on 9b71, ~2x on 9e7l and 9min.

Meanwhile the pooled iptm view shows nothing at all (hotspot 0.856 vs unguided 0.858).
The two measurements agree with the mechanism: hotspots steer *where* the nanobody binds,
not how confidently it folds, so an effect belongs in epitope recovery and not in iptm.
An effect appearing in exactly the metric it should, and absent from the one it should
not, is mild evidence it is real.

The amendment weakened this result, and that is reported rather than hidden. Under
the iptm-only rule it was 6-0-3 with a sign-test p ~ 0.03; dropping 9g3e, which was a
hotspot win, leaves 5-0-3 at two-sided p = 0.0625. The rule was applied consistently
to both findings, which is why one got stronger and one got weaker.

A sign test is not the prescribed analysis. B5 still owes the bootstrap over
nanobody/antigen clusters, and with 8 pairs the power is thin. Report as promising and
underpowered, never as established.

## B3 done: interface and geometry on all three structure sets

`src/evaluate/interface_metrics.py` (self-check included), output
`results/interface_metrics.csv`: 252 rows, 84 candidates x 3 structure sets, no
duplicates. The three sets are A's ProteinMPNN output in its generated pose
("generated"), A's RF2 prediction ("rf2"), and B's Boltz prediction from sequence
alone ("boltz").

### Two measurement decisions that the data forced

**The generated complex has no side chains.** A's ProteinMPNN output carries exactly
N, CA, C, O. RF2 and Boltz carry full side chains. So every contact count is computed
twice, on backbone atoms only (comparable across all three sets) and on all heavy atoms
(comparable between rf2 and boltz, and matching C's epitope definition). The gap is not
cosmetic: RF2's median full-epitope recovery reads 0.333 on backbone and 0.614 on heavy
atoms. `has_sidechains` marks which rows support the heavy-atom reading. The clash count
for "generated" is likewise not comparable to the other two: fewer atoms, fewer chances
to overlap.

**A's files continue antigen numbering after the nanobody.** Chain T runs 117..183 where
Boltz numbers chain B 1..67, so A's residues are shifted back by an offset read from the
structure and then checked against the manifest's antigen length and contiguity. Observed
offsets are exactly 116 for every A file and 0 for every Boltz file. A silently wrong
offset would move every lookup onto the wrong residue and still produce plausible numbers,
so the self-check covers it directly.

Rosetta InterfaceAnalyzer is deliberately absent, which PLAN.md permits. Contacts, clashes
and buried area cover the same ground without a licence or a build.

### Parser bug found and fixed; B2 results unaffected

`fold_metrics.read_heavy` was keeping hydrogens. RF2's PDBs have an empty element column,
so the parser fell back to the first character of the atom name, and PDB puts a numeric
prefix in column 13 for 4-character hydrogen names: "1HB" read as element "1". Per RF2
file that admitted 813 hydrogens of 2813 atoms as heavy atoms, inflating contacts and
clashes (one test candidate's clash count fell from 438 to 199 after the fix).

Checked on all three input types: Boltz outputs and C's HLT files both populate the
element column and contain zero hydrogens, and A's ProteinMPNN output has no hydrogens
either. Only RF2 triggers it, and RF2 was not scored in B2, so **no B2 number changes**.
Fixed by stripping leading digits before reading the element, with a regression assertion
in `fold_metrics.py --demo` built from a synthetic no-element-column "1HB" line.

### Finding 3: every RF2 interface is physically invalid, and Boltz's are not

Inter-chain heavy-atom pairs closer than 2.0 A:

| set | n | median | mean | max | with >0 clashes |
|---|---|---|---|---|---|
| rf2 | 84 | 96 | 131 | 711 | **84 of 84** |
| boltz | 84 | 0 | 0.5 | 4 | 34 of 84 |
| generated (backbone only, not comparable) | 84 | 0 | 0.0 | 1 | 4 of 84 |

Not a parsing artifact. Heavy-atom counts are identical between the rf2 and boltz files
for the same candidate (1445, 2060, 1698, 3228, 2006 on the five spot-checked), there are
no duplicate coordinates and no altlocs, and RF2's minimum inter-chain distance runs
0.30 to 0.57 A against Boltz's 1.86 to 2.29 A. The chains interpenetrate. No RF2 candidate
is clash-free; 50 of 84 Boltz predictions are, and Boltz's worst case is 4 marginal
contacts just under 2.0 A rather than overlap.

**This flatters the designs, which is the point.** RF2 reports median full-epitope
recovery of 0.614 on heavy atoms against Boltz's 0.086, because interpenetrating chains
"contact" everything near them. Its buried area tail is inflated the same way (max
2173 A^2 against Boltz's 1404), so BSA is not interpretable for rf2 at all. PLAN.md
already excluded RF2 as the final scorer because it sits inside the pipeline. This is an
independent, physical reason: its geometry is not valid, so any interface metric computed
on it reads high for the wrong reason.

### Finding 4: hotspot conditioning works at generation and does not survive

Hotspot recovery, backbone atoms, denominator 4 hotspots per target:

| set | arm | n | median | mean | hitting >=1 |
|---|---|---|---|---|---|
| generated | hotspot | 42 | 0.250 | 0.268 | **74%** |
| generated | unguided | 42 | 0.000 | 0.071 | **21%** |
| rf2 | hotspot | 42 | 0.250 | 0.339 | 76% |
| rf2 | unguided | 42 | 0.500 | 0.417 | 79% |
| boltz | hotspot | 42 | 0.000 | 0.083 | 24% |
| boltz | unguided | 42 | 0.000 | 0.065 | 19% |

In the generated pose the conditioning plainly works: 74% of hotspot-arm candidates touch
at least one hotspot against 21% of unguided-arm candidates. Through RF2 the separation is
gone, and the unguided arm even reads higher, which given Finding 3 is what an invalid
geometry would do. Under independent prediction from sequence alone both arms sit near
zero, 24% against 19%.

The hotspot signal is therefore a property of the generated backbone and its docked pose,
not of the designed sequence. That is consistent with Finding 1: the sequences do not
encode the intended interface strongly enough for an independent predictor to recover it.
It also explains why B5's arm effect on full-epitope recovery is real but small.

## B5 done: paired comparison, cluster bootstrap, failure taxonomy

`src/evaluate/b5_analysis.py` (self-check included). Outputs
`results/b5_paired.csv` (21 rows, one per target, the unit of analysis) and
`results/b5_failures.csv` (84 candidates with geometry and failure classes).

Clustering rule, which B5 owed. `data/processed/targets.csv` carries SAbDab's
`ab_ag_cluster`, and all 21 targets hold distinct values because C kept one
representative per cluster by resolution. Each target is therefore its own
antibody/antigen cluster, and a bootstrap over clusters is a bootstrap over targets.
The script asserts this instead of assuming it, so a future benchmark that puts two
targets in one cluster fails loudly rather than resampling them as independent.

Sequence duplicates, recomputed here: 84 candidates, 83 distinct. The one duplicate
is on 8qf4, so that target's arm comparison is partly a self-comparison.

### Finding 5: the arm effect on epitope recovery is real but marginal at n=8

Paired per target, hotspot median minus unguided median:

| scope | pairs | hotspot better | unguided better | tied | sign test |
|---|---|---|---|---|---|
| 8 included targets | 8 | 5 | 0 | 3 | p = 0.0625 |
| all 21 targets | 21 | 13 | 3 | 5 | p = 0.0213 |

On the included targets hotspot is never worse on any target, and the mean paired
difference is +0.059. Across all 21 the effect holds at +0.064 with a bootstrap 95% CI
of [+0.007, +0.124] and P(diff <= 0) = 0.013, and that bootstrap is meaningful because
three targets do favour unguided.

**The bootstrap on the 8 included targets is degenerate and must not be quoted.** With
5 positive differences, 3 ties and no negatives, a resampled mean can only land at or
below zero if every draw happens to be a tie, so the reported P(diff <= 0) = 0.000391
is exactly (3/8)^8, the chance of drawing only ties. It says nothing about the effect
size. `b5_analysis.py` now prints that warning itself whenever a scope has no negative
differences, so the number cannot be lifted out of the output by mistake. Report the
sign test, p = 0.0625, for the included set.

So: suggestive, consistent in direction, and short of conventional significance on the
targets we trust. The all-21 result is stronger but leans on targets whose native
control is unreliable. Both belong in the report, with the inclusion status attached.

No arm effect on hotspot contact at all (included 1-0-7, p = 1.0; all 21 5-2-14,
p = 0.45), which is exactly what Finding 4 predicts: the hotspot signal lives in the
generated pose, not in the sequence Boltz is given.

### Finding 6: the designs fold correctly and bind in the wrong place

Boltz prediction against A's generated complex, CA atoms. Thresholds are round numbers
fixed in advance and used only to label failures, never to select candidates.

| measurement | median (8 included) | median (all 84) | max |
|---|---|---|---|
| framework RMSD | 0.42 A | 0.40 A | 0.6 A |
| CDR RMSD, after framework superposition | 1.25 A | 1.35 A | 5.5 A |
| pose RMSD, antigen-superposed | 43.7 A | 43.7 A | 68.7 A |

| failure class | 8 included (n=32) | all 84 |
|---|---|---|
| wrong pose (>10 A) | 94% | 96% |
| clashes (>0 pairs <2 A) | 41% | 40% |
| confident but displaced (iptm>=0.7, recovery<0.1) | 47% | 33% |
| CDR refolding (>3 A) | 3% | 6% |
| framework distortion (>2 A) | 0% | 0% |
| no failure class | 1 candidate | 2 candidates |

This is the sharpest result in the project. The fixed-framework invariant holds almost
perfectly in prediction, at 0.4 A median and never worse than 0.6 A, and the designed
CDRs are largely reproduced at 1.35 A. Not one candidate shows framework distortion.
What fails is placement: the nanobody sits a median 44 A from where it was designed to
bind, on 96% of candidates.

The sequences therefore encode the designed fold and do not encode the designed
interface. That single sentence explains Finding 1 (7% epitope recovery), Finding 4
(hotspot conditioning dying between the pose and the sequence) and Finding 5 (a real
but small arm effect). The failure taxonomy PLAN.md anticipated turns out to be almost
entirely one class.

Only 2 of 84 candidates carry no failure class, and both are from the hotspot arm:
8qf4 hotspot nb1 (pose 9.2 A, iptm 0.91, epitope recovery 0.64, on an included target)
and 9ljj hotspot nb1 (pose 8.0 A, recovery 0.24, on an excluded one).

Cross-validated against A's independent measurement. A recorded
`target_aligned_antibody_rmsd` on its RF2 predictions; this is the same quantity on a
different predictor. Medians agree in scale, 43.7 A here against 34.9 A for A, and 77
of 84 candidates exceed 10 A in both while none is under 10 A in both. The per-candidate
correlation is weak (r = 0.27), which is expected rather than troubling: two predictors
that both place the nanobody somewhere wrong need not be wrong in the same direction.
The agreement that matters is the conclusion, and it is independent of the predictor.

## Smoke test on the 109 folds, NOT the RQ3 result

Pooled medians of `model_0` iptm from `runs/boltz_out5/`. This treats the 84 candidates
as independent observations, which B5 forbids. It is here only to show the batch
produced sane numbers. The real comparison is 9 paired target-level outcomes with
bootstrap over nanobody/antigen clusters, and it must wait for the floor sign-off.

| Scope | native | rfab | hotspot | unguided |
|---|---|---|---|---|
| All 21 targets | 0.726 (n=21) | 0.809 (n=84) | 0.806 (n=42) | 0.809 (n=42) |
| 9 iptm-stable targets (pre-amendment) | 0.920 (n=9) | 0.858 (n=36) | 0.856 (n=18) | 0.858 (n=18) |

Three things to carry into B5 as hypotheses, not findings:

1. On the stable targets designs land below their natives (0.858 vs 0.920), which is
   the sensible direction.
2. Hotspot and unguided are indistinguishable (0.856 vs 0.858 stable; 0.806 vs 0.809
   overall). If that survives pairing and bootstrapping it is a negative result for RQ3's
   premise, which is a legitimate finding and should be reported as one.
3. Pooled across all 21, `rfab` (0.809) appears to beat native (0.726). That is an
   artifact of the unstable targets dragging the native median down, and it is the
   clearest demonstration of why the stratification is not optional.

## Cost and time, actual vs estimated

| | Nov 16 estimate | Measured |
|---|---|---|
| Per fold | 2–5 min | ~22 s GPU, ~75 s wall (each job is its own `boltz predict`, so each pays model load) |
| 109 folds | 6–12 GPU-hours | ~2.5 h |
| Folding cost | $10–15 | ~$2.50 at $1.006/hr on-demand |
| Boltz cache | ~3 GB assumed | 7.6 GB (2.2 ckpt + 2.0 aff + 1.8 mols + 1.8 mols.tar) |

Credit was never the constraint and is even less of one now. `rm ~/.boltz/mols.tar`
reclaims 1.8 GB if the 43 GB volume gets tight.

## Traps that cost time, so they are not re-paid

- `python3-dev`. See above. Reads as a CUDA failure, is not one.
- The MSA guard is global. `run_boltz.sh` scans all `*/*.yaml` for missing MSAs
  regardless of `ONLY=`, so you cannot run one target until all 21 CSVs exist. Correct
  for the batch, awkward for one-offs — call `boltz predict` directly for those.
- `pgrep -f "boltz predict"` self-matches. A polling loop whose own ssh command line
  contains that string sees itself and reports RUNNING forever. Use a pidfile
  (`echo $! > /tmp/boltz.pid`, then `kill -0 $(cat ...)`).
- macOS rsync is old. `--info=stats2` is unsupported; use `--stats`.
- Don't search the AMI name in the console. Searching "Ubuntu Server 22.04" returns
  only a paid SQL Server bundle, and name searches land in the Community AMIs tab where
  every publisher, Amazon's included, shows as a bare account ID and looks untrusted.
  Use the Quick Start tile's dropdown, or resolve by SSM public parameter.
- Boltz caps at Python < 3.13, so 24.04 LTS (3.12) is right and 26.04 would fail.

## Still true, and still binding

The four non-negotiables are unchanged. Boltz has not touched generation or selection
and never sees a designed backbone, so #1 holds. #3 is why the native-iptm floor must be
written down before the design folds are read, not after.
