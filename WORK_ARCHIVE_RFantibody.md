# RFantibody generation: results and evaluation handoff

Prepared November 8, 2026. All paths below are relative to the repository root
unless explicitly stated otherwise.

## 1. Summary and current status

The small RFantibody benchmark run is complete. All 21 prepared antigen targets
were processed in two conditions: hotspot-guided and unguided. Each condition
generated two backbones, designed one sequence per backbone with ProteinMPNN,
and predicted each resulting complex with RF2.

| Deliverable | Completed |
|---|---:|
| Antigen targets | 21 |
| Target–condition runs | 42 / 42 |
| Accepted generated backbones | 84 |
| Designed sequence records | 84 |
| Best RF2 prediction structures | 84 |
| Candidates passing the output and lineage audit | 84 / 84 |
| Outstanding audit errors | 0 |
| Recorded model-process elapsed time, all attempts | 15,634.59 seconds / 4.34 hours |

**What this establishes:** the workflow runs end to end across the full target
set on the local NVIDIA RTX 3060 Ti (8 GB), and the accepted outputs are
structurally readable and consistent across stages.

**What it does not establish:** 84 distinct sequences, 84 good binders, successful
epitope targeting, or superiority of hotspot guidance. Sequence uniqueness,
interface quality, and condition comparisons still need analysis. No binding
experiments have been performed on these designs.

The next stage is evaluation of these existing candidates, not another
installation exercise or an automatic rerun of generation.

## 2. Exact run and experimental settings

Run directory: `runs/rfab_benchmark_1ce3344e4822/`.
In this document, `RUN/` means that directory; it is shorthand, not a shell variable.

| Setting | Value |
|---|---|
| Antigen source | Chain T from `data/processed/prepped_v2/*_H_T.pdb` |
| Native nanobody | Chain H excluded from generation inputs |
| Framework | Upstream `h-NbBCII10.pdb`; snapshot in `RUN/framework.pdb` |
| Conditions | `hotspot` and `unguided` |
| Backbones per target/condition | 2; deterministic design indices/seeds 0 and 1 |
| Designed loop lengths | H1: 7, H2: 6, H3: 8 |
| Diffusion steps | 50 |
| ProteinMPNN | One sequence per backbone; H1/H2/H3 only; temperature 0.1; deterministic mode |
| RF2 | One recycle; seed 37; hotspot-show proportion 0.0 |
| Upstream revision | `3e5a5bb1cb98befde9fc049ae3ff58cabb491f67` |

The RF2 setting of one recycle includes the initial prediction plus recycling;
it should not be described as the upstream default ten-recycle protocol.
The full configuration, model-weight hashes, input hashes, and framework hash
are in `RUN/config.json`; `configs/rfab_benchmark.json` is the report-facing copy.

Both conditions receive the same cropped antigen and independent framework.
Only the hotspot condition receives explicit hotspot residue labels.
**Unguided does not mean epitope-naive:** both conditions receive a crop derived
from the known native interface.

Hotspots are selected by closest native contact with a preference for 8 Å
separation between C-alpha atoms. Spacing is relaxed if necessary to select
four residues. This is a deterministic heuristic, not a burial calculation.
The selected hotspots, full epitope in T numbering, and spacing-relaxation flag
are recorded in `RUN/inputs/targets.json`.

## 3. Files the evaluation team should start with

| File | Purpose |
|---|---|
| `results/generation_status.csv` | One row per target/condition; completion, output counts, elapsed time, errors |
| `results/rfab_candidate_metrics.csv` | One row per accepted best RF2 prediction; confidence and structural-deviation scores |
| `RUN/designed_sequences.fasta` | Designed heavy-chain sequences for all 84 accepted candidates |
| `RUN/candidate_metrics.csv` and `.json` | Run-local versions of prediction metrics |
| `RUN/audit.json` and `audit.csv` | Validation results: 42 conditions, 84 candidates, no errors |
| `RUN/inputs/targets.json` | Hotspots, full epitope, input hashes, QC and input-conversion notes |
| `RUN/inputs/targets/*_T.pdb` | Exact antigen-only inputs archived for this run |
| `RUN/config.json` and `framework.pdb` | Run settings and actual starting framework |
| `RUN/<pdb_id>_<arm>/state.json` | Authoritative record of accepted stage outputs and all attempts |
| `RUN/<pdb_id>_<arm>/*_attempt_*/output.log` | Model output and diagnostics for each attempt |
| `src/generate/benchmark_rfantibody.py` | Active benchmark runner |
| `src/generate/audit_rfantibody.py` | Output and sequence-lineage audit |

The older generation shell scripts are reference material; use the Python
benchmark runner for this run. The standalone lysozyme installation pilot in
`runs/rfantibody_pilot_smoke/` is separate and must not be pooled into this dataset.
Likewise, `results/cdr_recovery.csv` describes the earlier native-backbone
sequence-recovery exercise, not these newly generated candidates.

## 4. Match each candidate across stages

Use `(pdb_id, arm, backbone index)` as the candidate identity. For example:

`pdb_00009mnz`, `hotspot`, `0` identifies one candidate.

Within `RUN/pdb_00009mnz_hotspot/`, the corresponding files are:

```text
rfdiffusion_attempt_1/outputs/nb_0.pdb
proteinmpnn_attempt_1/outputs/nb_0_dldesign_0.pdb
rf2_attempt_1/outputs/nb_0_dldesign_0_best.pdb
```

The FASTA header is `pdb_00009mnz_hotspot_nb_0_dldesign_0`.
The final `0` is the ProteinMPNN sequence sample; there is one sample per backbone.
The RF2 `_best` suffix identifies its selected best prediction, not an extra design.

**Do not assume every accepted attempt is attempt 1.** Read the `stages` entries
in each `state.json` for the accepted directories. An interrupted attempt was
retained during the overnight handoff. Recursively counting every PDB in the
run directory would count abandoned or duplicate files as additional candidates.

Metrics tables and state files contain absolute WSL paths beginning with
`/mnt/c/Users/adnja/Documents/GitHub/gen_nanobody/`. On another machine, replace
that repository prefix with the recipient's checkout location. Keep the rest of
the run-directory layout unchanged. The candidate identity should not depend on
the original absolute path.

## 5. What the existing metrics mean

These scores compare RF2 predictions with the generated/designed input complex,
not with the original native nanobody. The meanings below follow the installed
RF2 implementation in `rf2/modules/model_runner.py`.

| Column | Meaning and interpretation |
|---|---|
| `pred_lddt` | Aggregate RF2 predicted lDDT confidence, stored on a 0–1 scale here; higher means greater predicted local confidence. It is not CDR-only or interface-only confidence. |
| `interaction_pae` | Mean predicted alignment error over interchain residue pairs; lower indicates greater predicted relative-placement confidence. Units: Å. |
| `pae` | Mean of the saved model PAE tensor when serialized to the score table; not the full matrix. Units: Å. |
| `target_aligned_antibody_rmsd` | Antibody C-alpha displacement after aligning the target between design and prediction; reflects retention of the designed docking arrangement. Units: Å. |
| `target_aligned_cdr_rmsd` | CDR C-alpha displacement after target alignment; sensitive to both loop changes and movement of the nanobody relative to the antigen. Units: Å. |
| `framework_aligned_antibody_rmsd` | Antibody C-alpha displacement after framework alignment; focuses on monomer structural consistency. Units: Å. |
| `framework_aligned_cdr_rmsd` | CDR C-alpha displacement after framework alignment; measures loop-shape consistency independently of docking movement. Units: Å. |
| `framework_aligned_H1/H2/H3_rmsd` | The same framework-aligned comparison for each heavy-chain CDR. Units: Å. |
| `framework_aligned_L1/L2/L3_rmsd` | Not applicable: these nanobodies have no light chain. Blank CSV cells / JSON nulls are expected, not zeros. |

A confident nanobody fold can still move away from its designed interface.
Read confidence, framework-aligned RMSD, and target-aligned RMSD together;
do not rank candidates solely by `pred_lddt`.

The score tables inherit the upstream PDB score rounding to two decimal places.
Do not imply greater precision. The runner does **not** export full PAE or
per-residue confidence arrays, ipTM, experimental affinities, Rosetta interface
scores, contact recovery, or clash metrics. Those are not hidden columns in this
table; additional extraction or computation is required where appropriate.

## 6. Recommended next evaluation work

### A. Preserve the handoff and establish the evaluation table

1. Copy the entire `RUN/` directory, not just the small CSV summaries.
2. Include this document, the generation/audit scripts, configuration, and the
   benchmark metadata and mapping CSVs from `data/processed/`.
3. Verify the copied output with the audit command below.
4. Build one candidate table keyed by target, condition, and backbone index,
   joining the designed sequence, generated backbone path, predicted structure
   path, hotspots, and existing RF2 metrics.

The audit verifies file counts, finite coordinates and applicable scores,
H/T chains, CDR labels, input hashes, antigen sequence preservation, fixed
framework sequence during ProteinMPNN, and sequence preservation during RF2.
It does not validate binding, clash freedom, or correct epitope engagement.

### B. Evaluate structure and interface consistency

For each candidate, compute contact to the intended hotspots and full epitope
in both the generated and predicted complexes. Report the definition explicitly
(for example, minimum heavy-atom distance and a prespecified cutoff), including
the denominator used for any recovery fraction. Keep generated contacts,
predicted contacts, and retained-contact fractions separate.

Use target alignment for docking comparisons and framework alignment for loop
shape comparisons. Inspect representative overlays, including disagreement
between high local confidence and poor target-aligned consistency. Add steric
clash and interface-geometry checks; Rosetta scoring can provide a separate
assessment if available.

### C. Compare conditions without treating 84 candidates as independent targets

Summarize the two designs within each target/condition, then compare hotspot
versus unguided results as paired target-level outcomes. Report distributions
and paired differences for the 21 targets. Do not select only attractive examples
or claim statistical strength from 84 independent observations.

Keep related sequences/targets together in uncertainty estimates and follow the
project's cluster-aware evaluation policy. Two designs per condition are a
small feasibility sample; conclusions about diversity, sampling depth, or
general superiority need more evidence. Count sequence duplicates and define
any clustering criterion explicitly.

### D. Preserve data separation and independent final evaluation

These 21 targets come from the held-out structural benchmark. Do not tune
success thresholds against their observed scores and then call the result a
held-out pass rate. Fix scoring and thresholds using separate development data;
until then, report descriptive measurements rather than calibrated success claims.

RF2 is inside the design pipeline. If RF2 scores are used for candidate selection,
they are not independent final-scoring evidence. The team still needs to choose
and run a separate final predictor under the project protocol, or explicitly
document inability to meet that requirement. Its inputs should not supply the
generated CDR backbone or docked pose as a template. No such independent final
prediction is included in this handoff.

The binding-classification and affinity-ranking research questions also require
their separate experimental-label datasets; this generation run does not answer them.

## 7. Input and accounting caveats

- Antigen numbering is the cropped chain-T numbering. For original residue
  identifiers use `antigen_residue_mapping_v2.csv`, not the raw integer alone.
- Generated heavy-chain positions are not native IMGT residue numbers. Use each
  generated structure's CDR remarks and residue ordering for loop assignments.
- The MSE at 8rwf T85 was converted to MET for model parsing, with SE renamed
  SD and its element changed to sulfur. Coordinates and numbering were retained.
  The prepared antigen has 165 residues; an ATOM-only count of the canonical
  structure can omit the MSE. This conversion is recorded in the input manifest.
- The native 8yj8 complex has an H-chain geometry warning. That native H chain
  was removed; the QC flag remains documented rather than silently discarded.
- The 4.34-hour total sums recorded model-stage elapsed time, including the
  interrupted attempt. It is not measured active GPU time, energy consumption,
  peak-memory usage, or the elapsed interval between the first and last calendar
  timestamps. Preflight work and time between invocations are not included.
- Backbones on disk are not the same as candidates passing a scientific filter.
  Any yield-versus-compute figure must state which yield it counts.

## 8. Sharing and verification

`runs/` and generated `data/processed/rfab_inputs/` are gitignored. A Git push of
scripts and summary CSVs alone will **not** deliver the structures to a partner.
Share the run directory separately through the team's agreed storage and verify
the received copy. This document does not imply a backup or transfer has occurred.
Several new scripts and results are currently untracked and still need normal
repository review/commit if they are to be shared through Git.

From the repository root in WSL/Linux, validate the saved run without running
any models:

```bash
python3 src/generate/audit_rfantibody.py runs/rfab_benchmark_1ce3344e4822
```

Expected result: `passed: true`, `target_conditions: 42`,
`validated_candidates: 84`, and an empty `errors` list.
This command regenerates the short `RUN/HANDOFF.md` and
`results/RFANTIBODY_RUN.md`; this longer archive is maintained separately.

No generation command is needed to start analysis. If additional generation is
later authorized, preserve this frozen run and distinguish new settings and
outputs from this feasibility dataset.

## 9. Suggested report wording

We completed an end-to-end RFantibody feasibility run on 21 prepared antigen
targets using an independent nanobody framework. For each target, two backbones
were generated with explicit hotspot guidance and two without it; one sequence
per backbone was designed with ProteinMPNN and predicted with RF2. All 42
target–condition runs completed, yielding 84 candidates that passed file,
structure, and sequence-lineage checks. Recorded model-process elapsed time
across all attempts was 4.34 hours on the local RTX 3060 Ti system. These results
establish execution feasibility; interface quality, hotspot-guidance benefit,
independent structural assessment, and experimental binding remain separate
evaluation questions.
