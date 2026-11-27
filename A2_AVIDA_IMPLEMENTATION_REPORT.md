# A2 AVIDa-hIL6 implementation report

## Frozen input/output contract (November 26, 2026)

Source checkout: `RFantibody` at `44d4b710c6ec9530113ce106d41f4ac0ba4baefa`.
Protocol: `configs/avida_protocol.json`, written before splits/model metrics.
Source: official COGNANO/AVIDa-hIL6 revision
`cd783d573e015f83364e42287810ad597850aeed`, CC BY-NC 4.0.

Inputs are `AVIDa-hIL6.csv` (VHH_sequence, Ag_label, binary label, subject_species,
subject_name, subject_sex) and `antigen_sequences.csv` (Ag_label, Ag_sequence).
The release documents 573,891 pairs, 38,599 VHH sequences, 31 antigen variants,
and 20,980 positives; these are verification targets, not forced output counts.
Raw downloads remain under ignored `data/raw/avida_hil6/<revision>/` with pinned
URLs, byte counts and SHA-256 digests in a tracked source manifest.

Authoritative experimental pair identity is a content hash of dataset, exact VHH
sequence, antigen construct/sequence and released assay group. Original row and
subject metadata are retained. Antigens are joined many-to-one by Ag_label;
unmatched antigens, noncanonical sequences, nonbinary labels and conflicting pair
labels fail validation. Exact duplicate rows may be collapsed only with counts
reported. No generated RFantibody candidate is joined to these labels.

The split groups full VHH sequences using connected components of all pairwise
global edit-distance comparisons at distance/max(lengths) <= 0.10. This is a
precisely defined edit-similarity rule, not a claim of Chothia/CDR clonotyping.
Indels count as edits. All IL-6 variant measurements of a VHH and its component
remain in one split. Whole clusters are assigned deterministically without label
stratification; the intended train/development/test proportions are 70/15/15.
No split seed or threshold will be chosen using test performance. Antigen variants
can occur in multiple splits: this tests new VHH clusters against known variants,
not unseen-antigen generalization or the paper's antigen-held-out benchmark.

The prespecified sequence baseline is an additive logistic classifier using
train-only VHH trigram TF-IDF and train-only antigen variant one-hot encoding.
It has no structural channels and no antibody/variant interaction features.
No class balancing or evaluation subsampling is planned. Average precision
(step-integrated PR), AUROC and precision at 100 will use the original split
prevalence. Uncertainty resamples entire VHH components with all their variant
measurements. Model settings and seeds are fixed in the protocol, not tuned on test.

Planned outputs: reproducible download/curation/evaluation scripts; ignored
`data/processed/avida_*` full tables; tracked source/split/validation summaries,
`results/evaluation/avida_feasibility.csv`, sequence-baseline metrics and figures,
and a four-pair development-only pilot (two positives, two negatives, distinct
components). The pilot is for feasibility, not prevalence-adjusted performance.

Structural scoring is gated on exact experimental-pair structures, an appropriate
predictor and all claimed score channels. Missing channels remain missing. A
failed input/tool preflight is recorded as such, not as measured folding throughput.
No `rq1_scores.csv` with claimed structural comparisons will be created without
those measurements. This report will record actual outputs and validation below.

## Completed curation and split

The downloaded tables passed the release-count checks: **573,891 distinct pairs,
20,980 positives, 552,911 negatives, 38,599 VHH sequences and 31 variants**. No rows
were duplicated, dropped, relabeled, balanced or imputed. All antigen joins matched.
Released subject fields are retained. The assay group denotes the released
phage-display/biopanning/NGS label set; per-pair library, replicate, concentration
and affinity metadata are not present in these source CSVs and were not invented.
Nonbinder labels describe this screening assay, not universal nonbinding.

The exact computation evaluated 744,922,101 unordered sequence pairs and found
4,297,707 edges under the fixed 90% edit-similarity rule. Connected components
produced **10,796 clusters**. All components stay in one split; all measurements
of a VHH against different variants stay together. The largest component contains
15,650 VHHs and 180,555 experimental pairs. Components are sorted by descending
pair count, with seeded content hashes breaking ties, then assigned to the split
with the lowest current count divided by its target fraction. Labels are not used
to choose assignments. This is a new benchmark split, not the authors' split.

| Split | Pairs | Positive | Negative | Positive prevalence | VHHs | Clusters |
|---|---:|---:|---:|---:|---:|---:|
| Train | 401,723 | 20,356 | 381,367 | 5.0672% | 28,601 | 7,530 |
| Development | 86,084 | 318 | 85,766 | 0.3694% | 5,062 | 1,632 |
| Test | 86,084 | 306 | 85,778 | 0.3555% | 4,936 | 1,634 |
| Full release | 573,891 | 20,980 | 552,911 | 3.6557% | 38,599 | 10,796 |

**The prevalence shift is substantial.** Whole related-sequence components,
including the large component, were preserved rather than split to make the
labels look balanced. No seed/threshold search was performed after examining
these results. The test results apply to this one cluster-held-out split and do
not estimate performance at the full release prevalence, clinical prevalence,
unseen IL-6 variants, or a different alpaca's repertoire.

All released VHHs have a six-histidine tail and most share a leader-like prefix.
The full released sequences are preserved in both clustering and the baseline.
Shared flanks can inflate whole-sequence similarity; connected components can
also connect distant endpoints through intermediate sequences. These choices
conservatively group many related sequences but do not establish absence of all
biological homology or pretrained-model exposure. No undocumented trimming or
CDR numbering was introduced.

## Sequence-only baseline results

Train-only preprocessing learned 7,872 trigrams and the 31 known antigen categories.
The fixed SGD logistic model converged in seven iterations without warnings.
Training took 2.618 seconds in the recorded run; the complete baseline script,
including bootstrap statistics and figures, took 77.501 seconds. These are CPU
process elapsed measurements, not GPU-hours.

| Metric | Development | Test | Test 95% cluster-bootstrap interval |
|---|---:|---:|---:|
| Average precision (step-integrated PR area) | 0.068599 | 0.008754 | 0.004503–0.015501 |
| Prevalence / no-skill PR reference | 0.003694 | 0.003555 | 0.001610–0.006348 |
| AUROC | 0.642513 | 0.708768 | 0.647894–0.780808 |
| Precision at 100 experimental pairs | 0.21 | 0.00 | 0.00–0.00 |

The held-out top 100 contain **no positive pairs**. The modest ranking result is
not evidence of a useful high-precision screening system, and the large
development/test AP difference must not be hidden. An additive variant feature
does not model sequence-specific effects of antigen mutations. No structure,
interface-energy or combined-channel classifier was fitted.

Intervals use 1,000 bootstrap draws of whole VHH components, preserving every
member and variant measurement, with the fitted model held fixed. They do not
include retraining, split-selection or model-selection variability. Scores are
ranked descending with pair-ID order breaking ties independently of labels.
Precision-at-k counts experimental pairs; it is not a diversity-constrained
shortlist of unique nanobodies. The zero-width top-100 interval is an empirical
bootstrap outcome, not proof that the population quantity is known exactly.
No significance tests or model-comparison p-values were used.

## Files and reproducible commands

New implementation files:

- `configs/avida_protocol.json`: frozen source revision, clustering, split,
  model, evaluation and pilot-selection rules.
- `src/data/download_avida.py`: pinned official downloads, cache/checksum checks.
- `src/data/prepare_avida.py`: input validation, exact components, splits, class
  balance, source-row provenance and deterministic development-pilot selection.
- `src/evaluate/avida_baseline.py`: fixed sequence baseline, cluster-aware
  intervals, saved predictions/coefficients and PNG/PDF precision-recall figures.
- `src/evaluate/avida_pilot.py`: bounded four-pair RF2 pilot with explicit target
  mapping, measured subprocess time, confidence and geometry diagnostics.
- `src/evaluate/test_avida.py`, `test_avida_pilot.py`: validation/regression tests.
- `src/evaluate/requirements-avida.txt`: direct dependency versions;
  `requirements-avida-lock.txt`: complete environment used (CPython 3.10.22).

Shareable outputs under `results/evaluation/`:

- `avida_source_manifest.json`: source revision, URLs, license, hashes and sizes.
- `avida_curation_audit.json`: counts, component checks, timing and curated hashes.
- `avida_class_balance.csv`: 128 rows, split × variant plus totals; pairs,
  positives, negatives, prevalence, distinct VHHs and component counts.
- `avida_pilot_pairs.csv`: four exact development pairs with released sequences,
  labels, subject metadata, source rows, antigen constructs and cluster IDs.
- `avida_sequence_baseline_metrics.csv`: eight metric rows with point estimates,
  percentile intervals, denominators, bootstrap counts and scoring times.
- `avida_sequence_baseline_pr.png` / `.pdf`: sequence-baseline PR curves.
- `avida_baseline_run.json`: fixed settings reference, dependency versions,
  timing, warnings and model/prediction hashes.
- `avida_pilot_preflight.json`, `avida_feasibility.csv`,
  `avida_feasibility_decision.json`: measured modeling feasibility and gate decision.

Large/regenerable outputs remain ignored under existing repository rules:

- `data/raw/avida_hil6/<revision>/`: original source CSVs and dataset card.
- `data/raw/avida_pilot/1ALU.pdb`: hash-checked public IL-6 template.
- `data/processed/avida_pairs.csv.gz`: all 573,891 curated experimental pairs.
- `data/processed/avida_vhh_splits.csv.gz`: 38,599 sequence/cluster/split records.
- `data/processed/avida_sequence_predictions.csv.gz`: all 172,168 development/test
  pair IDs, labels, clusters, variants and sequence-only scores.
- `runs/avida_baseline/`: data-only model coefficients and feature vocabulary.
- `runs/avida_pilot/`: per-pair JSON inputs, model logs, accepted PDBs and attempt records.
- `.venv/avida/`: isolated analysis environment. `.gitignore` also excludes the
  incomplete `.venv-avida/` directory left by the initial unavailable-ensurepip setup.

No existing RFantibody installation, A1 files, generation archive, report TeX,
or supplied status brief was modified. No commit or remote push was made.

From the repository root in WSL/Linux:

```bash
uv venv --python 3.10 .venv/avida
uv pip install --python .venv/avida/bin/python -r src/evaluate/requirements-avida-lock.txt
.venv/avida/bin/python src/data/download_avida.py
.venv/avida/bin/python src/data/prepare_avida.py
.venv/avida/bin/python src/evaluate/avida_baseline.py
.venv/avida/bin/python -m unittest discover -s src/evaluate -p test_avida.py -v
```

The pilot requires the existing sibling RFantibody installation/weights and CUDA:

```bash
../RFantibody/.venv/bin/python src/evaluate/avida_pilot.py
../RFantibody/.venv/bin/python -m unittest discover -s src/evaluate -p test_avida_pilot.py -v
```

Pilot execution requires an empty output directory and refuses to overwrite a
previous pilot. Use `--out runs/avida_pilot_new` only for an intentional new run.
Tests validate saved outputs without rerunning models. The two environments are
separate so analysis dependencies do not alter the established RFantibody install.

## Sources and interpretation

- [Official AVIDa-hIL6 dataset](https://huggingface.co/datasets/COGNANO/AVIDa-hIL6),
  pinned revision and source checksums recorded locally, CC BY-NC 4.0.
- [Tsuruta et al., NeurIPS 2023](https://proceedings.neurips.cc/paper_files/paper/2023/file/8339dacd9df7ffe9623760f74169dd1e-Paper-Datasets_and_Benchmarks.pdf):
  dataset and experimental label-generation provenance. Cite this work when
  reporting these results. The labels derive from screening rather than measured
  affinities; the released labeled subset excludes unresolved observations.
- [RCSB PDB 1ALU](https://www.rcsb.org/structure/1ALU): antigen-only structural
  template, not an experimentally matched structure of any pilot VHH complex.

This implementation delivers the brief's curated-label, cluster-split,
sequence-baseline and modeling-feasibility outcomes. It does not answer whether
confidence/geometry/energy/combined scores distinguish experimental binders.

## Measured four-pair structural pilot and final gate

All four exact experimental pairs were selected from distinct development
clusters against wild-type IL-6 (two labeled positive, two negative). Selection
was deterministic by pair ID within label, before looking at structural outputs.
Only unlabeled pair sequences and the antigen structure were passed to RF2.

RF2 ran with the existing frozen upstream revision
`3e5a5bb1cb98befde9fc049ae3ff58cabb491f67` (clean working tree), saved weight hash,
seed 37, one recycle (initial prediction plus one recycling step), no hotspots,
and a 600-second timeout per pair, on the local RTX 3060 Ti (8 GB). No antibody
backbone or docking pose was supplied. The public 1ALU antigen-only structure
provided 157 observed backbone residues of the 218-residue released WT construct.
Every observed residue was checked against the released sequence using the PDB
DBREF mapping; 61 unobserved residues, including signal/tag regions, remained
explicitly masked. RF2's Python JSON input uses NaN placeholders for masked
coordinates. These are missing input data, not measured structural coordinates.

| Pilot pair suffix | Assay label | Model-process seconds | RF2 pred_lddt | Interchain PAE (Å) | Target contact residues, 5 Å | Interchain atom pairs below 2 Å |
|---|---:|---:|---:|---:|---:|---:|
| `1db17ffb76b6c06fe0935d79` | 1 | 75.268 | 0.77 | 23.80 | 108 | 1,230 |
| `baf6f2ce5cd6dc14755caedb` | 1 | 58.865 | 0.76 | 21.59 | 77 | 847 |
| `0019b29612c4ee7a0fbf6266` | 0 | 56.363 | 0.76 | 24.88 | 147 | 2,323 |
| `00d407f8924dcf95bc7c92f7` | 0 | 56.054 | 0.77 | 24.08 | 130 | 1,720 |

Total recorded model-process time: **246.550 seconds (4.109 minutes)**;
mean **61.638 seconds per pair**, including process/model startup. Input
preparation, weight hashing and downstream geometry calculation are excluded.
This does not establish throughput for a large batch, longer recycles, other
variants, different constructs or an independent predictor.

All outputs passed finite-coordinate and exact H/T sequence checks. `passed` in
the feasibility CSV means **execution/output validation**, not scientific quality.
Contacts use a minimum H-to-T heavy-atom distance <=5 Å; the severe-close-contact
diagnostic counts interchain heavy-atom pairs at distance <2 Å. These are raw
geometric diagnostics, not Rosetta energies or a calibrated success filter.
Target contact counts include the full released construct; the separate
`template_supported_target_contacts_5A` column restricts contacts to the 157
template-supported positions. H/T residue order is checked against the exact
pair sequence; no constant predicted-PDB numbering offset is assumed.

The very large close-contact counts make these structures unsuitable as evidence
of successful binding or as an unreviewed structural benchmark. Raw construct
tails, missing-template regions and the prediction protocol need further review;
this pilot does not identify which caused the geometry problems. Shared leader
or expression-related sequence features can also confound the sequence baseline.

**November 26 gate: NO-GO for full structure-based RQ1.** Local RF2 can execute
this small pilot, but it is not an independent final predictor; no configured
energy scorer was run, PyRosetta and Boltz were absent from the inspected RF2
environment, and the construct/geometry issues remain unresolved. No large
modeling campaign, balanced-pilot AUPRC, tuned combined score, or structural
`rq1_scores.csv` was produced. The valid completed deliverables are curation,
cluster split/prevalence, a sequence-only baseline, and measured feasibility.

## Final validation and artifact integrity

**12 A2 tests passed:** nine curation/baseline tests plus three pilot tests.
They cover transitive sequence components, the exact indel-distance boundary,
label-conflict and unmatched-antigen rejection, duplicate accounting,
label-independent deterministic split allocation, cluster-preserving bootstrap,
weighted-versus-explicit bootstrap replication, exact prediction/label joins,
and independent AP/AUROC/precision recomputation. Pilot tests verify missing
target masks, rejection of sequence mismatches, all four exact H/T outputs, and
contact/clash counts recomputed with a separate KD-tree distance implementation.
The PR figure was visually inspected. `git diff --check` passed.

`src/evaluate/audit_avida.py` and
`results/evaluation/avida_artifact_manifest.json` additionally record and verify
SHA-256 hashes and byte counts for implementation files, this report, source
downloads, curated data, baseline outputs, and the complete pilot artifacts.

```bash
# Create only after intentionally producing/reviewing new results:
.venv/avida/bin/python src/evaluate/audit_avida.py --write-manifest
# Verify the existing handoff without recomputation:
.venv/avida/bin/python src/evaluate/audit_avida.py
```

There are no unresolved inputs for the delivered fallback scope. Structural-channel
comparison remains explicitly unsupported. Generated-design results (RQ3) and
affinity ranking (RQ2) remain separate and unchanged.
