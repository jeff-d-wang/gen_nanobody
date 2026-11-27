# Task Brief — Final Week Execution Plan

**23 Nov 2026 · final report due 29 Nov · repository status checked 23 Nov**

---

## 1. The problem, restated

### Overall goal

Given an antigen and a specified surface patch, **can our pipeline generate a nanobody candidate that computationally engages that patch, and can its scoring components distinguish more plausible from less plausible designs?**

We have **no experimental binding measurements for the newly generated candidates**. Our report must separate (i) completed generation, (ii) structural/computational evaluation, and (iii) evaluation against **independent experimental labels**. These are not interchangeable forms of evidence.

### Three research questions

| | Question | Evidence required |
|---|---|---|
| **RQ1** | Can scores distinguish **experimentally screened binders** from nonbinders? | Experimental binding labels, corresponding pair-level scores, cluster-aware PR analysis |
| **RQ2** | Do scores track **measured nanobody affinity**? | Traceable experimental $K_D$, comparable target/assay groups, scores for those same pairs |
| **RQ3** | Does **explicit hotspot conditioning** improve computational engagement of the intended patch? | Paired hotspot/unguided candidate evaluation on the 21 targets |

The candidate-scoring hypothesis has three channels: **structural confidence, interface geometry, and interface energy**. A combined score is a hypothesis to test, not a result we can claim today.

### Critical qualification of RQ3

Both RFantibody arms use **the same antigen crop, derived from the native nanobody interface**. The unguided arm omits explicit hotspot labels; it is **not** a blind search across the entire antigen surface. RQ3 therefore measures the incremental effect of *explicit hotspot information within a reference-derived crop*, not de novo epitope discovery.

---

## 2. Where we stand: verified evidence versus missing evidence

### Repository snapshot

- **`main`**: the 21-target structural benchmark, native H/T PDBs, epitope and residue mappings, QC, native-backbone ProteinMPNN CDR recovery, and midterm figures are committed.
- **`RFantibody`**: two commits ahead of `main` at the time of inspection; generation runner, manifest/configuration, audit, handoff documentation, and summary CSVs are committed. They are **not yet merged into `main`**; no pull request was found at inspection time.
- **Raw run archive**: `runs/rfab_benchmark_1ce3344e4822/` is described in Aidan's handoff, but **its complete structures are not in GitHub** (the `runs/` directory is gitignored). The team must transfer and audit the archive before structure-level evaluation can be considered executable.

### Data that are available in GitHub

| Layer | File(s) | Meaning |
|---|---|---|
| Native reference | `data/processed/targets.csv`, `prepped_v2/*.pdb` | 21 curated native nanobody–antigen complexes; H = native nanobody, T = cropped antigen |
| Native epitope | `epitope_residues_v2.csv`, `antigen_residue_mapping_v2.csv` | 446 native-contact antigen residues (5 Å rule) and 3,055 original-to-T mapping rows |
| Quality | `structural_qc_final.csv` | 20 PASS, 1 CHECK (`8yj8` native H-chain warning) |
| Separate pipeline sanity check | `results/cdr_recovery.csv` | 1,344 rows = 21 targets × 16 sequence samples × 4 summaries (`cdr1`, `cdr2`, `cdr3`, `all`); **not** RFantibody-generated candidates |
| Generated-run summaries | `results/generation_status.csv`, `results/rfab_candidate_metrics.csv` on `RFantibody` | 42 completed target/condition runs; 84 candidate rows of RF2 confidence and design-to-reprediction RMSDs |
| Raw outputs, not yet shared through GitHub | Aidan's `runs/rfab_benchmark_1ce3344e4822/` archive | 84 reported RF2 best PDBs, corresponding generated backbones and sequence-designed structures, state/attempt metadata, FASTA, logs |

Aidan's audit **reports** 84/84 candidates passing file and sequence-lineage checks and 4.34 h of accumulated *model-process elapsed time*. A recipient must rerun the audit against the transferred archive. **4.34 h is not measured GPU-hours.**

### What we cannot yet claim

1. **RQ3 is not yet answered.** There is no committed 84-candidate epitope-contact comparison or target-paired hotspot/unguided analysis.
2. **Confidence has not been shown to fail or succeed at distinguishing a hit from a miss.** The hit/miss label and its relationship to RF2 confidence have not been computed.
3. **Candidate interface energy is not yet reported.** The native benchmark's contact analysis is not the generated candidates' geometry table.
4. **RQ1 and RQ2 are not yet evaluated.** Their experimental datasets and structure-linked scoring tables are absent.
5. A well-folded, confidently predicted candidate is **not** a demonstrated binder.

---

## 3. The dependency map

**Track G — Generated-candidate evaluation (RQ3):**

`Aidan's raw archive` → `Jeffery's 84-candidate contact/geometry table` → `21-target paired hotspot-versus-unguided analysis` → `RQ3 figure and limitations`

**Track E1 — Experimental binding labels (RQ1):**

`AVIDa-hIL6 labels + sequence-cluster split` → `feasible subset with pair-level structures/scores` → `confidence / geometry / energy / combined comparison` → `RQ1 metrics`

**Track E2 — Experimental affinity (RQ2):**

`curated experimental K_D records` → `matched target/assay groups and structural coverage audit` → `pair-level structures/scores` → `within-group affinity-ranking analysis`

**Do not join the 84 generated designs to AVIDa or affinity labels.** They are different biological pairs. Scoring *code* can be shared; scores and experimental labels can be joined only for the **same antibody–antigen pair**.

**Priority order:** secure the generated structures and finish RQ3 first. Treat full structure-based RQ1/RQ2 as **conditional** on data coverage and compute feasibility, not guaranteed final-week deliverables.

---

## 4. The work

### B — Jeffery: repository integration, interface evaluation, and RQ3

#### B0. Accept the archive and consolidate the repository — **23 Nov**

**Inputs:** `RFantibody` branch, `WORK_ARCHIVE_RFantibody.md`, Aidan's full `runs/rfab_benchmark_1ce3344e4822/` directory.

1. Review/merge `RFantibody` into `main` through a PR. A merge moves **code and tracked CSVs**, **not** the gitignored structures.
2. Coordinate with Aidan to transfer the raw archive through agreed team storage (or an explicitly configured artifact mechanism). Do **not** silently commit large files to ordinary Git.
3. Run `python3 src/generate/audit_rfantibody.py runs/rfab_benchmark_1ce3344e4822` on the received copy. Require `passed: true`, 42 conditions, 84 validated candidates, and zero errors.
4. Use each condition's `state.json` to locate **accepted** stage directories. Do not hardcode `attempt_1`; a rerun/interrupted attempt can coexist with the accepted one.

**Done when:** the team has a verified, reproducible way to access all 84 accepted final PDBs, their source backbones and designed structures, plus a manifest keyed by `(pdb_id, arm, backbone_index)`.

#### B1. Compute generated-interface contacts — **24–25 Nov; highest scientific priority**

**Inputs:** accepted RFdiffusion and RF2 PDBs, native epitope and residue mapping CSVs, hotspot list from Aidan's run manifest.

For each of the 84 candidates, define antigen contact residues **on chain T** using a prespecified minimum heavy-atom distance of **5 Å** to chain H. Use the native epitope (original numbering) only after mapping it to chain T. Let `P` be RF2-predicted contacts, `G` generated-backbone contacts, `E` native-epitope residues, and `S` the four specified hotspots.

- `epitope_coverage = |P ∩ E| / |E|` (the **primary RQ3 endpoint**)
- `epitope_precision = |P ∩ E| / |P|` (undefined if `P` is empty)
- `hotspot_coverage = |P ∩ S| / |S|` (evaluate for **both** arms)
- `contact_retention = |G ∩ P| / |G|` (undefined if `G` is empty)
- Number of interface contacts / contacting residues, severe clashes, and existing RF2 `pred_lddt` / `interaction_pae`

Keep generated and RF2-predicted contacts separate. Do not substitute target-aligned RMSD or global confidence for epitope engagement. Specify atom selection, missing-atom handling, and the clash cutoff in the output metadata.

**Output:** `results/evaluation/rfab_candidate_interfaces.csv` (**84 rows, one per accepted candidate**) and reproducible script(s) in `src/evaluate/`.

**Done when:** every candidate has an epitope-contact result or an explicit failure reason, and sample IDs match the accepted run state.

#### B2. Compare hotspot and unguided arms — **25–26 Nov**

Aggregate the **two candidates within each target × arm** (predefine mean and a secondary best-of-two view), then take **one paired difference per target**. Report the 21 paired differences, their distribution and uncertainty by **resampling targets**, not treating 84 candidates as independent antigen targets. Include representative successes **and** failures.

Evaluate whether `pred_lddt` and `interaction_pae` are associated with measured epitope coverage; only now can we examine whether confidence is informative about hit/miss. Avoid dichotomizing contacts unless thresholds were fixed independently.

**Outputs:** `results/evaluation/rfab_paired_targets.csv`, one comparison figure, one contact/confidence diagnostic figure, and a short method/limitations paragraph.

**Done when:** RQ3 has a transparent paired analysis. A null or mixed result is a valid result.

#### B3. Add interface energy/geometry — **26 Nov if feasible**

Use a consistent preparation protocol on the **same 84 RF2 final structures**. Compute Rosetta InterfaceAnalyzer terms such as `dG_separated` (Rosetta energy units, **not** measured thermodynamic ΔG), `dSASA_int`, optional shape complementarity, and explicit clash measures. Record failures and version/settings. Check sensitivity to antigen size and interface area; consider normalization and within-target comparisons.

**Output:** `results/evaluation/rfab_rosetta.csv`, keyed by candidate ID, with a method note. Do not block B1/B2 or RQ3 reporting if Rosetta cannot be installed or standardized by the cutoff.

#### B4. Secondary diagnostics — only after B1–B3

Classify plausible failure modes (wrong epitope, docking movement, CDR refolding, clashes, framework distortion), using separate target-aligned versus framework-aligned RMSDs. Investigate `9my8`'s native-backbone CDR3 recovery outlier as a **separate ProteinMPNN sanity-check result**; do not confuse those 336 redesigns with Aidan's 84 generated candidates.

---

### A — Aidan: raw-data handoff, compute accounting, and RQ1

#### A0. Deliver and validate the generated-run archive — **23 Nov**

Give Jeffery the **entire** `runs/rfab_benchmark_1ce3344e4822/` directory, including `state.json`, accepted PDBs, FASTA, `inputs/`, config, audit, and logs. Supply a checksum or manifest, plus the exact path/access instructions in a shared team location. Confirm that another teammate can execute the audit **without rerunning generation**.

**Done when:** B0 passes on the recipient's copy. This is a prerequisite for RQ3.

#### A1. Cost accounting — **24–25 Nov**

Use `results/generation_status.csv` and the run manifest to report counts of: completed conditions, generated backbones, unique sequences (after deduplication), audited RF2 predictions, and **candidates passing a stated computational filter** if such a filter is later fixed. Plot model-process **wall time** versus candidate count. Call it GPU-hours only if actual utilization/allocation data substantiate GPU-hours.

**Output:** `results/evaluation/rfab_compute_summary.csv` plus a clearly labeled figure/table.

#### A2. AVIDa-hIL6 label benchmark (RQ1) — **24–26 Nov, gated**

1. Download and version the experimental binding/nonbinding labels and paired sequences; retain target variant and assay metadata where available. Experimental nonbinder labels apply to that assay, not universal nonbinding.
2. Prevent leakage by splitting on **nanobody sequence clusters**, keeping measurements of related VHHs and target variants grouped appropriately. Record class prevalence; a balanced subsample must not be reported as population-prevalence AUPRC without adjustment.
3. Establish a **small throughput pilot first**, including whether an antibody–IL-6 complex can actually be modeled/scored using available compute. Start with a manageable paired positive/negative set; the original proposal's 1,000–2,000 folded pairs is **not** an established feasible sample size.
4. If paired structures and all relevant score channels are available, report AUPRC (with its prevalence baseline), AUROC, and precision-at-k; use cluster-aware intervals. Compare confidence, geometry, energy, and a prespecified combined score **on the same experimental pairs**. Freeze weights/rules on development data before final evaluation.
5. If structure-based scoring is infeasible by **26 Nov**, deliver the curated split, class balance, a sequence-only baseline where feasible, and an explicit limitation. **Do not invent structure-based RQ1 results.**

**Outputs:** `data/processed/avida_*` (small shareable tables or a download script), `results/evaluation/avida_feasibility.csv`, and, only if supported, `results/evaluation/rq1_scores.csv` / RQ1 figure.

---

### C — Charlie: affinity curation, experimental-pair matching, and RQ2

#### C0. Curate measured affinities — **23–25 Nov**

**Sources to examine:** ANDD, SAbDab-nano, sdAb-DB and their original cited experiments. These resources **overlap**; do not count database copies as independent observations.

**Minimum record fields:** `record_id`, `antibody_sequence`, `antigen_identity`, `antigen_construct`, `KD_value`, `KD_unit`, `KD_M`, `pKD`, `measurement_type`, `assay_method`, `assay_conditions`, `publication_or_DOI`, `source_database`, and `structure_id` if available. Preserve raw and normalized values separately.

- Keep **experimentally measured** affinities; exclude computational predictions.
- Do not treat censored observations (`<`, `>`, detection limit) as exact values. Store them with censoring flags or exclude from ordinary Spearman.
- Deduplicate by experimental record/publication and pair identity, while preserving legitimate replicate measurements.
- Define comparable groups by **antigen construct plus assay method/conditions**; do not rank unrelated antigens using raw $K_D$.
- For each group, count **distinct nanobody sequences**, usable exact affinities, and availability of experimental or realistically modelable pair-level complex structures.

**Outputs:** `data/processed/affinity_curated.csv`, `results/evaluation/affinity_coverage.csv`, and a provenance/selection note. These are useful deliverables even without a full RQ2 correlation analysis.

#### C1. RQ2 feasibility decision — **25 Nov**

Produce a table of candidate target/assay groups and score availability. Prefer groups with **at least five distinct experimentally measured nanobodies**, consistent assay context, and plausible structure access. This is a pragmatic screening threshold, not a universal statistical rule.

**Go** only if the exact experimental pairs can be assigned **their own** comparable structural scores in time. **No-go** if data are mostly single observations per antigen or lack workable complex structures; in that case report affinity coverage and the limitation, not a cross-dataset correlation.

**Important:** Jeffery's `rfab_rosetta.csv` scores the **84 generated sequences**, not the experimental nanobodies in the affinity table. It cannot directly provide RQ2 labels or features. Reuse Jeffery's **scoring method** only after generating/finding structures for the experimental pairs.

#### C2. Affinity-ranking test — **26 Nov, conditional**

Where pair-matched scores genuinely exist, calculate Spearman correlations between individual/combined scores and **$pK_D=-\log_{10}(K_D/\mathrm{M})$** **within each comparable target/assay group**; report group sample sizes and uncertainty. Keep this separate from the generated-design experiment. Weights for any combined score require independent development/calibration data; do not optimize them on the same group used to claim performance.

**Output if feasible:** `results/evaluation/rq2_pair_scores.csv` and a within-group rank-correlation figure/table. Otherwise output the C0/C1 data audit with a transparent no-go conclusion.

#### C3. Evaluation hygiene and reporting support — **26–28 Nov**

Check sequence/antigen overlap between new experimental tables and the 21 structural benchmark targets; distinguish exact sequence overlap, homologous-family similarity, and repeated measurements. Document split and source provenance. Integrate the validated RQ2 results or limitations into the report; help keep the three research questions' sample units separate.

---

### Shared scoring contract — everyone

Two **distinct** sample tables, one common metric definition:

- **Generated table (RQ3)** key: `(pdb_id, arm, backbone_index)`; native epitope is a **reference annotation**, not an experimental binding label.
- **Experimental-pair table (RQ1/RQ2)** key: `(dataset, experimental_pair_id, antigen_construct, assay_group)`; labels/$K_D$ apply to this exact sequence–antigen pairing.

Reusable metric columns can include `pred_lddt`, `interaction_pae`, `epitope_coverage` (generated benchmark only), `interface_contact_count`, `dSASA_int`, `dG_separated`, `clashes`, `scorer_version`, `structure_source`, and missingness flags. **A matching column name is not permission to join different biological pairs.**

Compare confidence-only, geometry-only, energy-only, and combined scoring **only within a dataset where all candidate-level scores and appropriate outcomes exist**. Missing channels stay missing; do not silently score them as zero. Do not use RF2's in-pipeline prediction as an independent final predictor.

---

## 5. Minimum viable outcomes and stop rules

| Track | Minimum acceptable final-report output | What would overclaim |
|---|---|---|
| **RQ3 (required)** | Audited 84-candidate lineage; contact/epitope metrics; 21 paired hotspot/unguided target summaries; caveats about reference-derived crop | “Hotspot guidance is ineffective/effective” without contact analysis, or “generated candidates bind” |
| **Interface scoring (high priority)** | Reproducible geometry metrics; Rosetta if feasible, with units and limitations | Equating Rosetta energy with measured $K_D$ or experimental binding free energy |
| **RQ1 (conditional)** | AVIDa split and prevalence; structure-modeling feasibility; scored classification if pair-level structures exist | Transferring RQ3 scores onto AVIDa labels; using a balanced AUPRC as real-world prevalence performance |
| **RQ2 (conditional)** | Deduplicated traceable measured-$K_D$ table; comparable assay groups and structure-coverage report; within-group ranking only if scores exist | Correlating 84 generated scores with unrelated published affinity values |

**Go/no-go checkpoint: 26 Nov.** After this date do not start another large structure-modeling campaign. Finish only analyses supported by existing, verified inputs and write the limitations clearly.

If the raw archive cannot be transferred, the team can still analyze the committed RF2 confidence/RMSD CSV descriptively, **but cannot claim contact-level RQ3 evaluation or Rosetta scoring of all 84 structures**.

---

## 6. Schedule and handoffs

| Date | Action / owner | Concrete end-of-day check |
|---|---|---|
| **23 Nov** | **A** delivers raw archive; **B** opens/reviews integration PR and verifies audit; **C** starts affinity source inventory | Archive accessible to B; archive audit status recorded; affinity sources identified |
| **24 Nov** | **B** implements contact evaluator and tests a few candidates; **A** creates AVIDa label split and small modeling pilot; **C** normalizes measured-$K_D$ records | Candidate-ID mapping works; AVIDa feasibility measured; first provenance table exists |
| **25 Nov** | **B** completes 84-candidate contact table and 21 target-pair summaries; **A** records compute costs/class prevalence; **C** completes affinity coverage and structure-link audit | RQ3 primary metrics reproducible; RQ1/RQ2 candidate coverage known |
| **26 Nov** | **B** completes RQ3 statistics and Rosetta if feasible; **A/C** decide go/no-go on experimental-pair structural scoring | Freeze analyses and scope; no unsupported RQ1/RQ2 promise remains |
| **27 Nov** | Everyone prepares figures, tables, error analysis, method detail | All plotted results trace to saved tables/scripts |
| **28 Nov** | Everyone writes/reviews final report; **C** integrates data/limitations | Full draft, provenance and limitations audited |
| **29 Nov** | Final checks and submission | Reproducibility, claims, figures, references and artifact links verified |

---

## 7. Non-negotiable interpretation and reporting rules

1. **No binding claim:** no generated candidate has been experimentally verified.
2. **No premature RQ3 verdict:** hotspot vs unguided outcome requires a defined contact endpoint and paired target analysis.
3. **No leakage:** thresholds/score weights selected on development data must not be tuned on the 21 held-out structural targets or the final experimental test set. Recent PDB deposition does not itself rule out pretrained-model exposure.
4. **Independent predictor:** RF2 was part of the generation/reprediction pipeline; if an excluded final predictor cannot be run, explicitly say so.
5. **Statistical units:** candidate designs from the same target are dependent. For RQ3 summarize at the 21-target level; for AVIDa cluster related nanobodies; for affinity respect target/assay groups.
6. **Compute units:** 4.34 h is summed recorded stage elapsed time, **not** active GPU-hours or energy cost.
7. **Full provenance:** record sample identities, successful/failed model attempts, preprocessing, residue numbering, source papers, assay conditions and scoring versions.
8. **Honest partial completion:** RQ1/RQ2 can legitimately end with a documented feasibility limit. Do not fill gaps with simulated or unmatched experimental labels.

---

## 8. Immediate messages to the team

**Aidan:** Please transfer the complete accepted RFantibody run archive to our shared location today and confirm Jeffery can audit it. In parallel, prepare AVIDa label/cluster splits and benchmark the actual cost of scoring a small number of experimental pairs before proposing a 1,000–2,000-structure run.

**Jeffery:** Please review/merge the `RFantibody` branch, verify the received archive, and prioritize the 84-candidate epitope-contact table and 21-target paired comparison. This is our shortest path to a defensible RQ3 result. Rosetta scoring follows if time permits.

**Charlie:** Please curate experimentally measured, traceable and deduplicated nanobody $K_D$ entries. The first decision is not “how many records exist?” but “how many *comparable experimental pairs with scorable structures* exist?” Report that coverage before attempting Spearman correlations.

**Team-wide decision:** By **26 Nov**, lock what has actually been measured, cut unsupported experiments, and devote the remaining time to error analysis and a reproducible report.

---

## 9. FOR AI: Repository-Grounded Execution Protocol

**Audience:** Any AI assistant implementing, reviewing, or extending a task in this brief. This is an execution contract, not a replacement for the scientific questions and acceptance criteria above. **Read the actual repository and available run archive before writing code or claiming that an analysis is complete.**

### 9.1. Establish current state, not the dated snapshot

1. Read the latest relevant branches and open PRs of [`jeff-d-wang/gen_nanobody`](https://github.com/jeff-d-wang/gen_nanobody), including `main` and `RFantibody`; inspect `PLAN.md`, current source code, configs, report sections, and the two `WORK_ARCHIVE_*.md` handoffs as applicable. Record the inspected branch, commit SHA, and date.
2. Inspect **existing analysis scripts and results before implementing duplicate work**. In particular, check whether 84-candidate contact metrics, paired RQ3 statistics, confidence-discrimination analyses, and their scripts have already been added or shared outside GitHub. If present, **audit and extend** them instead of restarting the analysis. This document's task status is a planning snapshot, not proof of present-day incompletion.
3. Determine whether the full RFantibody run archive is accessible. `runs/` is gitignored in the currently documented repository; GitHub summary CSVs are **not** a substitute for the accepted input/output PDBs, `state.json`, run manifests, and audit files. If the archive is unavailable, state which analyses are blocked and do not invent or regenerate their results.

### 9.2. Discover and validate the data interface

Before performing an analysis, write a short **input/output contract** identifying:

- **Inputs:** repository-relative source paths; any external/archive paths; exact columns and data types; row counts; units; required and optional fields; and the source commit or archive checksum.
- **Identifiers:** the authoritative candidate key (e.g., `pdb_id`, `arm`, backbone index), experimental pair key, and rules for joining tables. **Do not join generated candidates to unrelated experimental binding or affinity labels.**
- **Structure conventions:** H/T chain identities, residue numbering and insertion codes, antigen cropping, epitope/hotspot definitions, and explicit mappings between native, prepared, generated, and predicted structures. **Do not assume a constant RF2 T-chain numbering offset (such as 116) without validating it on the actual accepted structures.**
- **Accepted outputs:** resolve the active stage and accepted attempt through each target/arm's `state.json` and run manifest; do not count all `*_best.pdb` or `*_attempt_*` files indiscriminately. Preserve candidate lineage across RFdiffusion, ProteinMPNN, and RF2.
- **Outputs:** intended file paths under `src/evaluate/` and `results/`, schema, metric definitions, failure flags, and numerical acceptance checks.

Inspect actual CSV/PDB contents and the code that produces them. Test representative cases, including atypical residue numbering, insertion codes, missing atoms, modified residues, and interrupted/retried model stages. A convention observed on one structure is **not** automatically valid across 84 candidates. If the brief and repository disagree, use the repository to establish *existing behavior* and this brief to establish *scientific intent*; document and resolve the conflict before interpreting results.

### 9.3. Implement with explicit scientific definitions

- For RQ3, distinguish **native-epitope coverage** (fraction of reference epitope residues contacted by the RF2-predicted nanobody), **hotspot coverage**, **epitope precision**, and **design-to-reprediction contact retention**. These are different endpoints. State the heavy-atom cutoff, residue and chain mapping, the numerator and denominator of each metric, and how missing/zero-denominator cases are handled. Do not silently drop unmappable reference residues or change denominators.
- Interpret `target_aligned_antibody_rmsd` and `target_aligned_cdr_rmsd` as **design-to-reprediction pose consistency after target alignment**, not RMSD against the native crystallographic nanobody. Keep framework-aligned loop changes separate from docking changes.
- When testing hotspot versus unguided, aggregate the two candidate designs within each target/condition and compare the **21 paired targets**. Report per-target values, paired effect sizes, confidence intervals, prespecified tests, and any multiple-testing adjustment. Do not treat all 84 designs as independent targets. Both arms see a reference-derived antigen crop; “unguided” does **not** mean epitope-naive full-antigen search.
- For confidence-versus-contact analyses, avoid overinterpreting nonsignificant Spearman correlations as proof of no signal. Account for clustering within targets (e.g., target-level bootstrap/permutation); treat post-hoc median-based high/low-contact groups as exploratory, not calibrated hit/miss truth.
- For Rosetta scoring, document input preparation, scoring protocol and software version; `dG_separated` is a model energy estimate, **not** experimentally measured binding free energy or $K_D$. Do not use an outcome-derived feature (e.g., native-epitope overlap) to predict that same outcome and call it independent discrimination.
- For RQ1/RQ2, score **the experimental antibody–antigen pairs themselves** before comparing with binding labels or measured affinities. Deduplicate overlapping database sources, preserve assay context and units, respect original label prevalence or explain weighting, and keep cluster/target/assay groups intact in splits and uncertainty estimates.

### 9.4. Test, audit, and hand off

1. Run a small representative validation set before running the full analysis. Include at least one normal structure and one known edge case; verify chain/residue mapping and accepted-output selection against the archived artifacts.
2. Run reproducible checks on the complete outputs: expected target/arm/candidate counts when applicable (21/42/84), unique IDs, joins with no unexplained unmatched records, finite coordinates/metric values, valid score ranges, and independently recomputed summary statistics. **Expected counts are checks, not a reason to force results to match.**
3. Save scripts in `src/evaluate/`, analysis tables in `results/`, and a short methodological note with the exact command, input version, output schema, units, statistical test/adjustment, exceptions, and limitations. Avoid absolute machine-specific paths in team-facing deliverables.
4. Before declaring a task done, report **(a) verified completed work, (b) files/commits produced, (c) tests passed, (d) unresolved assumptions or blocked inputs, and (e) precisely which scientific claim the results support**. Do not infer task completion from the presence of a filename or from a prior assistant's summary.

**Stop rule:** If a required artifact, join key, structural mapping, experimental label, independent scoring condition, or validation check is missing or contradictory, stop the affected inference. Preserve any valid partial output and explain the blocker; **do not fabricate a substitute measurement, silently relax the criterion, or label exploratory evidence as a confirmed answer.**
