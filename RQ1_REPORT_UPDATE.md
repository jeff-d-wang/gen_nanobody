# RQ1 report integration

Updated November 27, 2026 against merged `main` at `6d4d1bf`. The report is now
maintained directly in `report/main.tex`, matching the merged single-file layout.
This supersedes the earlier integration into separate section files. No new
experiments were run and saved experiment outputs were not changed.

## Changes

- Wrote RQ1 directly into `report/main.tex` (Section 4.1). It reports
  the pinned AVIDa source, exact edit-distance clustering rule, label-independent
  split assignment, split counts/prevalence, fixed sequence-only model, bootstrap
  protocol, metric table, existing PR figure, four-pair structural pilot, scope
  decision, sources and reproducibility information.
- Updated the abstract, analysis and conclusion to describe
  completed RQ1 work and its limitations. Removed the unsupported abstract claim
  that AlphaFold 3 scoring was performed.
- Clarified that the pipeline diagram describes the intended screening protocol
  and that the structural confidence/geometry/energy comparison remains unfinished.
  RQ2/RQ3 protocols are distinguished from the completed RQ1 subsection.
- Preserved the merged midterm contribution, original schedule and original
  resource plan, explicitly labeling their status statements as historical.
  Added current RQ1 ownership/completion and measured-cost context.
- Scoped missing energy/independent-scoring statements specifically to the AVIDa
  experimental pairs. Newly merged Boltz and generated-candidate work concerns
  different pairs and is preserved rather than treated as missing globally.
- Reused `results/evaluation/avida_sequence_baseline_pr.pdf` directly rather than
  creating a separate copy of the plotted data or changing the experiment figure.
- Removed the unused bibliography backend; the report uses attributed inline
  source links, matching its existing citation approach. Placeholder bibliography
  entries are not printed or cited.
- Inlined the custom definitions previously imported from `pre.tex`, removed
  obsolete section-file boundary comments, and adjusted paragraph/table formatting
  and float placement. No manuscript `input`/`include` commands remain.
- Retained the existing ignored PDF output location. This revision does not
  alter the user's staged changes or recreate a `report/sec/` directory.

## Delivered files

- Main LaTeX entry point: `report/main.tex`.
- RQ1 text: directly in `report/main.tex`, Section 4.1.
- Compiled report: `output/pdf/gen_nanobody_report.pdf`, also refreshed at
  `report/main.pdf` so the report directory does not retain the obsolete PDF.
- All report prose and custom definitions are in `report/main.tex`; the existing
  document class, logos and result figure remain external assets.

## Verification

- Compiled with pdfLaTeX (TeX Live 2024), with a second pass for references.
- Final log contains no LaTeX warnings, unresolved references, overfull boxes or
  underfull boxes.
- Rendered all eleven PDF pages with Poppler and inspected the layout, including
  detailed views of the RQ1 metric table, PR figure and pilot discussion.
- Checked all 16 reported metric values/interval endpoints against the saved
  baseline CSV at the displayed six-decimal precision.
- Extracted the final PDF text and checked its eleven-page count and absence of
  unresolved-reference markers.
- Checked the existing A2 artifact manifest: 29 files match byte-for-byte, while
  19 text files fail the strict byte/hash check because this Windows checkout
  now uses CRLF. All 19 match the recorded hashes after CRLF-to-LF normalization;
  no other content difference was found. The original manifest and the user's
  staged experiment files were left untouched, so the strict audit still reports
  those newline-only mismatches.
- `git diff --check` passed (Git separately notes normal LF/CRLF conversion).

The report preserves the key negative RQ1 results: test precision at 100 is zero;
the pilot structures have many interchain close contacts; no energy channel or
independent final predictor was evaluated for those AVIDa pairs. It does not claim
that structural RQ1 is complete or transfer results from unrelated generated
candidates onto experimental labels. No models were rerun, and no commit, staging
operation or push was made.

## Rebuild

From the repository root in PowerShell:

```powershell
New-Item -ItemType Directory -Force runs/rq1_single_file_build, output/pdf | Out-Null
Set-Location report
pdflatex '-interaction=nonstopmode' '-halt-on-error' '-output-directory=../runs/rq1_single_file_build' main.tex
pdflatex '-interaction=nonstopmode' '-halt-on-error' '-output-directory=../runs/rq1_single_file_build' main.tex
Set-Location ..
Copy-Item -LiteralPath runs/rq1_single_file_build/main.pdf -Destination output/pdf/gen_nanobody_report.pdf
Copy-Item -LiteralPath runs/rq1_single_file_build/main.pdf -Destination report/main.pdf
```

Use a clean build directory if it contains auxiliary files from the old
bibliography-enabled version. Build logs and rendered previews are retained in
the ignored `runs/rq1_single_file_build/` directory.
