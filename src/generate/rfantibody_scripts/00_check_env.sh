#!/bin/bash
# Fails loudly before any GPU time is spent. Run from the repo root.
set -uo pipefail

REPO=${REPO:-$PWD}
RFAB=${RFAB:-}
RFAB_RUN=${RFAB_RUN:-python}
FRAMEWORK=${FRAMEWORK:-}

fail=0
say() { printf '%-10s %s\n' "$1" "$2"; }

if [ -z "$RFAB" ]; then
  say FAIL "RFAB is not set. export RFAB=/path/to/RFantibody"
  fail=1
else
  say ok "RFAB=$RFAB"
fi

INF=$RFAB/scripts/rfdiffusion_inference.py
if [ -n "$RFAB" ] && [ -f "$INF" ]; then
  say ok "inference script $INF"
else
  say FAIL "not found: $INF"
  fail=1
fi

if [ -z "$FRAMEWORK" ]; then
  for c in \
    "$RFAB/scripts/examples/example_inputs/h-NbBCII10.pdb" \
    "$RFAB/include/rfantibody/examples/example_inputs/h-NbBCII10.pdb"
  do
    [ -f "$c" ] && FRAMEWORK=$c && break
  done
fi
if [ -n "$FRAMEWORK" ] && [ -f "$FRAMEWORK" ]; then
  nh=$(awk '/^ATOM/ && substr($0,22,1)=="H"' "$FRAMEWORK" | wc -l)
  nr=$(grep -c 'PDBinfo-LABEL' "$FRAMEWORK" || true)
  say ok "framework $FRAMEWORK (H atoms $nh, CDR remarks $nr)"
  if [ "$nr" -eq 0 ]; then
    say FAIL "framework has no REMARK PDBinfo-LABEL lines; RFantibody needs HLT"
    fail=1
  fi
else
  say FAIL "no framework PDB found. export FRAMEWORK=/path/to/h-NbBCII10.pdb"
  fail=1
fi

for f in data/processed/epitope_residues_v2.csv \
         data/processed/antigen_residue_mapping_v2.csv \
         data/processed/prepped_v2.csv \
         data/processed/structural_qc_final.csv
do
  if [ -f "$REPO/$f" ]; then say ok "$f"; else say FAIL "missing $f (git pull?)"; fail=1; fi
done

n=$(ls "$REPO"/data/processed/prepped_v2/*.pdb 2>/dev/null | wc -l)
if [ "$n" -eq 21 ]; then say ok "prepped_v2 has 21 PDBs"; else say FAIL "prepped_v2 has $n PDBs, expected 21"; fail=1; fi

if command -v python3 >/dev/null; then say ok "python3 $(python3 -V 2>&1)"; else say FAIL "no python3"; fail=1; fi

echo
if [ "$fail" -eq 0 ]; then
  echo "ALL CHECKS PASSED"
else
  echo "CHECKS FAILED -- fix the FAIL lines above before Step 2"
  exit 1
fi
