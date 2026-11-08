#!/bin/bash
# Local replacement for 03_rfdiff_batch.sbatch when there is no Slurm.
#
# Usage: 03_rfdiff_local.sh [max_T_residues] [num_designs] [arms]
#   03_rfdiff_local.sh            # all 21 targets, 20 designs, both arms
#   03_rfdiff_local.sh 120 10     # only targets with T <= 120, 10 designs each
#   03_rfdiff_local.sh 999 20 hotspot
#
# Runs smallest antigen first and KEEPS GOING past a failure, because on a
# consumer GPU the large targets are expected to run out of memory and stopping
# there would throw away the small ones that did fit. Which sizes fit is itself
# a reportable number, so every attempt lands in runs/rfab/manifest.csv whether
# it succeeded or not.
set -uo pipefail

MAXT=${1:-100000}
NDES=${2:-20}
ARMS=${3:-"hotspot unguided"}

REPO=${REPO:-$PWD}
CSV=$REPO/data/processed/rfab_inputs/rfab_targets.csv

if [ ! -f "$CSV" ]; then
  echo "missing $CSV -- run 01_make_rfab_targets.py first" >&2
  exit 1
fi

if command -v nvidia-smi >/dev/null 2>&1; then
  echo "GPU:"
  nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | sed 's/^/  /'
else
  echo "GPU: no nvidia-smi found; RFdiffusion will fall back to CPU and be very slow"
fi
echo

mapfile -t PLAN < <(python3 - "$CSV" "$MAXT" <<'PY'
import csv, sys
rows = [r for r in csv.DictReader(open(sys.argv[1]))
        if int(r["T_residues"]) <= int(sys.argv[2])]
rows.sort(key=lambda r: int(r["T_residues"]))
for r in rows:
    print("%s %s" % (r["PDB_ID"], r["T_residues"]))
PY
)

echo "${#PLAN[@]} targets with T <= $MAXT, $NDES designs each, arms: $ARMS"
echo

ok=0; bad=0
for entry in "${PLAN[@]}"; do
  read -r ID TLEN <<< "$entry"
  for arm in $ARMS; do
    echo "### $ID (T=$TLEN) $arm"
    if bash "$REPO/src/generate/02_rfdiff_one.sh" "$ID" "$arm" "$NDES"; then
      ok=$((ok + 1))
    else
      bad=$((bad + 1))
      echo "!!! $ID $arm failed (exit $?) -- continuing with the next target"
    fi
    echo
  done
done

echo "=== done: $ok runs ok, $bad failed"
echo "next: python3 src/generate/04_report_status.py"
