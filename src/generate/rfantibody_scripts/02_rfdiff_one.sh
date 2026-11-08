#!/bin/bash
# Usage: 02_rfdiff_one.sh <PDB_ID> <hotspot|unguided> [num_designs]
#
# Looks the target path and hotspot string up in rfab_targets.csv at run time, so
# the hotspot numbering can never drift from the target file it was derived from.
# Idempotent: a target whose output dir already holds backbones is skipped, so a
# requeued array job resumes instead of redoing finished work.
set -uo pipefail

PDB_ID=${1:?usage: 02_rfdiff_one.sh <PDB_ID> <hotspot|unguided> [num_designs]}
ARM=${2:?arm must be hotspot or unguided}
NDES=${3:-20}

REPO=${REPO:-$PWD}
RFAB=${RFAB:?export RFAB=/path/to/RFantibody}
RFAB_RUN=${RFAB_RUN:-python}
FRAMEWORK=${FRAMEWORK:-$RFAB/scripts/examples/example_inputs/h-NbBCII10.pdb}

IN=$REPO/data/processed/rfab_inputs
CSV=$IN/rfab_targets.csv
OUT=$REPO/runs/rfab/${PDB_ID}_${ARM}
MANIFEST=$REPO/runs/rfab/manifest.csv

case "$ARM" in
  hotspot|unguided) ;;
  *) echo "arm must be hotspot or unguided, got '$ARM'" >&2; exit 2 ;;
esac

if ls "$OUT"/*.pdb >/dev/null 2>&1; then
  echo "skip $PDB_ID $ARM: $(ls "$OUT"/*.pdb | wc -l) backbones already in $OUT"
  exit 0
fi

# csv module, not cut: hotspot_res contains commas and is therefore quoted.
read -r TARGET_REL HOTSPOT < <(python3 - "$CSV" "$PDB_ID" <<'PY'
import csv, sys
path, pdb_id = sys.argv[1], sys.argv[2]
for r in csv.DictReader(open(path)):
    if r["PDB_ID"] == pdb_id:
        print(r["target_pdb"], r["hotspot_res"])
        break
else:
    sys.exit("no row for %s in %s" % (pdb_id, path))
PY
) || exit 1

TARGET=$IN/$TARGET_REL
if [ ! -f "$TARGET" ]; then
  echo "target missing: $TARGET -- run Step 2 first" >&2
  exit 1
fi

mkdir -p "$OUT" "$REPO/runs/rfab"

args=(
  --config-name antibody
  "antibody.target_pdb=$TARGET"
  "antibody.framework_pdb=$FRAMEWORK"
  "antibody.design_loops=[H1:,H2:,H3:5-13]"
  "inference.num_designs=$NDES"
  "inference.output_prefix=$OUT/nb"
)
[ "$ARM" = hotspot ] && args+=("ppi.hotspot_res=$HOTSPOT")

echo "=== $PDB_ID $ARM n=$NDES"
echo "    target    $TARGET"
echo "    framework $FRAMEWORK"
[ "$ARM" = hotspot ] && echo "    hotspots  $HOTSPOT"

start=$(date +%s)
$RFAB_RUN "$RFAB/scripts/rfdiffusion_inference.py" "${args[@]}" 2>&1 | tee "$OUT/rfdiff.log"
rc=${PIPESTATUS[0]}
end=$(date +%s)

n=$(ls "$OUT"/*.pdb 2>/dev/null | wc -l)

if [ ! -s "$MANIFEST" ]; then
  echo "utc,pdb_id,arm,num_designs,backbones_written,wall_seconds,exit_code,hotspot,job_id" > "$MANIFEST"
fi
printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$PDB_ID" "$ARM" "$NDES" "$n" \
  "$((end - start))" "$rc" "${HOTSPOT//,/;}" "${SLURM_JOB_ID:-interactive}" \
  >> "$MANIFEST"

echo "--- $PDB_ID $ARM: $n backbones, $((end - start))s, exit $rc"
exit "$rc"
