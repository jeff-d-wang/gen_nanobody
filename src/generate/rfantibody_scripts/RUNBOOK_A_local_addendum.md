# Addendum — running the generation lane locally (no PSC)

Read this next to `RUNBOOK_A_generation.md`. PSC is down, so Slurm is out. Everything
else in that runbook still applies; this file is only the delta.

## What changes

| Runbook step | Local |
|---|---|
| Step 0, pull and set `RFAB` / `RFAB_RUN` | same, but point at your local install |
| Step 1, write the scripts | skip if they are already on `main` |
| Step 2, `00_check_env.sh` | same, runs anywhere |
| Step 3, `01_make_rfab_targets.py` | same — CPU, stdlib, seconds |
| Step 4, pilot via `interact` | drop `interact`, run `02_rfdiff_one.sh` directly |
| **Step 5, `sbatch 03_rfdiff_batch.sbatch`** | **replaced by `03_rfdiff_local.sh`** |
| Step 6, `04_report_status.py` | same |

`02_rfdiff_one.sh` and `04_report_status.py` need no change: neither assumes Slurm, and
the manifest records `interactive` as the job id when `SLURM_JOB_ID` is unset.

## Do Steps 2 and 3 first, before touching the GPU

```bash
export REPO=$PWD
bash src/generate/00_check_env.sh
python3 src/generate/01_make_rfab_targets.py
git add data/processed/rfab_inputs && git commit -m "A: RFdiffusion target PDBs and hotspot table" && git push
```

No GPU involved. This produces the 21 target-only PDBs and the hotspot table, which is a
committed deliverable on its own and unblocks everything downstream — do it even if the
generation itself does not run tonight.

## The new driver

```bash
cat > src/generate/03_rfdiff_local.sh <<'EOF_LOCAL'
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
EOF_LOCAL
chmod +x src/generate/03_rfdiff_local.sh
```

Then start small and let it climb:

```bash
# one target first, to find out what your card can actually hold
bash src/generate/02_rfdiff_one.sh pdb_00009mnz hotspot 2

# then sweep, smallest antigen first, capped by antigen size
bash src/generate/03_rfdiff_local.sh 120 10
```

Arguments are `[max_T_residues] [num_designs] [arms]`. `03_rfdiff_local.sh 120 10` runs
only targets whose cropped antigen is at most 120 residues, 10 designs each, both arms.

It prints your GPU and its memory first, runs smallest antigen first, and **keeps going
past a failure**. On a consumer card the large targets are expected to run out of memory,
and stopping at the first one would throw away the small targets that did fit. Every
attempt lands in `runs/rfab/manifest.csv` with its exit code either way.

## Memory, concretely

RFdiffusion's cost grows with the square of total length, and total length here is the
cropped antigen plus the ~125-residue framework:

| rank | target | T | total ≈ |
|---:|---|---:|---:|
| 1 | 9mnz | 47 | 172 |
| 3 | 8qf4 | 67 | 192 |
| 6 | 9npt | 111 | 236 |
| 10 | 9ty2 | 140 | 265 |
| 15 | 8rwf | 164 | 289 |
| 21 | 9dc3 | 300 | 425 |

Do not guess where your card gives out — run rank 1, then raise `max_T` until it stops.
Where it stops is a number worth reporting, not a failure: RQ3 is candidate yield against
compute, and "hotspot-conditioned generation succeeded on targets up to N residues on a
<card>" is a measurement of exactly that axis.

## If RFantibody is not installed locally

The pilot ran on PSC, so a local install may not exist. RFantibody ships as a CUDA
container; standing that up is a multi-hour job, not a quick one, and it is the wrong
thing to start at this hour.

In that case tonight's generation deliverable is Steps 2 and 3 plus whatever the manifest
already holds, and the report says plainly that generation is blocked on compute with the
allocation expired. A stated constraint reads better than an empty figure.

Any backbones produced on PSC before it went down are on PSC disk and unreachable until
it returns. `04_report_status.py` reports only what is on the machine it runs on.
