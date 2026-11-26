# Runbook — A (generation): RFantibody on the 21-target benchmark

Owner: Aidan. Everything below runs from the repo root on PSC. Nothing here needs a
file opened in an editor: each block either writes a file or runs one.

`main` already carries the benchmark C built (`data/processed/prepped_v2/`, the epitope
and antigen-mapping tables, QC) and B's finished ProteinMPNN lane (`src/evaluate/`,
`results/cdr_recovery.csv`). `src/generate/` is empty — Step 1 fills it.

Two things RFantibody needs are not on `main` yet, and Step 3 produces both:

1. **A target-only PDB.** `prepped_v2/*.pdb` carries the native nanobody in chain `H`
   next to the antigen in chain `T`. RFdiffusion generates a new `H`, so the native one
   must not be in the file handed to it.
2. **`ppi.hotspot_res` in the target file's numbering.** `epitope_residues_v2.csv` is in
   *original* antigen numbering (chain `A`, residue 211). The crop renumbered the merged
   antigen to `T1..Tn`. `antigen_residue_mapping_v2.csv` holds the translation
   (`output_res_id`), but nothing on `main` applies it.

You do **not** need ANARCI, and you do **not** need to wait for `data/processed/hlt/`.
`src/evaluate/make_hlt.py` exists because B's CDR-recovery scoring had to be moved from
C's IMGT spans onto the Chothia spans RFantibody uses. That only affects a *framework*
PDB. You use RFantibody's shipped framework, which is already Chothia-numbered with
`REMARK PDBinfo-LABEL` lines. Step 2 verifies that rather than assuming it.

---

## Step 0 — pull and point at your install

```bash
cd $PROJECT/gen_nanobody
git pull

export REPO=$PWD
export RFAB=$PROJECT/RFantibody          # wherever your install is
export RFAB_RUN="python"                 # if you run RFantibody in a container, instead:
# export RFAB_RUN="apptainer exec --nv $PROJECT/rfantibody.sif poetry run python"
```

`RFAB_RUN` is however you invoked `rfdiffusion_inference.py` in your pilot. Everything
else is derived. Keep this shell: Step 5 submits from it and Slurm inherits these.

---

## Step 1 — write the five scripts

One paste. Creates `src/generate/`, nothing else.

```bash
mkdir -p "$REPO/src/generate" "$REPO/logs"
cd "$REPO"
```

### `src/generate/00_check_env.sh`

```bash
cat > src/generate/00_check_env.sh <<'EOF_00'
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
EOF_00
```

### `src/generate/01_make_rfab_targets.py`

```bash
cat > src/generate/01_make_rfab_targets.py <<'EOF_01'
#!/usr/bin/env python3
"""Build RFdiffusion target inputs from C's benchmark.

Two things RFantibody needs that are not on main:

  1. A target-only PDB. C's prepped_v2 files carry the native nanobody in chain
     H alongside the antigen in chain T. RFdiffusion generates a new H, so the
     native one must not be in the target file.
  2. ppi.hotspot_res in the numbering of that target file. C's
     epitope_residues_v2.csv is in ORIGINAL antigen numbering (chain A, residue
     211); the crop renumbered the merged antigen to T1..Tn. The translation
     lives in antigen_residue_mapping_v2.csv (output_res_id).

Hotspots are a subset, not the whole epitope: a 25-residue epitope passed whole
over-constrains the dock, and RFantibody's own examples use 3-6. Selection is
deterministic -- seed with the most buried epitope residue (smallest
min_distance_to_nanobody_A), then greedily add the next most buried residue
whose CA is at least HOTSPOT_MIN_SEP from every residue already chosen, so the
hotspots spread across the interface instead of clustering in one pocket. The
full epitope is still written to the epitope_all_T column.

stdlib only, so it runs under any python3 on PSC without an env.
"""
import argparse
import csv
import math
import os

HOTSPOT_MIN_SEP = 8.0


def parse_t_ca(pdb_path):
    coords = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[21] != "T" or line[12:16].strip() != "CA":
                continue
            coords[int(line[22:26])] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def dist(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def pick_hotspots(cand, coords, n):
    """cand is epitope residue ids in T numbering, most buried first."""
    chosen = [cand[0]]
    for resid in cand[1:]:
        if len(chosen) >= n:
            break
        if resid not in coords:
            continue
        seps = [dist(coords[resid], coords[c]) for c in chosen if c in coords]
        if not seps or min(seps) >= HOTSPOT_MIN_SEP:
            chosen.append(resid)
    # If the separation rule could not fill n slots, top up by burial order
    # rather than returning fewer hotspots than asked for.
    for resid in cand[1:]:
        if len(chosen) >= n:
            break
        if resid not in chosen:
            chosen.append(resid)
    return sorted(chosen[:n])


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed", default="data/processed")
    ap.add_argument("--out-dir", default="data/processed/rfab_inputs")
    ap.add_argument("--n-hotspots", type=int, default=4)
    args = ap.parse_args()

    proc = args.processed
    epi = read_rows(os.path.join(proc, "epitope_residues_v2.csv"))
    amap = read_rows(os.path.join(proc, "antigen_residue_mapping_v2.csv"))
    prep = {r["PDB_ID"]: r for r in read_rows(os.path.join(proc, "prepped_v2.csv"))}
    qc = {r["PDB_ID"]: r for r in read_rows(os.path.join(proc, "structural_qc_final.csv"))}

    def key(r):
        return (r["PDB_ID"], r["original_chain"], r["original_res_id"],
                r["original_insertion_code"] or "")

    to_out = {key(r): int(r["output_res_id"]) for r in amap}

    by_pdb = {}
    unmapped = 0
    for r in epi:
        out_id = to_out.get(key(r))
        if out_id is None:
            unmapped += 1
            continue
        by_pdb.setdefault(r["PDB_ID"], []).append(
            (float(r["min_distance_to_nanobody_A"]), out_id))

    os.makedirs(os.path.join(args.out_dir, "targets"), exist_ok=True)
    rows = []
    for pdb_id in sorted(by_pdb):
        src = os.path.join(proc, "prepped_v2", pdb_id + "_H_T.pdb")
        coords = parse_t_ca(src)

        dst = os.path.join(args.out_dir, "targets", pdb_id + "_T.pdb")
        with open(src) as fi, open(dst, "w") as fo:
            for line in fi:
                if line.startswith(("ATOM", "HETATM")) and line[21] != "T":
                    continue
                if line.startswith(("CRYST1", "ATOM", "HETATM", "TER", "END")):
                    fo.write(line)

        cand = [rid for _, rid in sorted(by_pdb[pdb_id])]
        hs = pick_hotspots(cand, coords, args.n_hotspots)
        q = qc.get(pdb_id, {})
        rows.append({
            "PDB_ID": pdb_id,
            "target_pdb": os.path.relpath(dst, args.out_dir),
            "T_residues": len(coords),
            "native_H_residues": prep.get(pdb_id, {}).get("nanobody_residues", ""),
            "epitope_residues": len(cand),
            "hotspot_res": "[" + ",".join("T%d" % i for i in hs) + "]",
            "epitope_all_T": " ".join("T%d" % i for i in cand),
            "QC": q.get("QC", ""),
            "QC_issues": q.get("issues", ""),
        })

    rows.sort(key=lambda r: r["T_residues"])
    for i, r in enumerate(rows, start=1):
        r["size_rank"] = i

    cols = ["size_rank", "PDB_ID", "target_pdb", "T_residues", "native_H_residues",
            "epitope_residues", "hotspot_res", "epitope_all_T", "QC", "QC_issues"]
    out_csv = os.path.join(args.out_dir, "rfab_targets.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    print("unmapped epitope residues: %d" % unmapped)
    print("targets written: %d" % len(rows))
    print()
    print("%4s %-14s %5s %5s  %-24s %s" %
          ("rank", "PDB", "T", "epi", "hotspot_res", "QC"))
    for r in rows:
        print("%4d %-14s %5d %5d  %-24s %s" %
              (r["size_rank"], r["PDB_ID"], r["T_residues"],
               r["epitope_residues"], r["hotspot_res"], r["QC"]))
    print()
    print("wrote %s" % out_csv)
    if unmapped:
        raise SystemExit("STOP: %d epitope residues did not map to a T residue id"
                         % unmapped)


if __name__ == "__main__":
    main()
EOF_01
```

### `src/generate/02_rfdiff_one.sh`

```bash
cat > src/generate/02_rfdiff_one.sh <<'EOF_02'
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
EOF_02
```

### `src/generate/03_rfdiff_batch.sbatch`

```bash
cat > src/generate/03_rfdiff_batch.sbatch <<'EOF_03'
#!/bin/bash
#SBATCH -J rfab
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH -t 08:00:00
#SBATCH -o logs/rfab_%A_%a.out
#SBATCH --array=0-20
#
# One array task per target, both arms, ordered smallest antigen first so the
# small ones land early. Submit with the allocation on the command line:
#   sbatch -A <allocation> src/generate/03_rfdiff_batch.sbatch
# Re-runnable: 02_rfdiff_one.sh skips finished targets, so requeueing after a
# timeout resumes.
set -uo pipefail

REPO=${REPO:-$SLURM_SUBMIT_DIR}
cd "$REPO" || exit 1
mkdir -p logs

# Carried from the submitting shell by --export=ALL (the default). If they are
# unset here, submit again from a shell where Step 1 was run.
: "${RFAB:?RFAB unset in the job environment}"
RFAB_RUN=${RFAB_RUN:-python}
NDES=${NDES:-20}

mapfile -t IDS < <(python3 - data/processed/rfab_inputs/rfab_targets.csv <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
rows.sort(key=lambda r: int(r["T_residues"]))
for r in rows:
    print(r["PDB_ID"])
PY
)

if [ "$SLURM_ARRAY_TASK_ID" -ge "${#IDS[@]}" ]; then
  echo "task $SLURM_ARRAY_TASK_ID beyond ${#IDS[@]} targets, nothing to do"
  exit 0
fi

ID=${IDS[$SLURM_ARRAY_TASK_ID]}
echo "task $SLURM_ARRAY_TASK_ID -> $ID"

rc=0
bash src/generate/02_rfdiff_one.sh "$ID" hotspot  "$NDES" || rc=$?
bash src/generate/02_rfdiff_one.sh "$ID" unguided "$NDES" || rc=$?
exit "$rc"
EOF_03
```

### `src/generate/04_report_status.py`

```bash
cat > src/generate/04_report_status.py <<'EOF_04'
#!/usr/bin/env python3
"""Generation status for the midterm report, read off runs/rfab/manifest.csv.

Writes results/generation_status.csv and prints the paragraph numbers. Counts
backbones on disk rather than trusting num_designs, so a partial run reports
what it actually produced. stdlib only.
"""
import argparse
import csv
import glob
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", default="results/generation_status.csv")
    args = ap.parse_args()

    man = os.path.join(args.repo, "runs", "rfab", "manifest.csv")
    if not os.path.exists(man):
        raise SystemExit("no manifest at %s -- nothing has run yet" % man)

    rows = list(csv.DictReader(open(man)))
    # Keep the last attempt per (target, arm); earlier ones were skipped or failed.
    latest = {}
    for r in rows:
        latest[(r["pdb_id"], r["arm"])] = r

    out_rows = []
    for (pid, arm), r in sorted(latest.items()):
        d = os.path.join(args.repo, "runs", "rfab", "%s_%s" % (pid, arm))
        on_disk = len([p for p in glob.glob(os.path.join(d, "*.pdb"))])
        out_rows.append({
            "pdb_id": pid,
            "arm": arm,
            "requested": int(r["num_designs"]),
            "backbones_on_disk": on_disk,
            "wall_seconds": int(r["wall_seconds"]),
            "exit_code": int(r["exit_code"]),
            "job_id": r["job_id"],
            "utc": r["utc"],
        })

    os.makedirs(os.path.dirname(os.path.join(args.repo, args.out)) or ".", exist_ok=True)
    cols = ["pdb_id", "arm", "requested", "backbones_on_disk", "wall_seconds",
            "exit_code", "job_id", "utc"]
    with open(os.path.join(args.repo, args.out), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)

    ok = [r for r in out_rows if r["exit_code"] == 0 and r["backbones_on_disk"] > 0]
    tot_bb = sum(r["backbones_on_disk"] for r in ok)
    tot_s = sum(r["wall_seconds"] for r in ok)
    targets = sorted({r["pdb_id"] for r in ok})
    arms = sorted({r["arm"] for r in ok})
    failed = [r for r in out_rows if r["exit_code"] != 0 or r["backbones_on_disk"] == 0]

    print("runs recorded          %d" % len(out_rows))
    print("runs with backbones    %d" % len(ok))
    print("targets with backbones %d / 21" % len(targets))
    print("arms represented       %s" % (", ".join(arms) or "none"))
    print("backbones on disk      %d" % tot_bb)
    print("GPU wall time          %.2f h (%d s)" % (tot_s / 3600.0, tot_s))
    if tot_bb:
        print("mean per backbone      %.1f s" % (tot_s / tot_bb))
    if failed:
        print("\nfailed or empty:")
        for r in failed:
            print("  %s %-8s exit %s, %d pdb" %
                  (r["pdb_id"], r["arm"], r["exit_code"], r["backbones_on_disk"]))
    print("\nwrote %s" % os.path.join(args.out))


if __name__ == "__main__":
    main()
EOF_04
```

```bash
chmod +x src/generate/*.sh
git add src/generate && git commit -m "A: RFantibody target prep, run wrapper, batch array, status reporter"
```

---

## Step 2 — check the environment before spending GPU time

```bash
bash src/generate/00_check_env.sh
```

Expect `ALL CHECKS PASSED`. It verifies the inference script exists, that the framework
PDB carries `REMARK PDBinfo-LABEL` lines (RFantibody silently misbehaves on a framework
without them), that all four of C's tables are present, and that `prepped_v2/` holds 21
PDBs. Any `FAIL` line names what to fix. Do not continue past a failure.

---

## Step 3 — build the targets and hotspots

CPU, under a minute, login node is fine.

```bash
python3 src/generate/01_make_rfab_targets.py
```

The first line must read `unmapped epitope residues: 0`. **If it is not 0, stop** — the
epitope table and the mapping table disagree and every hotspot downstream is wrong. The
script exits non-zero in that case.

It writes `data/processed/rfab_inputs/targets/*_T.pdb` (21 target-only PDBs) and
`rfab_inputs/rfab_targets.csv`, and prints the table it just built. The expected result:

```
rank PDB                T   epi  hotspot_res              QC
   1 pdb_00009mnz      47    11  [T31,T38,T43,T47]        PASS
   2 pdb_00009h77      57    17  [T2,T24,T35,T42]         PASS
   3 pdb_00008qf4      67    25  [T3,T9,T16,T35]          PASS
   4 pdb_00009g2g      88    12  [T30,T35,T39,T58]        PASS
   5 pdb_00009min     100    22  [T45,T82,T89,T95]        PASS
   6 pdb_00009npt     111    25  [T16,T70,T74,T77]        PASS
   7 pdb_00009yc3     115    17  [T4,T8,T11,T109]         PASS
   8 pdb_00009g3e     130    18  [T71,T75,T83,T111]       PASS
   9 pdb_00008idi     131    18  [T2,T6,T82,T85]          PASS
  10 pdb_00009ty2     140    19  [T29,T56,T100,T108]      PASS
  11 pdb_00009e7l     142    15  [T81,T86,T113,T140]      PASS
  12 pdb_00009b71     145    17  [T19,T44,T58,T99]        PASS
  13 pdb_00009rye     158    21  [T26,T31,T136,T141]      PASS
  14 pdb_00009my8     161    25  [T4,T15,T143,T152]       PASS
  15 pdb_00008rwf     164    30  [T54,T100,T109,T112]     PASS
  16 pdb_00009ljj     169    21  [T14,T65,T86,T91]        PASS
  17 pdb_00009rxi     174    22  [T25,T74,T101,T121]      PASS
  18 pdb_00008yj8     190    19  [T14,T97,T171,T182]      CHECK
  19 pdb_00009rmk     225    28  [T31,T35,T111,T218]      PASS
  20 pdb_00009q6n     240    29  [T72,T86,T161,T223]      PASS
  21 pdb_00009dc3     300    35  [T132,T181,T255,T259]    PASS
```

Four hotspots, not the whole epitope: a 25-residue epitope passed whole over-constrains
the dock, and RFantibody's own examples use 3–6. Selection is deterministic — seed with
the most buried epitope residue, then greedily add the next most buried residue whose CA
is ≥ 8 Å from every residue already chosen, so the hotspots spread across the interface
instead of clustering in one pocket. The full epitope is kept in the `epitope_all_T`
column if you want a different policy later.

```bash
git add data/processed/rfab_inputs && git commit -m "A: RFdiffusion target PDBs and hotspot table"
```

---

## Step 4 — pilot on the three smallest targets

```bash
interact -p GPU-shared --gres=gpu:v100-32:1 -t 01:00:00
cd $PROJECT/gen_nanobody
export REPO=$PWD RFAB=$PROJECT/RFantibody RFAB_RUN="python"   # re-export in the new shell

bash src/generate/02_rfdiff_one.sh pdb_00009mnz hotspot 5
```

`num_designs=5` first — confirm it runs before scaling. Then the other two and the
control arm:

```bash
bash src/generate/02_rfdiff_one.sh pdb_00009h77 hotspot 5
bash src/generate/02_rfdiff_one.sh pdb_00008qf4 hotspot 5
bash src/generate/02_rfdiff_one.sh pdb_00009mnz unguided 5
nvidia-smi --query-gpu=memory.used --format=csv
```

The unguided arm is the same command with `unguided`; the script drops `ppi.hotspot_res`
and nothing else changes. That is the matched control the proposal's RQ3 needs.

Record the peak memory and the per-design wall time from the `---` summary lines. Those
two numbers are what let you state in the report how much of the sweep fits before the
allocation ends.

Each run appends a row to `runs/rfab/manifest.csv`. A target that already has backbones
on disk is skipped, so re-running is safe.

---

## Step 5 — submit the full sweep

From the login-node shell of Step 0, not from inside `interact`:

```bash
sbatch -A <your_allocation_id> src/generate/03_rfdiff_batch.sbatch
squeue -u $USER
```

21 array tasks, one per target, ordered smallest antigen first, both arms per task,
`num_designs=20` each. Override with `NDES=40 sbatch ...` if the pilot timings say there
is room. The array is re-runnable: a requeue after a timeout resumes rather than redoing
finished targets.

`pdb_00009dc3` (T=300, plus a ~125-residue framework ≈ 425 residues) is the last task
because attention cost grows with the square of length and it is the heaviest in the set.
Leave it last rather than promoting it.

---

## Step 6 — the number that goes in the report

```bash
python3 src/generate/04_report_status.py
git add results/generation_status.csv src/generate && git commit -m "A: generation status"
```

Prints targets covered, backbones on disk, GPU wall time, mean seconds per backbone, and
anything that failed or came back empty; writes `results/generation_status.csv`.

It counts PDBs on disk rather than trusting `num_designs`, so a partially finished sweep
reports what it actually produced.

For the midterm, report what is on disk plus the queued job IDs — not a projection. A
real small number is worth more than a planned large one.

---

## Four things that will otherwise cost you time

**Use `prepped_v2/`, not `prepped_mpnn/`.** Step 3 already reads the right one. The `T`
numbering happens to be identical in both, so this would not actually break — but
`prepped_mpnn/` exists to work around a ProteinMPNN insertion-code parser bug and has
nothing to do with targets.

**`8yj8` is flagged `CHECK`, and you should still run it.** The 2.194 Å C–N anomaly is in
chain `H` (`H_CN_anomalies=1`, `T_CN_anomalies=0`), and chain `H` is discarded for
generation — you are producing a new one. Its target is clean. Reading
`structural_qc_final.csv` without checking which chain the anomaly is in is the easy way
to drop a good target.

**Use RFantibody's shipped framework, not the native nanobody.** The proposal commits to
withholding "the reference nanobody's CDR sequences, coordinates, and binding pose." A
native framework leaks the answer and is a different experiment.

**Submit with `sbatch`, not `interact`.** Interactive sessions die with the terminal and
the allocation ends Monday Nov 9. Anything not queued as a batch job by Sunday is not in
the final report.
