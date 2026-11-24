#!/usr/bin/env bash
# Fold every job in runs/boltz_in with Boltz, on one EC2 GPU instance.
#
# Converted from psc_boltz.sbatch when PSC access ended. The Slurm directives are
# the only part that died; everything else here was earned on PSC and still earns
# its keep on EC2, for reasons noted at each step.
#
# Boltz is the EXCLUDED FINAL SCORER: it never touches generation or candidate
# selection, and it folds from sequence only, with no template and no designed
# backbone. Both are requirements from the proposal, not stylistic choices.
#
#   cd ~/gen_nanobody/runs/boltz_in && nohup bash ../../src/evaluate/run_boltz.sh > ../fold.log 2>&1 &
#
# Run it from the INPUT ROOT (the dir holding <target>_H_T/ and msa/), because the
# manifests reference MSAs by relative path so the tree can be rsynced anywhere.

set -uo pipefail          # not -e: one bad fold must not kill the batch

IN_ROOT=${IN_ROOT:-.}
OUT=${OUT:-../boltz_out}
MSA=${MSA:-cached}        # cached | server | empty  (must match make_fold_inputs.py)
RECYCLE=${RECYCLE:-3}
SAMPLES=${SAMPLES:-1}
ONLY=${ONLY:-}            # substring filter, e.g. ONLY=pdb_00008qf4 for one target

cd "$IN_ROOT" || { echo "no such dir: $IN_ROOT"; exit 1; }
mkdir -p "$OUT"

MSA_FLAG=""
[ "$MSA" = "server" ] && MSA_FLAG="--use_msa_server"

# Fail before the GPU is spent, not 90s into each fold. On PSC this caught a CUDA 13
# torch wheel on a V100 (sm_70 dropped), which surfaced only AFTER the MSA fetch, once
# per fold. The same late failure exists on EC2 with a mismatched wheel; an A10G is
# sm_86. Cheap insurance against burning a whole batch.
python3 - <<'PY' || { echo "ABORT: torch cannot run on this GPU"; exit 2; }
import sys, torch
if not torch.cuda.is_available():
    print("torch sees no CUDA device"); sys.exit(1)
cap = torch.cuda.get_device_capability()
sm, archs = "sm_%d%d" % cap, torch.cuda.get_arch_list()
print("gpu         %s (%s)" % (torch.cuda.get_device_name(0), sm))
print("torch       %s, built for %s" % (torch.__version__, " ".join(archs)))
if sm not in archs:
    # Not fatal on its own: cubins are forward compatible across minor versions
    # within a major, so sm_86 code runs on an sm_89 L40S. The matmul decides.
    print("note: %s absent from the arch list; relying on minor-version compat" % sm)
try:
    x = torch.randn(256, 256, device="cuda")
    float((x @ x).sum().cpu())
except Exception as e:
    print("GPU UNUSABLE: %s: %s" % (type(e).__name__, e))
    print("  %s vs torch arch list %s" % (sm, " ".join(archs)))
    sys.exit(1)
print("gpu check   ok")
PY

echo "host        $(hostname)"
echo "boltz       $(pip show boltz 2>/dev/null | awk '/^Version:/{print $2}')"
echo "cache       ${BOLTZ_CACHE:-$HOME/.boltz}"
echo "msa mode    $MSA"
echo "inputs      $(pwd)"
echo "outputs     $OUT"
echo

# Fold in manifest order, which puts the native control first for every target. If the
# predictor cannot reproduce a known interface, the design numbers are not interpretable
# anyway, and we want to know that in minute five rather than hour nine.
jobs=()
for man in */manifest.csv; do
  tgt=$(dirname "$man")
  while IFS=, read -r job _rest; do
    [ "$job" = "job" ] && continue
    [ -n "$ONLY" ] && [[ "$tgt/$job" != *"$ONLY"* ]] && continue
    [ -f "$tgt/$job.yaml" ] || { echo "WARN  $tgt/$job.yaml in manifest but not on disk"; continue; }
    jobs+=("$tgt/$job")
  done < "$man"
done
[ ${#jobs[@]} -gt 0 ] || { echo "no jobs found under $(pwd)"; exit 1; }

if [ "$MSA" = "cached" ]; then
  missing=$(grep -ho 'msa: msa/[^ ]*' */*.yaml | sort -u | sed 's/msa: //' | while read -r f; do
    [ -f "$f" ] || echo "$f"; done)
  if [ -n "$missing" ]; then
    echo "ABORT: MSA=cached but these MSA files are missing:"; echo "$missing"
    echo "harvest them once per antigen (fold each native with --use_msa_server and copy"
    echo "boltz_results_native/msa/native_1.csv), or rerun with MSA=server"; exit 3
  fi
fi

echo "${#jobs[@]} jobs queued"
echo

done_n=0; fail=0; skip=0
t_all=$SECONDS
for j in "${jobs[@]}"; do
  dest="$OUT/$j"
  # Resume, not restart. On Slurm this guarded against the 2h walltime kill; here it
  # guards against a dropped ssh session, an overnight interruption or a spot reclaim.
  # find, not a ** glob: globstar is off by default, so "**" would match one level
  # only and never see boltz_results_*/predictions/*/.
  if [ -d "$dest" ] && [ -n "$(find "$dest" -name '*model_0.pdb' -print -quit 2>/dev/null)" ]; then
    skip=$((skip + 1)); continue
  fi
  mkdir -p "$dest"
  echo "=== $j ==="
  t0=$SECONDS
  # --output_format pdb so fold_metrics.py reuses the plain PDB reader the rest of
  # src/evaluate already has; no mmCIF dependency anywhere in the project.
  if boltz predict "$j.yaml" \
        --out_dir "$dest" \
        --output_format pdb \
        --diffusion_samples "$SAMPLES" \
        --recycling_steps "$RECYCLE" \
        $MSA_FLAG > "$dest/boltz.log" 2>&1; then
    done_n=$((done_n + 1))
    echo "--- ok in $((SECONDS - t0))s   [$done_n done, $fail failed, $skip skipped]"
  else
    fail=$((fail + 1))
    echo "--- FAILED after $((SECONDS - t0))s, see $dest/boltz.log"
    tail -5 "$dest/boltz.log" | sed 's/^/      /'
  fi
done

echo
echo "$done_n folded, $fail failed, $skip already done, of ${#jobs[@]} in $((SECONDS - t_all))s"
echo "next: rsync $OUT back, then per target:"
echo "  python3 src/evaluate/fold_metrics.py report --target <pdb_id>"
exit $((fail > 0))
