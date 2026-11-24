#!/usr/bin/env bash
# Fold the end-to-end target with Boltz on PSC Bridges-2.
#
# SUPERSEDED Nov 16 2026: PSC access ended. The live runner is src/evaluate/run_boltz.sh,
# which carries over the GPU preflight, --output_format pdb, native-first ordering and
# the resume check. This file is kept only as the record of the run config we had on
# PSC; it has no working submission path any more.
#
# Boltz is the EXCLUDED FINAL SCORER: it never touches generation or candidate
# selection, and it folds from sequence only, with no template and no designed
# backbone. Both are requirements from the proposal, not stylistic choices.
#
#   sbatch --export=ALL,IN=$PROJECT/gen_nanobody/boltz_in/pdb_00008qf4_H_T,OUT=$PROJECT/gen_nanobody/boltz_out/pdb_00008qf4_H_T \
#          src/evaluate/psc_boltz.sbatch
#
# Account cis260266p is inferred from the PSC paths C recorded in targets.csv
# (/ocean/projects/cis260266p/cyang28). Confirm with `projects` before submitting.
#
#SBATCH -J boltz_nb
#SBATCH -p GPU-shared
#SBATCH -N 1
#SBATCH --gres=gpu:h100-80:1
#SBATCH --ntasks-per-node=5
#SBATCH -t 02:00:00
#SBATCH -o boltz_%j.out
#SBATCH -e boltz_%j.err
#SBATCH -A cis260266p

set -euo pipefail

IN=${IN:?set IN=<dir of .yaml files>}
OUT=${OUT:?set OUT=<output dir>}
MSA=${MSA:-server}            # server | empty
ENVNAME=${ENVNAME:-boltz}
RECYCLE=${RECYCLE:-3}
SAMPLES=${SAMPLES:-1}

module load anaconda3
# shellcheck disable=SC1091
source activate "$ENVNAME"

# Weights are multi-GB and /jet/home is small, so keep the cache on /ocean
# unless the caller already pointed BOLTZ_CACHE somewhere.
export BOLTZ_CACHE=${BOLTZ_CACHE:-$PROJECT/.boltz}
mkdir -p "$BOLTZ_CACHE" "$OUT"

MSA_FLAG="--use_msa_server"
[ "$MSA" = "empty" ] && MSA_FLAG=""

# FIX_TORCH=1 installs a build that still carries Volta (sm_70) cubins. The
# default wheel is a CUDA 13 build, and CUDA 13 dropped Volta, so on a V100 it
# fails only after the MSA, once per fold.
if [ "${FIX_TORCH:-0}" = "1" ]; then
  echo "installing torch ${TORCH_VER:-2.5.1} (cu121, has sm_70)..."
  pip install --quiet --index-url https://download.pytorch.org/whl/cu121 \
      "torch==${TORCH_VER:-2.5.1}" || { echo "ABORT: torch install failed"; exit 3; }
fi

# Fail before the GPU is spent, not 90s into each fold. torch ships cubins for
# a fixed set of architectures; a V100 is sm_70 and CUDA 13 builds dropped it,
# which surfaces only as "no kernel image is available for execution on the
# device" after the MSA has already been fetched.
python - <<'PY' || { echo "ABORT: torch cannot run on this GPU"; exit 2; }
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
# the arch list can also lie if the install is half-broken, so actually use the GPU
try:
    x = torch.randn(256, 256, device="cuda")
    float((x @ x).sum().cpu())
except Exception as e:
    print("GPU UNUSABLE: %s: %s" % (type(e).__name__, e))
    print("  %s vs torch arch list %s" % (sm, " ".join(archs)))
    print("  pick a supported GPU (--gres=gpu:h100-80:1), or resubmit with FIX_TORCH=1")
    sys.exit(1)
print("gpu check   ok")
PY

echo "host        $(hostname)"
echo "boltz       $(pip show boltz 2>/dev/null | awk '/^Version:/{print $2}')"
echo "cache       $BOLTZ_CACHE"
echo "msa mode    $MSA"
echo "inputs      $IN"
echo

shopt -s nullglob
# native first: it is the control, and if the predictor cannot reproduce the
# known interface the design numbers are not interpretable anyway.
yamls=()
[ -f "$IN/native.yaml" ] && yamls+=("$IN/native.yaml")
for y in "$IN"/*.yaml; do
  [ "$(basename "$y")" = "native.yaml" ] || yamls+=("$y")
done
[ ${#yamls[@]} -gt 0 ] || { echo "no .yaml in $IN"; exit 1; }

fail=0
for y in "${yamls[@]}"; do
  job=$(basename "$y" .yaml)
  dest="$OUT/$job"
  # find, not a ** glob: globstar is off by default, so "**" would match one
  # level only and never see boltz_results_*/predictions/*/
  if [ -d "$dest" ] && [ -n "$(find "$dest" -name '*model_0.pdb' -print -quit 2>/dev/null)" ]; then
    echo "skip  $job (already done)"
    continue
  fi
  mkdir -p "$dest"
  echo "=== $job ==="
  t0=$SECONDS
  # --output_format pdb so fold_metrics.py can reuse the plain PDB reader the
  # rest of src/evaluate already uses; no mmCIF dependency anywhere.
  if boltz predict "$y" \
        --out_dir "$dest" \
        --output_format pdb \
        --diffusion_samples "$SAMPLES" \
        --recycling_steps "$RECYCLE" \
        $MSA_FLAG; then
    echo "--- $job ok in $((SECONDS - t0))s"
  else
    echo "--- $job FAILED after $((SECONDS - t0))s"
    fail=$((fail + 1))
  fi
done

echo
echo "$fail failed of ${#yamls[@]}"
echo "next: rsync the outputs back, then"
echo "  python3 src/evaluate/fold_metrics.py report --target $(basename "$IN" _H_T)"
exit $fail
