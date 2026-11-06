#!/usr/bin/env bash
# ProteinMPNN CDR redesign over C's prepped targets -> results/cdr_recovery.csv
#
#   src/evaluate/run_mpnn_cdr.sh data/processed/hlt
#
# CDRs come from the remarks chothia2HLT.py writes into the HLT files, the same
# definition A's RFdiffusion/proteinmpnn uses, so recovery scores the real
# pipeline stage.
#
# Env overrides: NSEQ (16) TEMP (0.1) SEED (37) HEAVY (H) RUN (runs/mpnn_<seed>)
# Runs on CPU if no GPU; ProteinMPNN picks the device itself.
set -euo pipefail

HLT_DIR=${1:?usage: $0 <hlt_dir>}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MPNN=$REPO/third_party/ProteinMPNN
NSEQ=${NSEQ:-16}; TEMP=${TEMP:-0.1}; SEED=${SEED:-37}; HEAVY=${HEAVY:-H}
RUN=${RUN:-$REPO/runs/mpnn_seed$SEED}

[ -d "$MPNN" ] || { echo "missing $MPNN - git clone https://github.com/dauparas/ProteinMPNN.git third_party/ProteinMPNN"; exit 1; }
mkdir -p "$RUN"

# 1. ProteinMPNN parses the PDBs. Everything downstream indexes off THIS file.
python3 "$MPNN/helper_scripts/parse_multiple_chains.py" \
  --input_path "$HLT_DIR" --output_path "$RUN/parsed.jsonl"

# 2. Designed chain = nanobody only; antigen stays visible context.
python3 "$MPNN/helper_scripts/assign_fixed_chains.py" \
  --input_path "$RUN/parsed.jsonl" --output_path "$RUN/chains.jsonl" --chain_list "$HEAVY"

# 3. Free the CDRs, fix the framework. Written here rather than with
#    make_fixed_positions_dict.py because that helper applies ONE position list
#    to every entry, and each nanobody's CDRs sit at different positions.
python3 "$REPO/src/evaluate/mpnn_cdr.py" --heavy "$HEAVY" --hlt "$HLT_DIR" mask \
  --parsed "$RUN/parsed.jsonl" --out "$RUN/fixed_positions.jsonl"

# 4. Redesign.
python3 "$MPNN/protein_mpnn_run.py" \
  --jsonl_path "$RUN/parsed.jsonl" \
  --chain_id_jsonl "$RUN/chains.jsonl" \
  --fixed_positions_jsonl "$RUN/fixed_positions.jsonl" \
  --out_folder "$RUN" \
  --num_seq_per_target "$NSEQ" --sampling_temp "$TEMP" --seed "$SEED" --batch_size 1

# 5. Per-CDR recovery vs native.
python3 "$REPO/src/evaluate/mpnn_cdr.py" --heavy "$HEAVY" --hlt "$HLT_DIR" recover \
  --parsed "$RUN/parsed.jsonl" --seqs "$RUN/seqs" \
  --out "$REPO/results/cdr_recovery.csv"

echo
echo "run dir: $RUN"
echo "next:    python3 src/evaluate/fig3_recovery.py"
