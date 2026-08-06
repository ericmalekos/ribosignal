#!/usr/bin/env bash
#SBATCH --job-name=rs_mamba
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/mamba_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/mamba_%j.err
set -euo pipefail
# Bidirectional-Mamba mixer variant of orf_v2_attn_onehot -- a controlled swap vs the transformer
# (baseline results/loto/orf_v2_attn_onehot_holdout_Hepatocytes and results/ablation/attn4_onehot_...).
# ONLY the mixer changes: 2 transformer layers -> 2 bidirectional Mamba blocks (--mixer mamba). Runs in
# the Orthrus SIF (mamba_ssm 1.2.0.post1, torch 2.2). Smoke via MAXTX/EPOCHS/EVALCAP/RUN_NAME overrides.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
ORFTRACK=$NEW/data/packed/orf_track_v2.npy
MAXTX=${MAXTX:-6000}; EPOCHS=${EPOCHS:-30}; EVALCAP=${EVALCAP:-2500}; RUN_NAME=${RUN_NAME:-mamba_onehot}
OUT=$NEW/results/ablation/${RUN_NAME}_holdout_Hepatocytes
mkdir -p "$OUT" "$NEW/logs/train"
echo "mamba-mixer run=$RUN_NAME maxtx=$MAXTX epochs=$EPOCHS host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-?} start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" --holdout Hepatocytes \
    --emb_backend onehot --input_mode both --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 --max_tx_per_tissue "$MAXTX" --resume \
    --use_orf_track --n_attn_layers 2 --mixer mamba \
    --epochs "$EPOCHS" --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 --eval_cap "$EVALCAP" --num_workers 8 --device cuda --seed 0

echo "-- eval_extra --"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" || echo "(eval_extra failed, non-fatal)"
echo "mamba-mixer $RUN_NAME done=$(date -Iseconds)"
