#!/usr/bin/env bash
#SBATCH --job-name=ribo_backend
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/backend_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/backend_%A_%a.err
set -euo pipefail

# RiNALMo vs Orthrus-4track vs concat, fold 0, on the identical gene-disjoint split.
# Training only reads pre-extracted per-token .npy + numpy, so the RiNALMo CUDA SIF (torch)
# serves all three backends; the Orthrus SIF is only needed for extraction.
MODES=(rinalmo orthrus concat)
BACKEND=${MODES[$SLURM_ARRAY_TASK_ID]}

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
OUT=$NEW/results/backend_cmp/${BACKEND}_f0
mkdir -p "$NEW/logs/train" "$OUT"

echo "backend=$BACKEND host=$(hostname) start=$(date -Iseconds) gpu=${CUDA_VISIBLE_DEVICES:-NA}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" \
    --emb_backend "$BACKEND" \
    --test_fold 0 \
    --epochs 40 --patience 6 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 24000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0

echo "backend=$BACKEND done=$(date -Iseconds)"

# Extra metrics (magnitude + posture-B proxy) on the same run, same SIF.
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || true
echo "backend=$BACKEND extra_done=$(date -Iseconds)"
