#!/usr/bin/env bash
#SBATCH --job-name=ribo_train_f0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/fold0_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/fold0_%j.err
set -euo pipefail

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
OUT=$NEW/results/fold0_rinalmo
mkdir -p "$NEW/logs/train" "$OUT"

echo "host=$(hostname) start=$(date -Iseconds) gpu=$CUDA_VISIBLE_DEVICES"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# same-path bind so the scripts' hardcoded /private/... paths resolve read-write inside.
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" \
    --epochs 40 --patience 6 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 24000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0

echo "done=$(date -Iseconds)"
