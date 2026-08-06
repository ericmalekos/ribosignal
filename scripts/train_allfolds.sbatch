#!/usr/bin/env bash
#SBATCH --job-name=ribo_train_cv
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/cv_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/cv_%A_%a.err
set -euo pipefail

# 5-fold cross-validation: array task k holds out fold k as test. Launch AFTER the fold-0
# baseline (train_fold0.sbatch) confirms the approach.
FOLD=${SLURM_ARRAY_TASK_ID}
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
OUT=$NEW/results/fold${FOLD}_rinalmo
mkdir -p "$NEW/logs/train" "$OUT"

echo "host=$(hostname) fold=$FOLD start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" --test_fold "$FOLD" \
    --epochs 40 --patience 6 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 24000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "fold=$FOLD done=$(date -Iseconds)"
