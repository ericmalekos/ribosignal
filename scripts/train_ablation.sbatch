#!/usr/bin/env bash
#SBATCH --job-name=ribo_ablation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/ablate_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/ablate_%A_%a.err
set -euo pipefail

# fold-0 input ablation: task 0 = embedding-only, task 1 = coverage-only. Compared against
# the both-inputs fold-0 baseline (results/fold0_rinalmo) to quantify each input's value.
MODES=(emb cov)
MODE=${MODES[$SLURM_ARRAY_TASK_ID]}
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
OUT=$NEW/results/fold0_rinalmo_${MODE}only
mkdir -p "$NEW/logs/train" "$OUT"

echo "host=$(hostname) mode=$MODE start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" --test_fold 0 --input_mode "$MODE" \
    --epochs 40 --patience 6 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 24000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "mode=$MODE done=$(date -Iseconds)"
