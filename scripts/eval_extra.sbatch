#!/usr/bin/env bash
#SBATCH --job-name=ribo_eval_extra
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/eval_extra_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/eval_extra_%j.err
set -euo pipefail

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model

RUNS="fold0_rinalmo fold1_rinalmo fold2_rinalmo fold3_rinalmo fold4_rinalmo \
      fold0_rinalmo_embonly fold0_rinalmo_covonly"

echo "host=$(hostname) start=$(date -Iseconds)"
for R in $RUNS; do
  if [ ! -f "$NEW/results/$R/best.pt" ]; then
    echo "SKIP $R (no best.pt yet)"; continue
  fi
  echo "=== eval_extra $R ==="
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/eval_extra.py" --run "$NEW/results/$R" --device cuda
done
echo "done=$(date -Iseconds)"
