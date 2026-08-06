#!/usr/bin/env bash
#SBATCH --job-name=ribo_gpu_smoke
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/gpu_smoke_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/gpu_smoke_%j.err
set -euo pipefail

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
mkdir -p "$NEW/logs/train"
echo "host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# exact real-run invocation, but pointed at the 1500-tx dry pack (RIBO_PACK_DIR) for a
# fast GPU end-to-end validation of the SIF + CUDA + dataloader-in-SIF path.
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_PACK_DIR="$NEW/data/packed_dryrun" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$NEW/results/gpu_smoke" \
    --epochs 2 --patience 5 \
    --channels 128 --n_blocks 6 --dropout 0.1 \
    --budget 16000 --eval_cap 150 --num_workers 4 \
    --max_train 400 --device cuda --seed 0

echo "done=$(date -Iseconds)"
