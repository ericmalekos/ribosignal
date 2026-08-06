#!/usr/bin/env bash
#SBATCH --job-name=orthrus_tok
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=1-8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orthrus/tok_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orthrus/tok_%A_%a.err
set -euo pipefail
i=${SLURM_ARRAY_TASK_ID}
REPO=/private/groups/carpenterlab/emalekos/RNAZoo_meta
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
WEIGHTS_HOST=/private/groups/carpenterlab/emalekos/orthrus6t_weights
NEW=$REPO/RNAZoo/experiments/riboseq_signal_model
mkdir -p "$NEW/data/orthrus4t_token_emb/chunk_${i}" "$NEW/logs/orthrus"

# Same universe/chunks as RiNALMo (data/rinalmo_chunks/chunk_${i}.fa) and the same
# 10000-nt cap, so the Orthrus per-token set matches the RiNALMo one 1:1 and the
# train/val/test split is identical across embedding backends.
INPUT="$NEW/data/rinalmo_chunks/chunk_${i}.fa"
OUTDIR="$NEW/data/orthrus4t_token_emb/chunk_${i}"

echo "chunk=$i host=$(hostname) start=$(date -Iseconds)"; nvidia-smi -L
singularity exec --nv --no-home --cleanenv \
  --env PYTHONNOUSERSITE=1 --env HF_HUB_OFFLINE=1 --env TRANSFORMERS_OFFLINE=1 \
  --bind ${REPO}:${REPO} --bind ${WEIGHTS_HOST}:/weights \
  "$SIF" /opt/conda/envs/orthrus/bin/python \
    $NEW/scripts/extract_orthrus4t_pertoken.py \
    --input "$INPUT" \
    --output "$OUTDIR" \
    --weights /weights/orthrus-large-4-track \
    --max-length 10000 \
    --checkpoint-every 500 \
    --dtype float16
echo "chunk=$i done=$(date -Iseconds)"
n=$(ls "$OUTDIR"/*_tokens.npy 2>/dev/null | wc -l)
echo "chunk=$i token files: $n"
