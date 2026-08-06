#!/usr/bin/env bash
#SBATCH --job-name=rinalmo_tok
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=1-8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rinalmo/tok_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rinalmo/tok_%A_%a.err
set -euo pipefail
i=${SLURM_ARRAY_TASK_ID}
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
mkdir -p "$NEW/data/rinalmo_token_emb/chunk_${i}" "$NEW/logs/rinalmo"

# /private is read-only inside the SIF; pass paths under the /work bind mount.
WORKROOT=/work/RNAZoo/experiments/riboseq_signal_model
INPUT="$WORKROOT/data/rinalmo_chunks/chunk_${i}.fa"
OUTDIR="$WORKROOT/data/rinalmo_token_emb/chunk_${i}"

echo "chunk=$i host=$(hostname) start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta:/work \
  "$SIF" \
  python3 /work/Embeddings/rinalmo_lncrna/scripts/extract_lncrna_rinalmo.py \
    --input "$INPUT" \
    --output "$OUTDIR" \
    --max-length 10000 \
    --checkpoint-every 200 \
    --rebuild-pooled \
    --dtype float16
echo "chunk=$i done=$(date -Iseconds)"
n=$(ls "$NEW/data/rinalmo_token_emb/chunk_${i}"/*_tokens.npy 2>/dev/null | wc -l)
echo "chunk=$i token files: $n"
