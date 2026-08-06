#!/usr/bin/env bash
#SBATCH --job-name=hydra_emb
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=0-11
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hydrarna/emb_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hydrarna/emb_%A_%a.err
set -euo pipefail
# HydraRNA per-token embeddings for the whole Fibroblast universe (36,668 tx), 12 GPU chunks -> per-tx .npy
# in data/hydrarna_token_emb/. Resume-safe (skips existing). Build the tx_index after with build_embedding_index.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-hydrarna-latest.img
FASTA=$NEW/data/fibroblast_universe.fa
OUT=$NEW/data/hydrarna_token_emb
mkdir -p "$OUT" "$NEW/logs/hydrarna"
N=$(grep -c ">" "$FASTA")
NCHUNK=12
PER=$(( (N + NCHUNK - 1) / NCHUNK ))
START=$(( SLURM_ARRAY_TASK_ID * PER ))
END=$(( START + PER ))
# Triton (Mamba kernel compiler) needs a writable $HOME/.triton -> node-local scratch (per array task).
SCRATCH=/data/tmp/emalekos/hydra_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCRATCH/.triton"
echo "chunk $SLURM_ARRAY_TASK_ID: tx [$START,$END) of $N  host=$(hostname) start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env HOME="$SCRATCH" --env TRITON_CACHE_DIR="$SCRATCH/.triton" \
  -B "$BASE" -B /data/tmp/emalekos "$SIF" \
  python3 "$NEW/scripts/extract_hydrarna_pertoken.py" \
    -i "$FASTA" -o "$OUT" --start "$START" --end "$END" --log-every 200
echo "chunk $SLURM_ARRAY_TASK_ID done=$(date -Iseconds)"
