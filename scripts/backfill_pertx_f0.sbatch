#!/usr/bin/env bash
#SBATCH --job-name=ribo_pertx_f0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/pertx_f0_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/pertx_f0_%j.err
set -euo pipefail

# Backfill pertx.tsv (per-tx profile-Pearson dump, newly added to eval_extra) for the 5 REUSED
# fold-0 improve runs, so the multifold pooled aggregation can include their fold-0 transcripts.
# Eval-only (loads best.pt); no retraining. baseline fold-0 comes from the fresh baseline_f0 run
# in the multifold array, not here. orf_v2* must point RIBO_ORF_TRACK at orf_track_v2.npy so the
# reloaded model sees the same ORF features it was trained on.
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model

for r in orf_f0 attn_f0 orf_attn_f0 orf_v2_f0 orf_v2_attn_f0; do
  ORFTRACK=$NEW/data/packed/orf_track.npy
  case $r in orf_v2_f0|orf_v2_attn_f0) ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;; esac
  echo "=== pertx $r (orf=$ORFTRACK) $(date -Iseconds) ==="
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env RIBO_ORF_TRACK="$ORFTRACK" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/eval_extra.py" --run "$NEW/results/improve/$r" --eval_cap 0 \
    || echo "WARN pertx failed for $r"
done
echo "=== pertx_f0 done $(date -Iseconds) ==="
