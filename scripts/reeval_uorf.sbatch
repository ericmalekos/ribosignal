#!/usr/bin/env bash
#SBATCH --job-name=ribo_reeval
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=0-29%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/reeval_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/reeval_%A_%a.err
set -euo pipefail

# Re-run the UPDATED eval_extra (now also emits the uORF 5'UTR + dORF 3'UTR profile-Pearson
# blocks and the utr columns in pertx.tsv) on the 30 completed improvement-matrix cells.
# Eval-only: loads best.pt, no retraining. eval_cap 0 = all test tx, so pertx.tsv is refreshed
# in full (same full-transcript profile_pearson as before + the new UTR columns). RIBO_ORF_TRACK
# must match each run's training track or the reloaded model sees wrong ORF features.
TASKS=(
  baseline:0 baseline:1 baseline:2 baseline:3 baseline:4
  orf:0 orf:1 orf:2 orf:3 orf:4
  orf_v2:0 orf_v2:1 orf_v2:2 orf_v2:3 orf_v2:4
  attn:0 attn:1 attn:2 attn:3 attn:4
  orf_attn:0 orf_attn:1 orf_attn:2 orf_attn:3 orf_attn:4
  orf_v2_attn:0 orf_v2_attn:1 orf_v2_attn:2 orf_v2_attn:3 orf_v2_attn:4
)
ENTRY=${TASKS[$SLURM_ARRAY_TASK_ID]}
CFG=${ENTRY%%:*}
FOLD=${ENTRY##*:}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track.npy
case $CFG in orf_v2|orf_v2_attn) ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;; esac
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/improve/${CFG}_f${FOLD}

if [ ! -f "$OUT/best.pt" ]; then echo "SKIP $CFG f$FOLD: no best.pt"; exit 0; fi
echo "reeval cfg=$CFG fold=$FOLD orftrack=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || echo "WARN reeval failed $CFG f$FOLD"
echo "reeval cfg=$CFG fold=$FOLD done=$(date -Iseconds)"
