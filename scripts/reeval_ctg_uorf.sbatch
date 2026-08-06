#!/usr/bin/env bash
#SBATCH --job-name=ribo_reeval_ctg
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=0-9%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/reeval_ctg_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/reeval_ctg_%A_%a.err
set -euo pipefail

# Re-eval the 10 CTG cells (orf_ctg + orf_ctg_attn x 5 folds) with the updated eval_extra so they
# get the uORF 5'UTR / dORF 3'UTR blocks + utr columns in pertx.tsv, matching the other configs.
# Their in-job eval ran with the pre-uORF eval_extra. v3 track (orf_track_v3.npy) so the reloaded
# model sees the ATG+CTG features it was trained on. Gate this afterok on the CTG training array.
TASKS=(
  orf_ctg:0 orf_ctg:1 orf_ctg:2 orf_ctg:3 orf_ctg:4
  orf_ctg_attn:0 orf_ctg_attn:1 orf_ctg_attn:2 orf_ctg_attn:3 orf_ctg_attn:4
)
ENTRY=${TASKS[$SLURM_ARRAY_TASK_ID]}
CFG=${ENTRY%%:*}
FOLD=${ENTRY##*:}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v3.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/improve/${CFG}_f${FOLD}

if [ ! -f "$OUT/best.pt" ]; then echo "SKIP $CFG f$FOLD: no best.pt"; exit 0; fi
echo "reeval-ctg cfg=$CFG fold=$FOLD orftrack=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || echo "WARN reeval failed $CFG f$FOLD"
echo "reeval-ctg cfg=$CFG fold=$FOLD done=$(date -Iseconds)"
