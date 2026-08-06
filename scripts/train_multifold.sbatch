#!/usr/bin/env bash
#SBATCH --job-name=ribo_mfold
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=0-24%12
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/mfold_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/mfold_%A_%a.err
set -euo pipefail

# Multi-fold confirmation of the fold-0 improvement matrix (results.md Task 8). Extends every
# config to all 5 gene-disjoint chromosome folds, so the lncRNA deltas rest on ~1,239 pooled
# transcripts instead of the 264 in fold 0. Fold-0 cells for the 5 non-baseline configs are
# REUSED from the existing results/improve/*_f0 runs; this array only trains the missing cells:
#   baseline: folds 0-4 (fresh, improve lineage; baseline_f0 double-checks backend_cmp 0.612 pc)
#   orf, attn, orf_attn, orf_v2, orf_v2_attn: folds 1-4 (fold 0 reused)
# Settings match backend_cmp/rinalmo_f0 (the fold-0 matrix baseline) exactly: emb_backend
# rinalmo, budget 24000 (16000 for attention, O(L^2) memory), seed 0, identical CNN hypers.
# orf_v2* point RIBO_ORF_TRACK at the non-AUG-aware track; RIBO_ORF_TRACK is a no-op for configs
# that do not pass --use_orf_track (dataset.py only reads it when the ORF track is enabled).
TASKS=(
  baseline:0 baseline:1 baseline:2 baseline:3 baseline:4
  orf:1 orf:2 orf:3 orf:4
  attn:1 attn:2 attn:3 attn:4
  orf_attn:1 orf_attn:2 orf_attn:3 orf_attn:4
  orf_v2:1 orf_v2:2 orf_v2:3 orf_v2:4
  orf_v2_attn:1 orf_v2_attn:2 orf_v2_attn:3 orf_v2_attn:4
)
ENTRY=${TASKS[$SLURM_ARRAY_TASK_ID]}
CFG=${ENTRY%%:*}
FOLD=${ENTRY##*:}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track.npy
case $CFG in
  baseline)    EXTRA="";                                  BUDGET=24000 ;;
  orf)         EXTRA="--use_orf_track";                   BUDGET=24000 ;;
  attn)        EXTRA="--n_attn_layers 2";                 BUDGET=16000 ;;
  orf_attn)    EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000 ;;
  orf_v2)      EXTRA="--use_orf_track";                   BUDGET=24000; ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;;
  orf_v2_attn) EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000; ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/improve/${CFG}_f${FOLD}
mkdir -p "$NEW/logs/train" "$OUT"

echo "cfg=$CFG fold=$FOLD budget=$BUDGET extra='$EXTRA' orftrack=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" \
    --emb_backend rinalmo --test_fold "$FOLD" \
    $EXTRA \
    --epochs 40 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget "$BUDGET" --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0

echo "cfg=$CFG fold=$FOLD train_done=$(date -Iseconds)"

# count-head magnitude + posture-B + per-tx pertx.tsv (eval_extra caps its budget at cfg budget)
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || echo "WARN eval_extra failed $CFG f$FOLD"
echo "cfg=$CFG fold=$FOLD extra_done=$(date -Iseconds)"
