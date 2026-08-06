#!/usr/bin/env bash
#SBATCH --job-name=ribo_ctg
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=0-9%10
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/ctg_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/ctg_%A_%a.err
set -euo pipefail

# ATG+CTG ORF track (orf_track_v3.npy): occupancy over ATG-or-CTG-initiated ORFs + 2-level start
# channel (ATG=1.0, CTG=0.5). Tests whether adding the single dominant alternative start (CTG) to
# the ORF definition helps, vs ATG-only (orf, +0.056 lncRNA) and the full 10-codon near-cognate
# set (orf_v2, +0.033). Run at all 5 gene-disjoint folds so it is directly comparable to those on
# the pooled ~1,239-lncRNA basis. Same improve-lineage settings (budget 24000 / 16000 attn, seed
# 0, identical CNN hypers). Two configs: orf_ctg (track only) and orf_ctg_attn (+ 2 attn layers).
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
case $CFG in
  orf_ctg)      EXTRA="--use_orf_track";                   BUDGET=24000 ;;
  orf_ctg_attn) EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000 ;;
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

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || echo "WARN eval_extra failed $CFG f$FOLD"
echo "cfg=$CFG fold=$FOLD extra_done=$(date -Iseconds)"
