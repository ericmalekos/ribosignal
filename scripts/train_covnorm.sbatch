#!/usr/bin/env bash
#SBATCH --job-name=ribo_covnorm
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=0-9%10
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/covnorm_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/covnorm_%A_%a.err
set -euo pipefail

# Depth-normalized RNA-seq input (--cov_norm global_mean: divide coverage by the dataset global
# mean per-nt depth before log1p, removing absolute sequencing depth so the input transfers across
# datasets). Retrain the baseline and the best config (orf_v2_attn) at all 5 folds and compare to
# the existing raw-depth runs (results/improve/*), which are cov_norm=raw. Everything else matches
# those runs (seed 0, same hypers, same split). Output to results/covnorm/ so the raw runs stay.
TASKS=(
  baseline:0 baseline:1 baseline:2 baseline:3 baseline:4
  orf_v2_attn:0 orf_v2_attn:1 orf_v2_attn:2 orf_v2_attn:3 orf_v2_attn:4
)
ENTRY=${TASKS[$SLURM_ARRAY_TASK_ID]}
CFG=${ENTRY%%:*}
FOLD=${ENTRY##*:}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track.npy
case $CFG in
  baseline)     EXTRA="";                                  BUDGET=24000 ;;
  orf_v2_attn)  EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000; ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/covnorm/${CFG}_f${FOLD}
mkdir -p "$NEW/logs/train" "$OUT"

echo "cfg=$CFG fold=$FOLD cov_norm=global_mean budget=$BUDGET extra='$EXTRA' orftrack=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" \
    --emb_backend rinalmo --test_fold "$FOLD" \
    --cov_norm global_mean \
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
