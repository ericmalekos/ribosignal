#!/usr/bin/env bash
#SBATCH --job-name=rs_loto_ablate
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_ablate_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_ablate_%j.err
set -euo pipefail
# One LOTO ablation run (Hepatocytes holdout), replicating the exact config of
# results/loto/orf_v2_attn_onehot_holdout_Hepatocytes and changing ONE variable via --export:
#   RUN_NAME, BACKEND, INPUT_MODE (both|emb|cov), N_ATTN, CHANNELS.
# Everything else is pinned to that orf_v2_attn config so each run is a clean single-variable contrast.
# Trains then self-evaluates with eval_extra.py (pc/lncRNA/uORF metrics). Distinct from the within-tissue
# train_ablation.sbatch (that one is train.py fold0; this is train_loto.py cross-tissue).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
ORFTRACK=$NEW/data/packed/orf_track_v2.npy
RUN_NAME=${RUN_NAME:?set RUN_NAME}
BACKEND=${BACKEND:-onehot}
INPUT_MODE=${INPUT_MODE:-both}
N_ATTN=${N_ATTN:-2}
CHANNELS=${CHANNELS:-256}
OUT=$NEW/results/ablation/${RUN_NAME}_holdout_Hepatocytes
mkdir -p "$OUT" "$NEW/logs/train"
echo "loto-ablation run=$RUN_NAME backend=$BACKEND input_mode=$INPUT_MODE n_attn=$N_ATTN channels=$CHANNELS host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-?} start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" --holdout Hepatocytes \
    --emb_backend "$BACKEND" --input_mode "$INPUT_MODE" --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 --max_tx_per_tissue 6000 --resume \
    --use_orf_track --n_attn_layers "$N_ATTN" --n_heads 8 \
    --epochs 30 --patience 6 --warmup 2 \
    --channels "$CHANNELS" --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 --device cuda --seed 0

echo "-- eval_extra --"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" || echo "(eval_extra failed, non-fatal)"
echo "loto-ablation $RUN_NAME done=$(date -Iseconds)"
