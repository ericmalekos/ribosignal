#!/usr/bin/env bash
#SBATCH --job-name=loto_peaky
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_peaky_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_peaky_%j.err
set -euo pipefail
# O3 track 2: the SAME model / config / holdout as train_loto_noBrain.sbatch (nokozak mm1), PLUS the
# anti-smoothing entropy-gap penalty (--peakiness_weight). This trains a variant IN PARALLEL with the
# already-trained current model (results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes)
# so the two are directly comparable. Hypothesis: pushing the model to predict structured/peaked profiles
# (instead of the smooth hedge) makes the O3 shape features discriminative and reduces RiboCode's novel-ORF
# over-calling. KOZAK is never set (no-Kozak ORF track), same as every retrain.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
TAG=${TAG:-mm1}         # coverage provenance; mm1 = unique-mapper packs (the final-model coverage)
PK=${PK:-0.2}           # entropy-gap weight; PKTAG replaces the dot so out dirs never collide across weights
PKTAG=${PK/./p}
OUT=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_${TAG}_peaky${PKTAG}_holdout_Hepatocytes
mkdir -p "$NEW/logs/train" "$OUT"
echo "TAG=$TAG peakiness_weight=$PK holdout=Hepatocytes train_tissues=7(noBrain) host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes \
    --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC \
    --emb_backend onehot --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    --use_orf_track --n_attn_layers 2 \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 --peakiness_weight "$PK" \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "train_done=$(date -Iseconds) -> $OUT"
