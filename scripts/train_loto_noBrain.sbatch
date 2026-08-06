#!/usr/bin/env bash
#SBATCH --job-name=loto_noBrain
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_noBrain_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_noBrain_%j.err
set -euo pipefail
# Retrain the deployed one-hot orf_v2_attn WITHOUT Brain (noise; period_obs 0.044), holding out Hepatocytes
# (same eval fold as the deployed model, for a clean comparison). Trains on the existing mm20 coverage packs
# (COV_TAG default) or the mm1 packs (COV_TAG=_mm1) if/when the mm1 re-align lands. 7 tissues:
# Fibroblast, VSMC, ES, Fat, HA_EC, HCAEC, HUVEC.
# KOZAK: NEVER set the Kozak factor. The Task-20 ablation showed the hand-picked heuristic is redundant with
# what the model learns from raw sequence (mildly harmful cross-tissue), so all retrains use the NO-KOZAK ORF
# track (orf_track_v2_nokozak.npy == build_orf_track.py --mode ext --kozak none), row-aligned to the same
# Fibroblast universe as the heuristic v2 track. (The deployed model stays heuristic-pinned; that is historical
# only -- do not build new heuristic tracks.)
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
TAG=${TAG:-mm20}    # provenance tag in the out dir name; set TAG=mm1 for the mm1-coverage retrain
OUT=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_${TAG}_holdout_Hepatocytes
mkdir -p "$NEW/logs/train" "$OUT"
echo "TAG=$TAG holdout=Hepatocytes train_tissues=7(noBrain) host=$(hostname) start=$(date -Iseconds)"
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
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "train_done=$(date -Iseconds) -> $OUT"
