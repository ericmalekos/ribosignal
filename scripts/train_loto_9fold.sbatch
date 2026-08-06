#!/usr/bin/env bash
#SBATCH --job-name=loto_9fold
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --array=0-8%3
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto9_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto9_%A_%a.err
set -euo pipefail
# FULL leave-one-tissue-out: each of the 9 Chothani tissues held out once, trained on the other 8,
# using orf_v2_attn (the pilot winner). Gives the per-tissue cross-tissue spread. Throttled %3 so it
# does not monopolize the shared GPU partition (9 folds x ~12h -> ~3 waves). Backend selectable:
#   sbatch scripts/train_loto_9fold.sbatch                       # default rinalmo
#   sbatch --export=ALL,BACKEND=orthrus scripts/train_loto_9fold.sbatch
# STAGED -- launch only after the pilot 4-way table validates and the backend is chosen.
BACKEND=${BACKEND:-rinalmo}
TISSUES=(Fibroblast VSMC ES Fat HA_EC Brain HCAEC HUVEC Hepatocytes)
HOLDOUT=${TISSUES[$SLURM_ARRAY_TASK_ID]}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2.npy   # orf_v2_attn uses the non-AUG track
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/loto_9fold/${BACKEND}_orf_v2_attn_holdout_${HOLDOUT}
mkdir -p "$NEW/logs/train" "$OUT"

echo "backend=$BACKEND holdout=$HOLDOUT host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout "$HOLDOUT" \
    --emb_backend "$BACKEND" --cov_norm global_mean \
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

echo "backend=$BACKEND holdout=$HOLDOUT train_done=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 --device cuda || echo "WARN eval_extra failed $HOLDOUT"
echo "backend=$BACKEND holdout=$HOLDOUT extra_done=$(date -Iseconds)"
