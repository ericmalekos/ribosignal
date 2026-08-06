#!/usr/bin/env bash
#SBATCH --job-name=ribo_improve_v2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/improve_v2_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/improve_v2_%A_%a.err
set -euo pipefail

# Non-AUG-aware ORF track (orf_track_v2.npy, channel 3 = graded start_ext over AUG + 9 near-
# cognates). Same fold-0 split/backend/hyperparameters as the ATG-only v1 improve runs; only
# the ORF track file changes (via RIBO_ORF_TRACK), so v1(atg) vs v2(ext) isolates the value of
# non-AUG start information. attn-only is shared with the v1 array (track-independent).
MODES=(orf_v2 orf_v2_attn)
CFG=${MODES[$SLURM_ARRAY_TASK_ID]}
case $CFG in
  orf_v2)      EXTRA="--use_orf_track";                   BUDGET=24000 ;;
  orf_v2_attn) EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000 ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2.npy
OUT=$NEW/results/improve/${CFG}_f0
mkdir -p "$NEW/logs/train" "$OUT"

echo "cfg=$CFG budget=$BUDGET extra='$EXTRA' orf=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" \
    --emb_backend rinalmo --test_fold 0 \
    $EXTRA \
    --epochs 40 --patience 6 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget $BUDGET --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0

echo "cfg=$CFG done=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || true
echo "cfg=$CFG extra_done=$(date -Iseconds)"
