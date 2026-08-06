#!/usr/bin/env bash
#SBATCH --job-name=ribo_improve
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/improve_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/improve_%A_%a.err
set -euo pipefail

# Improvement A/B vs the RiNALMo fold-0 baseline (results/backend_cmp/rinalmo_f0 = 0.604):
#   orf      = + sequence-derived ORF-candidate track (annotation-free; uORF/dORF/lncRNA safe)
#   attn     = + 2 self-attention layers after the conv body (full-transcript context)
#   orf_attn = both
# Same fold-0 split, backend, seed, and CNN hyperparameters as the baseline; only the input
# track and/or the attention stack change. Evaluate pc AND lncRNA (not pc-only).
MODES=(orf attn orf_attn)
CFG=${MODES[$SLURM_ARRAY_TASK_ID]}
case $CFG in
  orf)      EXTRA="--use_orf_track";                    BUDGET=24000 ;;
  attn)     EXTRA="--n_attn_layers 2";                  BUDGET=16000 ;;
  orf_attn) EXTRA="--use_orf_track --n_attn_layers 2";  BUDGET=16000 ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
OUT=$NEW/results/improve/${CFG}_f0
mkdir -p "$NEW/logs/train" "$OUT"

echo "cfg=$CFG budget=$BUDGET extra='$EXTRA' host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
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
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || true
echo "cfg=$CFG extra_done=$(date -Iseconds)"
