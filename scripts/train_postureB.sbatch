#!/usr/bin/env bash
#SBATCH --job-name=postureB
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=16:00:00
#SBATCH --array=0-3
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/postureB_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/postureB_%A_%a.err
set -euo pipefail
# Posture-B multimap sensitivity check (results.md Task 21). 2x2: {baseline, orf_v2_attn} x
# {A = full posture-A universe, B = highest-expressed-isoform-per-gene posture-B universe}. Onehot,
# fold-0 held-out-chrom, identical code/seed -- the ONLY variable is the transcript universe (the
# posture-B include-set via RIBO_REPRESENTATIVE_TX). Both orf configs use the no-Kozak track (the new
# default, results.md Task 20) so A-vs-B isolates the multimap posture, not the Kozak choice. Retraining
# BOTH postures fresh (rather than reusing an old posture-A run) keeps the comparison airtight.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
REP=$NEW/data/fibroblast_representative_tx.txt
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
TASKS=(baseline_A baseline_B orf_v2_attn_A orf_v2_attn_B)
CFG=${TASKS[$SLURM_ARRAY_TASK_ID]}
case $CFG in
  baseline_A)    EXTRA="";                                BUDGET=24000; REPENV="" ;;
  baseline_B)    EXTRA="";                                BUDGET=24000; REPENV="$REP" ;;
  orf_v2_attn_A) EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000; REPENV="" ;;
  orf_v2_attn_B) EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000; REPENV="$REP" ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac
OUT=$NEW/results/postureB/${CFG}_onehot_f0
mkdir -p "$NEW/logs/train" "$OUT"
# ORF track env only matters when --use_orf_track is passed; harmless otherwise.
ORFENV=(--env RIBO_ORF_TRACK="$ORFTRACK")
REPENVARG=(); [ -n "$REPENV" ] && REPENVARG=(--env RIBO_REPRESENTATIVE_TX="$REPENV")

echo "cfg=$CFG posture=$([ -n "$REPENV" ] && echo B || echo A) budget=$BUDGET extra='$EXTRA' host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  "${ORFENV[@]}" "${REPENVARG[@]}" -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train.py" \
    --out "$OUT" \
    --emb_backend onehot --test_fold 0 \
    $EXTRA \
    --epochs 40 --patience 6 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget $BUDGET --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "cfg=$CFG train_done=$(date -Iseconds)"

# eval_extra on the run's OWN test set (posture-A -> full fold-0 test; posture-B -> representative
# fold-0 test), so the pooled medians are directly comparable per the roadmap.
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  "${ORFENV[@]}" "${REPENVARG[@]}" -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || true
echo "cfg=$CFG extra_done=$(date -Iseconds)"
