#!/usr/bin/env bash
#SBATCH --job-name=union_ablate
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/union_ablate_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/union_ablate_%A_%a.err
set -euo pipefail
# INPUT ABLATION against the DEPLOYED recipe (FIGURES_PLAN "what else would strengthen this", item 1).
#
# THE CLAIM UNDER TEST. The whole proteogenomics argument rests on the model being CELL-TYPE-SPECIFIC:
# it is handed a query cell type's own RNA-seq and is supposed to predict that cell type's
# translation. If the model could reach the same accuracy with the coverage channel zeroed, then it
# is a sequence-only ORF prior, its output is the same for every cell type, and every "cell-type-
# specific discovery" claim in the manuscript collapses. This job measures that directly.
#
#   arm 0  input_mode=emb   sequence + ORF track, COVERAGE ZEROED   -> "no RNA-seq"
#   arm 1  input_mode=cov   coverage only, SEQUENCE + ORF ZEROED    -> the symmetric counterpart
#
# The `both` arm is not retrained: it already exists as the deployed run
#   results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
# and re-running it would only add seed noise to the contrast.
#
# WHY NOT scripts/train_loto_ablation.sbatch, which already takes INPUT_MODE. That script is pinned to
# the OLD Fibroblast-universe config: data/packed/orf_track_v2.npy, no --kozak none, no RIBO_PACK_SUFFIX.
# Ablating there would contrast against a model that is not the one shipped, so the result would not
# speak to the deployed claim. This script clones train_loto_union.sbatch line for line and changes
# ONLY --input_mode.
#
#   sbatch --array=0-1 scripts/train_union_inputablation.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
# The deployed model is one-hot + transformer on the UNION universe. The union track was BUILT with
# --kozak none, so data/packed_union/orf_track_v2_nokozak.npy is a SYMLINK to orf_track_v2.npy --
# same bytes, two names. train_loto_union.sbatch passes the plain name, so this passes the plain name
# too: the ablation must differ from the deployed run in --input_mode and nothing else, including the
# string recorded in the run config.
ORFTRACK=$NEW/data/packed_union/orf_track_v2.npy
[ -s "$ORFTRACK" ] || { echo "missing $ORFTRACK" >&2; exit 1; }

MODES=(emb cov)
MODE=${MODES[$SLURM_ARRAY_TASK_ID]}
OUT=$NEW/results/ablation/orf_v2_attn_onehot_union_nokozak_input_${MODE}_holdout_Hepatocytes
mkdir -p "$OUT" "$NEW/logs/train"
echo "union input-ablation mode=$MODE orftrack=$(basename "$ORFTRACK") host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-?} start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_PACK_SUFFIX=union \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes \
    --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC \
    --emb_backend onehot --input_mode "$MODE" --cov_norm global_mean \
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

echo "-- eval_extra --"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_PACK_SUFFIX=union --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" || echo "(eval_extra failed, non-fatal)"
echo "union input-ablation $MODE done=$(date -Iseconds)"
