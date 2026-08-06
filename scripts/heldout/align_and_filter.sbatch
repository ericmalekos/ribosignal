#!/usr/bin/env bash
#SBATCH --job-name=ho_align
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/align_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/align_%A_%a.err
# Array runner: one task per Ribo-seq SRR of $DATASET. Submit with:
#   sbatch --array=0-<N-1> --export=ALL,DATASET=<dataset> align_and_filter.sbatch
set -euo pipefail
HELDOUT_DIR=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/heldout
export DATASET="${DATASET:?export DATASET=<dataset>}"
source "$HELDOUT_DIR/heldout_config.sh"
srrs=($RIBO_SRRS)
SRR="${srrs[$SLURM_ARRAY_TASK_ID]:?array index out of range}"
echo "task $SLURM_ARRAY_TASK_ID -> $DATASET $SRR on $(hostname)"
bash "$HELDOUT_DIR/align_and_filter.sh" "$DATASET" "$SRR"
