#!/usr/bin/env bash
#SBATCH --job-name=ho_rnacov
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/rnacov_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/rnacov_%A_%a.err
# Array runner: one task per RNAseq SRR of $DATASET -> per-nt transcriptome coverage hd5. Submit with:
#   sbatch --array=0-<N-1> --export=ALL,DATASET=<dataset> rnaseq_coverage_heldout.sbatch
set -euo pipefail
HELDOUT_DIR=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/heldout
export DATASET="${DATASET:?export DATASET=<dataset>}"
# shellcheck source=/dev/null
source "$HELDOUT_DIR/heldout_config.sh"
srrs=($RNA_SRRS)
SRR="${srrs[$SLURM_ARRAY_TASK_ID]:?array index out of range}"
echo "task $SLURM_ARRAY_TASK_ID -> $DATASET RNA $SRR on $(hostname)"
bash "$HELDOUT_DIR/rnaseq_coverage_heldout.sh" "$DATASET" "$SRR"
