#!/usr/bin/env bash
#SBATCH --job-name=ho_ribocode
#SBATCH --partition=medium
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/ribocode_%A.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/ribocode_%A.err
# Pooled RiboCode P-site build for $DATASET. Submit AFTER its align array, e.g.:
#   sbatch --dependency=afterok:<align_jobid> --export=ALL,DATASET=<dataset> ribocode_heldout.sbatch
set -euo pipefail
HELDOUT_DIR=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/heldout
export DATASET="${DATASET:?export DATASET=<dataset>}"
bash "$HELDOUT_DIR/ribocode_heldout.sh" "$DATASET"
