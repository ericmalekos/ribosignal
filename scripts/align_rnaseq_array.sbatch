#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_cov
#SBATCH --partition=medium
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --array=1-32%12
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaseq_star/cov_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaseq_star/cov_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
LIST=$NEW/data/fibroblast_rnaseq_srr.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -n "$SRR" ] || { echo "no SRR at row ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
bash "$NEW/scripts/align_rnaseq_sample.sh" "$SRR"
