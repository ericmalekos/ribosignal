#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_bamcov
#SBATCH --partition=medium
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=1-32%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaseq_star/bamcov_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaseq_star/bamcov_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
LIST=$NEW/data/fibroblast_rnaseq_srr.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -n "$SRR" ] || { echo "no SRR at row ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
bash "$NEW/scripts/align_rnaseq_sample_bamcov.sh" "$SRR"
