#!/usr/bin/env bash
#SBATCH --job-name=rna_mm1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=1-55%12
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rna_mm1_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rna_mm1_%A_%a.err
set -euo pipefail
# mm1 re-align of the re-downloaded Chothani RNA-seq: align_rnaseq_sample.sh already uses
# --outFilterMultimapNmax 1 -> per_sample/<SRR>_coverage.hd5 (fresh mm1; none exist so no skip).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
mkdir -p $NEW/logs
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NEW/data/chothani_mm1_srr.txt")
[ -n "$SRR" ] || { echo "no SRR at row ${SLURM_ARRAY_TASK_ID}"; exit 1; }
bash "$NEW/scripts/align_rnaseq_sample.sh" "$SRR"
