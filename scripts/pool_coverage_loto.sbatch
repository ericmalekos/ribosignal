#!/usr/bin/env bash
#SBATCH --job-name=loto_pool_cov
#SBATCH --partition=medium
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --array=0-7%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pool_cov_loto_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pool_cov_loto_%A_%a.err
set -euo pipefail
# Pool each non-Fibroblast tissue's per-sample RNAseq coverage into one tissue file.
# Runs after align_rnaseq_loto (afterok): the per-sample coverage hd5 must all exist.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
TISSUES=(Brain ES Fat HA_EC HCAEC HUVEC Hepatocytes VSMC)
T=${TISSUES[$SLURM_ARRAY_TASK_ID]}
echo "tissue=$T host=$(hostname) start=$(date -Iseconds)"
$PY "$NEW/scripts/pool_rnaseq_coverage_tissue.py" "$T"
echo "tissue=$T done=$(date -Iseconds)"
