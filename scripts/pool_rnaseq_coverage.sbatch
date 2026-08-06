#!/usr/bin/env bash
#SBATCH --job-name=pool_rnaseq_cov
#SBATCH --partition=medium
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pool_rnaseq_cov_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pool_rnaseq_cov_%j.err
set -euo pipefail
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
echo "host=$(hostname) start=$(date -Iseconds)"
$PY scripts/pool_rnaseq_coverage.py
echo "done=$(date -Iseconds)"
