#!/usr/bin/env bash
#SBATCH --job-name=fib_psite_target
#SBATCH --partition=medium
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/build_target_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/build_target_%j.err
set -euo pipefail
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
echo "host=$(hostname) start=$(date -Iseconds)"
$PY scripts/build_fibroblast_psite_target.py
echo "done=$(date -Iseconds)"
