#!/usr/bin/env bash
#SBATCH --job-name=ribo_postpool
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/postpool_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/postpool_%j.err
set -euo pipefail
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
echo "host=$(hostname) start=$(date -Iseconds)"
# 1) validate the pooled coverage (axis, controls, salmon + P-site correlations)
"$PY" "$NEW/scripts/validate_coverage.py"
# 2) pack target + coverage for the universe -> data/packed/ (no args = pooled + full)
"$PY" "$NEW/scripts/pack_target_coverage.py"
echo "done=$(date -Iseconds)"
