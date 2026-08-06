#!/usr/bin/env bash
#SBATCH --job-name=repl_concord
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/replicate_concordance_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/replicate_concordance_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
echo "host=$(hostname) start=$(date -Iseconds)"
"$PY" "$NEW/scripts/replicate_concordance.py" --seeds 5 --min_pooled 50 --min_half 10
echo "done=$(date -Iseconds)"
