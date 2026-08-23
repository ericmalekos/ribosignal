#!/usr/bin/env bash
#SBATCH --job-name=canon_score
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/canon_score_%j.out
set -euo pipefail
# Re-score the 3x3 and re-run codon occupancy on the CANONICAL substrate.
# The off-recipe equivalents are archived at results/_archive_offrecipe_2026_08_13/.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
echo "host=$(hostname) start=$(date -Iseconds)"

echo "=== 3x3 scoring (canonical) ==="
"$PY" "$NEW/scripts/score_liver3x3.py" \
  --results "$NEW/results/mouse_liver_3x3_canon" \
  --out "$NEW/results/mouse_liver_3x3_canon/scored"

echo "=== codon occupancy (canonical) ==="
bash "$NEW/scripts/run_codon_occupancy_all.sh"
"$PY" "$NEW/scripts/codon_compare.py" --dir "$NEW/results/codon_occupancy_canon" \
  --out "$NEW/results/codon_occupancy_canon/_compare_asite.json" --site asite
echo "done=$(date -Iseconds)"
