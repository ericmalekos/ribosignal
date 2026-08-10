#!/usr/bin/env bash
#SBATCH --job-name=jan_pois002
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/jan_pois002_%j.out
set -euo pipefail
# The Janich-driven model arm's Poisson sweep has theta 0.1/0.2/0.35 but not 0.02, which is the
# CDS-anchored operating point used for the Wang arm. Using a different theta per arm would make the
# two model rows non-comparable in the three-way table, so generate the missing point.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
OUT=$NEW/results/o2_liver/janich_wanguni/dropin_sweep_poisson/theta_0.02; mkdir -p "$OUT"
"$RPY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$NEW/results/o2_liver/janich_wanguni/pred_profiles.npz" \
  --annot "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot" --variant pred_preddepth --out "$OUT" \
  --min_aa 5 --pval 0.05 --pred_poisson --pred_scale 0.02
echo "done=$(date -Iseconds): $(($(wc -l < "$OUT"/pred_preddepth_collapsed.txt) - 1)) calls"
