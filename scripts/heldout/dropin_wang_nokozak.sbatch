#!/usr/bin/env bash
#SBATCH --job-name=dropin_wang_nk
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/dropin_wang_nk_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/dropin_wang_nk_%A_%a.err
set -euo pipefail
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
VARIANTS=(real pred_obsdepth pred_preddepth); V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
OUT=$NEW/results/o2_liver/wang_nokozak/dropin; mkdir -p "$OUT"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$NEW/results/o2_liver/wang_nokozak/pred_profiles.npz" \
  --annot "$ECH/data/ribocode_annot_mouse" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
echo "wang dropin $V -> $OUT/${V}_collapsed.txt"
