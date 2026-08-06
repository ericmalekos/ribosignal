#!/usr/bin/env bash
#SBATCH --job-name=dropin_janich_nk
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/dropin_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/dropin_%A_%a.err
set -euo pipefail
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
VARIANTS=(real pred_obsdepth pred_preddepth); V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
OUT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/o2_liver/janich_nokozak/dropin; mkdir -p "$OUT"
export PYTHONNOUSERSITE=1
"$PY" /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/ribocode_dropin.py   --profiles /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/o2_liver/janich_nokozak/pred_profiles.npz   --annot $ECH/data/ribocode_annot_mouse --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
