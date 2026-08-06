#!/usr/bin/env bash
#SBATCH --job-name=sweep_pval
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-9
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_pval_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_pval_%A_%a.err
set -euo pipefail
RUN=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
PVALS=(0.05 1e-4 1e-8 1e-15 1e-25 1e-40 1e-60 1e-90 1e-130 1e-200)
P=${PVALS[$SLURM_ARRAY_TASK_ID]}
OUT=$RUN/dropin_sweep_pval/theta_${P}; rm -rf "$OUT"; mkdir -p "$OUT"
export PYTHONNOUSERSITE=1
# lever (d): fixed predicted depth (pred_scale 1.0), tighten the RiboCode p-value threshold
"$PY" /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/ribocode_dropin.py --profiles "$RUN/dropin/pred_profiles.npz"   --annot /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human/data/ribocode_annot --variant pred_preddepth --pred_scale 1.0   --out "$OUT" --min_aa 5 --pval "$P"
