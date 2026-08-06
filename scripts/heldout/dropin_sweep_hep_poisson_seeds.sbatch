#!/usr/bin/env bash
#SBATCH --job-name=sweep_pois_seed
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-19
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_pois_seed_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_pois_seed_%A_%a.err
set -euo pipefail
RUN=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
THETAS=(0.02 0.05 0.1 0.2 0.35 0.5 0.7 1.0 1.5 2.5)
SEED=$(( 1 + SLURM_ARRAY_TASK_ID / 10 ))
TH=${THETAS[$(( SLURM_ARRAY_TASK_ID % 10 ))]}
OUT=$RUN/dropin_sweep_poisson_seed${SEED}/theta_${TH}; rm -rf "$OUT"; mkdir -p "$OUT"
export PYTHONNOUSERSITE=1
"$PY" /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/ribocode_dropin.py --profiles "$RUN/dropin/pred_profiles.npz"   --annot $ECH/data/ribocode_annot --variant pred_preddepth --pred_scale "$TH" --pred_poisson   --seed "$SEED" --out "$OUT" --min_aa 5 --pval 0.05
