#!/usr/bin/env bash
#SBATCH --job-name=sweep_gse_pois
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-19
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_gse_pois_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_gse_pois_%A_%a.err
set -euo pipefail
G=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/gse120762
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot
CONDS=(nt lps); THETAS=(0.02 0.05 0.1 0.2 0.35 0.5 0.7 1.0 1.5 2.5)
COND=${CONDS[$(( SLURM_ARRAY_TASK_ID / 10 ))]}
TH=${THETAS[$(( SLURM_ARRAY_TASK_ID % 10 ))]}
OUT=$G/pred_${COND}/dropin_poisson_sweep/theta_${TH}; rm -rf "$OUT"; mkdir -p "$OUT"
export PYTHONNOUSERSITE=1
"$PY" /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/ribocode_dropin.py --profiles "$G/pred_${COND}/pred_profiles.npz"   --annot "$ANNOT" --variant pred_preddepth --pred_scale "$TH" --pred_poisson --seed 0   --out "$OUT" --min_aa 5 --pval 0.05
