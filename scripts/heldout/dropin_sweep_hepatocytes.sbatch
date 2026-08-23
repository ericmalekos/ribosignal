#!/usr/bin/env bash
#SBATCH --job-name=dropin_sweep
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-9
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_hep_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/sweep_hep_%A_%a.err
set -euo pipefail
# Check 2: CDS-anchored operating curve on the high-depth Hepatocytes holdout (mm1 no-kozak model).
# Sweep the effective-depth knob theta = ribocode_dropin --subsample (theta<1 thins the predicted
# density -> lower effective depth -> fewer, higher-confidence calls; theta=1.0 = raw over-calling).
# Thinning is Poisson-like, so it also injects realistic detection noise into the smooth profile.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary          # human (Hepatocytes)
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin"
THETAS=(0.02 0.05 0.1 0.2 0.35 0.5 0.7 1.0 1.5 2.5)
TH=${THETAS[$SLURM_ARRAY_TASK_ID]}
OUT=$RUN/dropin_sweep/theta_${TH}
rm -rf "$OUT"; mkdir -p "$OUT"
echo "theta(pred_scale)=$TH start=$(date -Iseconds)"
export PYTHONNOUSERSITE=1
# pred_preddepth (standalone) at effective depth = predicted_total * theta (the CDS-anchored knob)
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$RUN/dropin/pred_profiles.npz" --annot "$ANNOT" --variant pred_preddepth \
  --pred_scale "$TH" --seed 0 --out "$OUT" --min_aa 5 --pval 0.05
echo "theta=$TH done -> $OUT/pred_preddepth_collapsed.txt"
