#!/usr/bin/env bash
#SBATCH --job-name=dropin_gse_nk
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --array=0-5
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dropin_gse_nk_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dropin_gse_nk_%A_%a.err
set -euo pipefail
# RiboCode drop-in on the GSE120762 predicted profiles from the mm1 NO-KOZAK model, per condition
# x 3 variants (real / pred_obsdepth / pred_preddepth). Mouse annotation.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
G=$NEW/proteogenomics/data/gse120762
ANNOT=$ECH/data/ribocode_annot_mouse
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin"
CONDS=(nt lps); VARIANTS=(real pred_obsdepth pred_preddepth)
COND=${CONDS[$((SLURM_ARRAY_TASK_ID / 3))]}
V=${VARIANTS[$((SLURM_ARRAY_TASK_ID % 3))]}
PROFILES=$G/pred_${COND}_nokozak/pred_profiles.npz
OUT=$G/pred_${COND}_nokozak/dropin
mkdir -p "$OUT"
echo "cond=$COND variant=$V (no-kozak model) start=$(date -Iseconds)"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
echo "cond=$COND variant=$V done -> $OUT/${V}_collapsed.txt"
