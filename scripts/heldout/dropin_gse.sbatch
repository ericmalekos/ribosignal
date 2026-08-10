#!/usr/bin/env bash
#SBATCH --job-name=dropin_gse
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --array=0-5
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dropin_gse_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dropin_gse_%A_%a.err
set -euo pipefail
# RiboCode drop-in on the GSE120762 predicted profiles (deployed heuristic model), per condition x
# 3 density variants. 'real' = RiboCode on the universe-restricted OBSERVED P-sites (the fair
# reference, same tx set as the predicted); 'pred_obsdepth' = predicted SHAPE * real depth;
# 'pred_preddepth' = fully standalone (predicted shape * count-head total, no real Ribo data).
# Compare the resulting ORF-type composition (predicted vs real) per condition and NT-vs-LPS.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
G=$NEW/proteogenomics/data/gse120762
ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"            # mouse: GSE120762 is BMDM
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin"

CONDS=(nt lps)
VARIANTS=(real pred_obsdepth pred_preddepth)
COND=${CONDS[$((SLURM_ARRAY_TASK_ID / 3))]}
V=${VARIANTS[$((SLURM_ARRAY_TASK_ID % 3))]}
PROFILES=$G/pred_${COND}/pred_profiles.npz
OUT=$G/pred_${COND}/dropin
mkdir -p "$OUT"
echo "cond=$COND variant=$V host=$(hostname) start=$(date -Iseconds)"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
echo "cond=$COND variant=$V done=$(date -Iseconds)  -> $OUT/${V}_collapsed.txt"
