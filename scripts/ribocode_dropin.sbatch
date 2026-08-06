#!/usr/bin/env bash
#SBATCH --job-name=ribo_dropin
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dropin_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dropin_%A_%a.err
set -euo pipefail
# Task 14 step 2: run RiboCode's detectORF on the substituted density (3 variants, one per array task).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=${RUN:-results/loto/orf_v2_attn_holdout_Hepatocytes}
PROFILES=$NEW/$RUN/dropin/pred_profiles.npz
ANNOT=$ECH/data/ribocode_annot
OUT=$NEW/$RUN/dropin
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin" "$OUT"

VARIANTS=(real pred_obsdepth pred_preddepth)
V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
echo "variant=$V host=$(hostname) start=$(date -Iseconds)"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
echo "variant=$V done=$(date -Iseconds)"
