#!/usr/bin/env bash
#SBATCH --job-name=ribo_dropfloor
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/floor_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/floor_%A_%a.err
set -euo pipefail
# Task 14 FP-reduction: re-run the drop-in with a per-position density floor (zero predicted positions
# below floor_mult x uniform) to kill the spurious 3'UTR dORF over-calls. Pred variants only (real is
# untouched). Output tagged <variant>_floor<mult>.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=${RUN:-results/loto/orf_v2_attn_holdout_Hepatocytes}
PROFILES=$NEW/$RUN/dropin/pred_profiles.npz
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
OUT=$NEW/$RUN/dropin
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin" "$OUT"

RUNS=("pred_obsdepth 0.5" "pred_obsdepth 1.0" "pred_preddepth 0.5" "pred_preddepth 1.0")
read -r V F <<< "${RUNS[$SLURM_ARRAY_TASK_ID]}"
echo "variant=$V floor_mult=$F host=$(hostname) start=$(date -Iseconds)"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$OUT" \
  --min_aa 5 --pval 0.05 --floor_mult "$F"
echo "variant=$V floor_mult=$F done=$(date -Iseconds)"
