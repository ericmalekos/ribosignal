#!/usr/bin/env bash
#SBATCH --job-name=xsp_dropin
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/dropin_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/dropin_%A_%a.err
#
#   sbatch --array=0-2 --export=ALL,DS=xsp_yeast,ARCH=attn,SP=yeast scripts/xspecies/dropin_xspecies.sbatch
#
# RiboCode run on the MODEL'S predicted profile instead of the observed one, so predicted
# and observed ORF calls come out of the same caller under the same statistics. Three
# variants, per scripts/heldout/dropin_*.sbatch:
#   real            the observed profile (the reference call set)
#   pred_obsdepth   predicted SHAPE scaled to the observed depth -- isolates whether the
#                   model puts the ribosomes in the right places
#   pred_preddepth  predicted shape AND the model's own count head -- the fully de novo call
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
DS=${DS:?set DS}; ARCH=${ARCH:?set ARCH}; SP=${SP:?set SP}
VARIANTS=(real pred_obsdepth pred_preddepth); V=${VARIANTS[${SLURM_ARRAY_TASK_ID:-0}]}
PROF=$NEW/results/loto/orf_v2_${ARCH}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin_${DS}/pred_profiles.npz
ANNOT=$NEW/data/ribocode_annot_${SP}
OUT=$NEW/results/xspecies_dropin/${DS}_${ARCH}; mkdir -p "$OUT"
[[ -s "$PROF" ]] || { echo "missing $PROF" >&2; exit 1; }
[[ -d "$ANNOT" ]] || { echo "missing $ANNOT" >&2; exit 1; }
export PYTHONNOUSERSITE=1
echo "ds=$DS arch=$ARCH species=$SP variant=$V start=$(date -Iseconds)"
"$PY" "$NEW/scripts/ribocode_dropin.py" --profiles "$PROF" --annot "$ANNOT" \
  --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
n=$(( $(wc -l < "$OUT/${V}_collapsed.txt" 2>/dev/null || echo 1) - 1 ))
echo "=== $DS/$ARCH/$V -> $n calls  DONE=$(date -Iseconds)"
