#!/usr/bin/env bash
#SBATCH --job-name=ho_finalize
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:55:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/finalize_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/finalize_%A_%a.err
set -euo pipefail
# Merge the dump shards for one (DATASET, BACKEND), then run the RiboCode drop-in (3 density variants
# + compare against the held-out calls). Chained afterok the shard array. DATASET/BACKEND/NSHARDS/ANNOT.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PYC=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PYR=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
DATASET=${DATASET:?}; BACKEND=${BACKEND:?}; NSHARDS=${NSHARDS:?}
ANNOT=${ANNOT:-/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary}
OFFICIAL=$NEW/data/heldout_psites/$DATASET/${DATASET}_collapsed.txt
DROPIN=$NEW/results/heldout/$DATASET/$BACKEND/dropin
PROFILES=$DROPIN/pred_profiles.npz
export PYTHONNOUSERSITE=1
echo "finalize dataset=$DATASET backend=$BACKEND nshards=$NSHARDS host=$(hostname) start=$(date -Iseconds)"

"$PYC" "$NEW/scripts/heldout/merge_pred_shards.py" "$DROPIN" "$NSHARDS"
[ -f "$PROFILES" ] || { echo "merge failed, no $PROFILES"; exit 1; }

for V in real pred_obsdepth pred_preddepth; do
  echo "-- variant $V --"
  "$PYR" "$NEW/scripts/ribocode_dropin.py" \
    --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$DROPIN" --min_aa 5 --pval 0.05
done
echo "-- compare --"
"$PYR" "$NEW/scripts/compare_dropin_calls.py" \
  --dropin_dir "$DROPIN" --official "$OFFICIAL" --profiles "$PROFILES" \
  --key genomic --min_len 90 --min_enrichment 0.5 --out "$DROPIN"
echo "finalize $DATASET/$BACKEND done=$(date -Iseconds)"
