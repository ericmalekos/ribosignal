#!/usr/bin/env bash
#SBATCH --job-name=ho_dropin
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/dropin_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/dropin_%A_%a.err
set -euo pipefail
# RiboCode drop-in on the external held-out predicted profiles: 3 density variants (real /
# pred_obsdepth / pred_preddepth) through detectORF, then compare predicted calls against the
# held-out RiboCode calls. One array task per backend. DATASET via --export=ALL,DATASET=<name>.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
DATASET=${DATASET:?set DATASET=<name> via --export}
ANNOT=${ANNOT:-/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary}                    # v49 human; mouse -> ribocode_annot_mouse
OFFICIAL=$NEW/data/heldout_psites/$DATASET/${DATASET}_collapsed.txt
export PYTHONNOUSERSITE=1

BACKENDS=(onehot orthrus rinalmo)
BK=${BACKENDS[$SLURM_ARRAY_TASK_ID]}
OUT=$NEW/results/heldout/$DATASET/$BK
DROPIN=$OUT/dropin
PROFILES=$DROPIN/pred_profiles.npz
echo "dataset=$DATASET backend=$BK host=$(hostname) start=$(date -Iseconds)"
[ -f "$PROFILES" ] || { echo "missing $PROFILES (predict stage did not run?)"; exit 1; }

for V in real pred_obsdepth pred_preddepth; do
  echo "-- variant $V --"
  "$PY" "$NEW/scripts/ribocode_dropin.py" \
    --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$DROPIN" --min_aa 5 --pval 0.05
done

echo "-- compare vs held-out RiboCode calls --"
"$PY" "$NEW/scripts/compare_dropin_calls.py" \
  --dropin_dir "$DROPIN" --official "$OFFICIAL" --profiles "$PROFILES" \
  --key genomic --min_len 90 --min_enrichment 0.5 --out "$DROPIN"
echo "dataset=$DATASET backend=$BK done=$(date -Iseconds)"
