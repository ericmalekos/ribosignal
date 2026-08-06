#!/usr/bin/env bash
#SBATCH --job-name=ho_predict
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/predict_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/predict_%A_%a.err
set -euo pipefail
# Apply a Chothani-trained orf_v2_attn model to an EXTERNAL held-out dataset: dump predicted
# per-nt profiles (for the drop-in) + run the localization eval against the held-out RiboCode calls.
# One array task per sequence backend. Pass DATASET via --export=ALL,DATASET=<name>.
# RIBO_ORF_TRACK -> the shared orf_track_v2.npy (held-out pack offsets are identical to the ref pack,
# so it aligns 1:1). RIBO_ONEHOT_FASTA optionally overrides the one-hot universe (cross-species).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
DATASET=${DATASET:?set DATASET=<name> via --export}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2.npy}
FASTA_OVERRIDE=${FASTA_OVERRIDE:-}          # empty for human (reuse fibroblast_universe.fa)
TX2CDS=${TX2CDS:-}                           # empty for human (default v49 tx2cds.tsv); mouse -> vM38

# backend -> trained LOTO-Hepatocytes orf_v2_attn run
BACKENDS=(onehot orthrus rinalmo)
RUNS=(
  results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
  results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes
  results/loto/orf_v2_attn_holdout_Hepatocytes
)
BK=${BACKENDS[$SLURM_ARRAY_TASK_ID]}
RUN=${RUNS[$SLURM_ARRAY_TASK_ID]}
OUT=$NEW/results/heldout/$DATASET/$BK
mkdir -p "$OUT/dropin" "$NEW/logs/heldout"
echo "dataset=$DATASET backend=$BK run=$RUN host=$(hostname) start=$(date -Iseconds)"

FASTA_ARG=()
ENV_FASTA=()
if [ -n "$FASTA_OVERRIDE" ]; then
  FASTA_ARG=(--fasta "$FASTA_OVERRIDE")
  ENV_FASTA=(--env RIBO_ONEHOT_FASTA="$FASTA_OVERRIDE")
fi
ENV_TX2CDS=()
if [ -n "$TX2CDS" ]; then
  ENV_TX2CDS=(--env RIBO_TX2CDS="$TX2CDS")
fi

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT/dropin" "${FASTA_ARG[@]}" --device cuda

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" "${ENV_TX2CDS[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT" --min_orf_nt 30 --eval_cap 0 "${FASTA_ARG[@]}" --device cuda

echo "dataset=$DATASET backend=$BK done=$(date -Iseconds)"
