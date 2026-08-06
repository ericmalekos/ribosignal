#!/usr/bin/env bash
#SBATCH --job-name=ho_predict_cpu
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/predictcpu_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/predictcpu_%A_%a.err
set -euo pipefail
# CPU fallback for the held-out predict stage (GPU partition 48/48 booked cluster-wide; the model
# is ~5M params so CPU is viable, ~3.5h/dump + ~3.5h/localization). Same env interface as the GPU
# runner: DATASET (required), ORFTRACK, FASTA_OVERRIDE, TX2CDS.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
DATASET=${DATASET:?set DATASET=<name> via --export}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2.npy}
FASTA_OVERRIDE=${FASTA_OVERRIDE:-}
TX2CDS=${TX2CDS:-}

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
echo "dataset=$DATASET backend=$BK run=$RUN host=$(hostname) CPU start=$(date -Iseconds)"

FASTA_ARG=(); ENV_FASTA=()
if [ -n "$FASTA_OVERRIDE" ]; then
  FASTA_ARG=(--fasta "$FASTA_OVERRIDE"); ENV_FASTA=(--env RIBO_ONEHOT_FASTA="$FASTA_OVERRIDE")
fi
ENV_TX2CDS=()
if [ -n "$TX2CDS" ]; then ENV_TX2CDS=(--env RIBO_TX2CDS="$TX2CDS"); fi

# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT/dropin" "${FASTA_ARG[@]}" --device cpu
echo "dump done=$(date -Iseconds)"

# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" "${ENV_TX2CDS[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT" --min_orf_nt 30 --eval_cap 0 "${FASTA_ARG[@]}" --device cpu
echo "dataset=$DATASET backend=$BK done=$(date -Iseconds)"
