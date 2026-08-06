#!/usr/bin/env bash
#SBATCH --job-name=ho_prelim
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=00:50:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/prelim_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/prelim_%A_%a.err
set -euo pipefail
# Fast preliminary: capped one-hot localization on the SHORT partition (1h, backfills into the
# reserved window while the full 4-way eval waits in the medium queue). Gives an early cross-study
# non-canonical AUROC read on ~EVAL_CAP held-out transcripts. DATASET + EVAL_CAP + optional
# ORFTRACK/FASTA_OVERRIDE/TX2CDS via --export.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
DATASET=${DATASET:?set DATASET=<name>}
EVAL_CAP=${EVAL_CAP:-4000}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2.npy}
FASTA_OVERRIDE=${FASTA_OVERRIDE:-}
TX2CDS=${TX2CDS:-}
RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
OUT=$NEW/results/heldout/$DATASET/onehot_prelim
mkdir -p "$OUT" "$NEW/logs/heldout"
echo "PRELIM dataset=$DATASET cap=$EVAL_CAP host=$(hostname) start=$(date -Iseconds)"

FASTA_ARG=(); ENV_FASTA=(); ENV_TX2CDS=()
if [ -n "$FASTA_OVERRIDE" ]; then FASTA_ARG=(--fasta "$FASTA_OVERRIDE"); ENV_FASTA=(--env RIBO_ONEHOT_FASTA="$FASTA_OVERRIDE"); fi
if [ -n "$TX2CDS" ]; then ENV_TX2CDS=(--env RIBO_TX2CDS="$TX2CDS"); fi

# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" "${ENV_TX2CDS[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT" --min_orf_nt 30 --eval_cap "$EVAL_CAP" "${FASTA_ARG[@]}" --device cpu
echo "PRELIM done=$(date -Iseconds)"
