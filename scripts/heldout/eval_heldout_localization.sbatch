#!/usr/bin/env bash
#SBATCH --job-name=ho_localiz
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/localiz_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/localiz_%A_%a.err
set -euo pipefail
# Full-set localization eval on an EXTERNAL held-out dataset (the localization-only half of
# eval_heldout_predict.sbatch -- the predicted profiles already exist, and eval_localization.py re-runs
# inference itself, so the dump step is redundant). One array task per sequence backend. Pass DATASET
# via --export=ALL,DATASET=<name>. Writes results/heldout/<DATASET>/<backend>/localization_metrics.json.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
DATASET=${DATASET:?set DATASET=<name> via --export}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2.npy}
FASTA_OVERRIDE=${FASTA_OVERRIDE:-}          # empty for human (reuse fibroblast_universe.fa)
TX2CDS=${TX2CDS:-}                           # empty for human (default v49 tx2cds.tsv); mouse -> vM38

BACKENDS=(onehot orthrus rinalmo)
RUNS=(
  results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
  results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes
  results/loto/orf_v2_attn_holdout_Hepatocytes
)
BK=${BACKENDS[$SLURM_ARRAY_TASK_ID]}
RUN=${RUNS[$SLURM_ARRAY_TASK_ID]}
OUT=$NEW/results/heldout/$DATASET/$BK
mkdir -p "$OUT" "$NEW/logs/heldout"
echo "dataset=$DATASET backend=$BK run=$RUN host=$(hostname) start=$(date -Iseconds)"

FASTA_ARG=(); ENV_FASTA=()
if [ -n "$FASTA_OVERRIDE" ]; then
  FASTA_ARG=(--fasta "$FASTA_OVERRIDE"); ENV_FASTA=(--env RIBO_ONEHOT_FASTA="$FASTA_OVERRIDE")
fi
ENV_TX2CDS=()
if [ -n "$TX2CDS" ]; then ENV_TX2CDS=(--env RIBO_TX2CDS="$TX2CDS"); fi

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" "${ENV_TX2CDS[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT" --min_orf_nt 30 --eval_cap 0 "${FASTA_ARG[@]}" --device cuda
echo "dataset=$DATASET backend=$BK done=$(date -Iseconds)"
