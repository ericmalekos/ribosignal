#!/usr/bin/env bash
#SBATCH --job-name=ho_dumpshard
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=00:55:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/dumpshard_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/dumpshard_%A_%a.err
set -euo pipefail
# One shard of the held-out predicted-profile dump on the short partition (CPU). The full dump of
# an orf_v2_attn model over all test tx exceeds short's 1h wall (the 2-layer transformer's O(L^2)
# attention on long transcripts is the cost); interleaved sharding makes each piece fit, then
# merge_pred_shards.py stitches them. array index = shard. DATASET + BACKEND + NSHARDS via --export.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
DATASET=${DATASET:?set DATASET=<name>}
BACKEND=${BACKEND:?set BACKEND=onehot|orthrus|rinalmo}
NSHARDS=${NSHARDS:?set NSHARDS}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2.npy}
FASTA_OVERRIDE=${FASTA_OVERRIDE:-}
case "$BACKEND" in
  onehot)  RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes ;;
  orthrus) RUN=results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes ;;
  rinalmo) RUN=results/loto/orf_v2_attn_holdout_Hepatocytes ;;
  *) echo "bad BACKEND $BACKEND"; exit 1 ;;
esac
OUT=$NEW/results/heldout/$DATASET/$BACKEND
mkdir -p "$OUT/dropin" "$NEW/logs/heldout"
echo "dump shard $SLURM_ARRAY_TASK_ID/$NSHARDS dataset=$DATASET backend=$BACKEND host=$(hostname) start=$(date -Iseconds)"

FASTA_ARG=(); ENV_FASTA=()
if [ -n "$FASTA_OVERRIDE" ]; then FASTA_ARG=(--fasta "$FASTA_OVERRIDE"); ENV_FASTA=(--env RIBO_ONEHOT_FASTA="$FASTA_OVERRIDE"); fi

# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT/dropin" "${FASTA_ARG[@]}" \
    --nshards "$NSHARDS" --shard "$SLURM_ARRAY_TASK_ID" --device cpu
echo "dump shard $SLURM_ARRAY_TASK_ID done=$(date -Iseconds)"
