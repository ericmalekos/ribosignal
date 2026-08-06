#!/usr/bin/env bash
#SBATCH --job-name=ribo_loc_cpu
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/loccpu_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/loccpu_%A_%a.err
set -euo pipefail
# CPU fallback for the onehot localization eval (GPU saturated cluster-wide; tiny model -> CPU OK).
RUNS=(
  results/loto/baseline_onehot_holdout_Hepatocytes
  results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
)
RUN=${RUNS[$SLURM_ARRAY_TASK_ID]}
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
mkdir -p "$NEW/logs/eval"
echo "run=$RUN host=$(hostname) CPU start=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env OMP_NUM_THREADS=16 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$NEW/$RUN" --min_orf_nt 30 --eval_cap 0 --device cpu
echo "run=$RUN done=$(date -Iseconds)"
