#!/usr/bin/env bash
#SBATCH --job-name=ribo_dump_cpu
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dump_cpu_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dump_cpu_%j.err
set -euo pipefail
# CPU fallback for dump_pred_profiles (Task 14 step 1) -- the GPU partition is saturated cluster-wide
# and the onehot model is tiny (~5M params), so a CPU forward over the held-out test tx is feasible.
# Same script, no --nv, --device cpu.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
RUN=${RUN:-results/loto/orf_v2_attn_onehot_holdout_Hepatocytes}
ORFTRACK=$NEW/data/packed/orf_track_v2.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
mkdir -p "$NEW/logs/dropin"
echo "run=$RUN host=$(hostname) CPU start=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env OMP_NUM_THREADS=16 \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --device cpu
echo "done=$(date -Iseconds)"
ls -la "$NEW/$RUN/dropin/"
