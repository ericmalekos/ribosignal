#!/usr/bin/env bash
#SBATCH --job-name=ribo_dumpprof
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=01:30:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dump_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dump_%j.err
set -euo pipefail
# Task 14 step 1: dump per-nt PREDICTED profiles (+ real counts + count-head predicted total) for the
# held-out-Hepatocytes orf_v2_attn LOTO model, so RiboCode can be run on them (ribocode_dropin.py).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
RUN=${RUN:-results/loto/orf_v2_attn_holdout_Hepatocytes}
ORFTRACK=$NEW/data/packed/orf_track_v2.npy
SIF=${SIF:-/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img}
mkdir -p "$NEW/logs/dropin"

echo "run=$RUN host=$(hostname) start=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" \
    --run "$NEW/$RUN" --device cuda
echo "done=$(date -Iseconds)"
ls -la "$NEW/$RUN/dropin/"
