#!/usr/bin/env bash
#SBATCH --job-name=orf_track_v2
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/build_v2_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/build_v2_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/logs/orf"
# non-AUG-aware variant -> orf_track_v2.npy (does NOT touch orf_track.npy that v1 runs use).
# --kozak heuristic is PINNED here: this file trains the deployment orf_v2_attn model, so it must keep
# reproducing the original heuristic-Kozak track byte-for-byte. The build_orf_track.py DEFAULT is now
# --kozak none (results.md Task 20: the Kozak factor is redundant/mildly harmful), but flipping the
# default must NOT silently rebuild the deployed model's track -- hence the explicit heuristic here.
$PY "$NEW/scripts/build_orf_track.py" --mode ext --kozak heuristic --out "$NEW/data/packed/orf_track_v2.npy"
echo "orf v2 build done=$(date -Iseconds)"
