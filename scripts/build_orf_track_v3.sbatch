#!/usr/bin/env bash
#SBATCH --job-name=orf_track_v3
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/build_v3_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/build_v3_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/logs/orf"
# ATG+CTG variant: occupancy over ATG-or-CTG-initiated ORFs + 2-level start channel (ATG=1.0,
# CTG=0.5). Tests the single dominant alternative start. Writes orf_track_v3.npy only (does not
# touch orf_track.npy / orf_track_v2.npy that the v1 / v2 runs use).
$PY "$NEW/scripts/build_orf_track.py" --mode atgctg --out "$NEW/data/packed/orf_track_v3.npy"
echo "orf v3 build done=$(date -Iseconds)"
