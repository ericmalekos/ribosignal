#!/usr/bin/env bash
#SBATCH --job-name=orf_wang_nk
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/orf_wang_nk_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/orf_wang_nk_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PACK=$NEW/data/packed_heldout_mouse_wang_liver
FASTA=$NEW/data/heldout_refs/mouse_wang_universe.fa
$PY $NEW/scripts/build_orf_track.py --mode ext --kozak none --pack "$PACK" --fasta "$FASTA" \
  --out "$PACK/orf_track_v2_nokozak.npy"
echo "wang no-kozak orf track -> $PACK/orf_track_v2_nokozak.npy"
