#!/usr/bin/env bash
#SBATCH --job-name=orf_tcell_nk
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/orf_tcell_nk_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/orf_tcell_nk_%j.err
set -euo pipefail
# NO-KOZAK ORF track for the GSE155087 T-cell held-out universe, matching the deployed mm1 no-kozak model's
# training input. Run AFTER build_heldout_pack (needs the pack's tx_order/offsets to align the track).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/logs/orf"
PACK=$NEW/data/packed_heldout_mouse_gse155087_tcell
FASTA=$NEW/data/heldout_refs/mouse_gse155087_tcell_universe.fa
"$PY" "$NEW/scripts/build_orf_track.py" --mode ext --kozak none --pack "$PACK" --fasta "$FASTA" \
  --out "$PACK/orf_track_v2_nokozak.npy"
echo "no-kozak orf track -> $PACK/orf_track_v2_nokozak.npy"
