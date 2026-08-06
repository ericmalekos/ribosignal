#!/usr/bin/env bash
#SBATCH --job-name=dump_janich_nk
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/dump_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/dump_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RUN=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
PACK=$NEW/data/packed_heldout_mouse_janich_liver
REF=$NEW/data/heldout_refs
# NOTE: mouse_janich_universe_tx.txt is the REBUILT expressed universe (10,431 tx, TPM>=1),
# not the broken 2,434-tx one -- see data/_quarantine/janich_stale_universe_2026-07-30/README.md
OUT=$NEW/results/o2_liver/janich_nokozak; mkdir -p "$OUT"
# The model is the NO-KOZAK mm1 model, so it must be fed the NO-KOZAK ORF track: train/inference must
# match (standing rule). This previously pointed at orf_track_v2.npy (the Kozak-heuristic track), which
# was a silent mismatch -- the nokozak track for this pack did not exist until 2026-07-30.
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2_nokozak.npy" \
  --env RIBO_ONEHOT_FASTA="$REF/mouse_janich_universe.fa" \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta \
  "/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img" \
  python3 $NEW/scripts/dump_pred_profiles.py --run "$RUN" \
    --heldout mouse_janich_liver --tx_list "$REF/mouse_janich_universe_tx.txt" \
    --fasta "$REF/mouse_janich_universe.fa" --out "$OUT" --device cuda
echo "janich dump -> $OUT/pred_profiles.npz"
