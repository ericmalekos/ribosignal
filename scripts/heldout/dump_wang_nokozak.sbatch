#!/usr/bin/env bash
#SBATCH --job-name=dump_wang_nk
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/dump_wang_nk_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/dump_wang_nk_%j.err
set -euo pipefail
# Predict on the EXPRESSED Wang universe (13,041 tx), not the raw 221,835-tx annotation.
#
# `--tx_list` bypasses the `min_signal` gate that every other eval applies. Before 2026-07-31 this pointed
# at mouse_wang_universe.txt, which is essentially all of vM38 unfiltered, because Wang had no salmon quant
# to filter with (its RNA was turned into coverage and the fastqs deleted). The model therefore called ORFs
# across ~205k transcripts the reference Ribo-seq never covered, and every one counted as a false positive:
# over-call 6.1x, and O2 non-canonical read 13% of ceiling instead of ~50%. Wang RNA was re-downloaded and
# quantified on 2026-07-31 (91-93% mapping), so a real expressed universe now exists. The pack still holds
# all 221,835 tx; --tx_list just selects the expressed subset to predict on.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RUN=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
PACK=$NEW/data/packed_heldout_mouse_wang_liver
REF=$NEW/data/heldout_refs
TXLIST=$REF/mouse_wang_expressed_universe_tx.txt      # 13,041 expressed (TPM>=1); NOT the raw annotation
OUT=$NEW/results/o2_liver/wang_nokozak; mkdir -p "$OUT"
[ -s "$TXLIST" ] || { echo "ABORT: missing $TXLIST (run quant_wang_rna.sbatch first)"; exit 1; }
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2_nokozak.npy" \
  --env RIBO_ONEHOT_FASTA="$REF/mouse_wang_universe.fa" \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta "/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img" \
  python3 $NEW/scripts/dump_pred_profiles.py --run "$RUN" \
    --heldout mouse_wang_liver --tx_list "$TXLIST" \
    --fasta "$REF/mouse_wang_universe.fa" --out "$OUT" --device cuda
echo "wang dump (no-kozak model) -> $OUT/pred_profiles.npz"
