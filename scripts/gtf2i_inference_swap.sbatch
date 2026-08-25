#!/usr/bin/env bash
#SBATCH --job-name=gtf2i_swap
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/gtf2i_swap_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/gtf2i_swap_%j.err
set -euo pipefail
# INPUT-SWAP TEST, no retraining. Take the EXISTING mm1-trained mamba4 model and run it twice on the
# 35 GTF2I transcripts, changing only which coverage the pack supplies (mm1 vs mm10).
#
# This decides the memorisation question directly and for free. If the model fills the 1,617 nt gap
# when handed mm10 coverage, then it reads its RNA input and the mm1 blank was a corrupted input, not
# a memorised output -- and the fix works at INFERENCE time, without the 14 h retrain.
# If it still paints the gap empty on mm10 coverage, the retrain is doing something else entirely.
#
# mamba_ssm is CUDA-only (causal_conv1d_fwd asserts x.is_cuda), hence a GPU even for 35 transcripts.
R=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
RUN=${RUN:-$R/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes}
G=$R/data/mm25_diagnostic/gtf2i_inference
for SUF in union mm10cov; do
  O=$G/mamba4_dump_$SUF; mkdir -p "$O"
  echo "=== mamba4 mm1-trained, coverage=$SUF  $(date -Iseconds) ==="
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env RIBO_PACK_SUFFIX=$SUF \
    --env RIBO_ORF_TRACK="$R/data/packed_union/orf_track_v2.npy" \
    --env RIBO_ONEHOT_FASTA="$R/data/union_universe.fa" \
    -B "$BASE" "$SIF" \
    python3 "$R/scripts/dump_pred_profiles.py" --run "$RUN" --out "$O" \
      --tx_list "$G/gtf2i_tx.txt" --fasta "$R/data/union_universe.fa" \
      --min_signal 0 --device cuda
done
echo "GTF2I SWAP DONE $(date -Iseconds)"
