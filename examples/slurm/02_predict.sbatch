#!/usr/bin/env bash
#SBATCH --job-name=ribosignal_predict
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
# Predict the per-nucleotide profile on the GPU. Both architectures also run on CPU (mamba4
# through the pure-PyTorch reference, 6.3x slower than attn there). RIBO_MAMBA_IMPL=cuda makes
# a missing mamba_ssm an error rather than a silent switch to that far slower reference path.
: "${RIBO_IMAGE:?}"; : "${RIBO_WORK:?}"
ARCH=${ARCH:-mamba4}
singularity exec --nv --cleanenv -B "$RIBO_WORK" \
  --env RIBO_ONEHOT_FASTA="$RIBO_WORK/ref/transcripts.fa" \
  --env RIBO_ORF_TRACK="$RIBO_WORK/pack/demo/orf_track_v2_nokozak.npy" \
  --env RIBO_MAMBA_IMPL=cuda \
  "$RIBO_IMAGE" python /opt/riboseq-model/scripts/dump_pred_profiles.py \
    --run /opt/riboseq-model/weights --arch "$ARCH" \
    --pack "$RIBO_WORK/pack/demo" --out "$RIBO_WORK/pred/$ARCH" --device cuda
