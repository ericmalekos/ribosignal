#!/usr/bin/env bash
#SBATCH --job-name=eval_noBrain
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/eval_noBrain_mm20_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/eval_noBrain_mm20_%j.err
set -euo pipefail
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env RIBO_ORF_TRACK="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/packed/orf_track_v2.npy" \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta "/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img" \
  python3 "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/eval_extra.py" --run "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_noBrain_mm20_holdout_Hepatocytes" --eval_cap 0 --device cuda
