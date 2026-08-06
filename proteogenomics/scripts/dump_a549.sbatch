#!/usr/bin/env bash
#SBATCH --job-name=dump_a549
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_a549_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_a549_%j.err
set -euo pipefail
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes   # deployment onehot orf_v2_attn (trained on orf_track_v2)
ORFTRACK=$NEW/data/packed/orf_track_v2.npy                # the track this model trained with
echo "host=$(hostname) start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$RUN" \
    --heldout human_a549_encode --tx_list "$NEW/proteogenomics/data/A549_pilot/A549_score_tx.txt" \
    --out "$NEW/proteogenomics/data/A549_pilot/pred" --device cuda
echo "done=$(date -Iseconds)"
