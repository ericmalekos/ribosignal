#!/usr/bin/env bash
#SBATCH --job-name=dump_hbl1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_hbl1_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_hbl1_%j.err
set -euo pipefail
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
PILOT=$NEW/proteogenomics/data/HBL1_pilot
ORFTRACK=$NEW/data/packed_heldout_human_hbl1/orf_track_v2.npy   # fresh, heuristic (matches model training)
echo "host=$(hostname) start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$PILOT/HBL1_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$RUN" \
    --heldout human_hbl1 --tx_list "$PILOT/HBL1_universe_tx.txt" \
    --fasta "$PILOT/HBL1_universe.fa" --out "$PILOT/pred" --budget 16000 --device cuda
echo "done=$(date -Iseconds)"
