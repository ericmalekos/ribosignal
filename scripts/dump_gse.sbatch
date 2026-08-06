#!/usr/bin/env bash
#SBATCH --job-name=dump_gse
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/dump_gse_%x_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/dump_gse_%x_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
G=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/gse120762
PACK=$NEW/data/packed_heldout_mouse_gse120762_${COND:?}
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2.npy" --env RIBO_ONEHOT_FASTA="$G/gse120762_${COND}_universe.fa" \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta "/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img" \
  python3 $NEW/scripts/dump_pred_profiles.py --run "$NEW/results/loto/orf_v2_attn_onehot_holdout_Hepatocytes" \
    --heldout mouse_gse120762_${COND} --tx_list "$G/gse120762_${COND}_universe_tx.txt" \
    --fasta "$G/gse120762_${COND}_universe.fa" --out "$G/pred_${COND}" --device cuda
echo "dump $COND done -> $G/pred_${COND}/pred_profiles.npz"
