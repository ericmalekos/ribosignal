#!/usr/bin/env bash
#SBATCH --job-name=eval_extra_union
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval_extra_union_%j.out
set -euo pipefail
# eval_extra on the DEPLOYED union model, to give figure B4 its per-stratum model numbers on the
# model that actually ships. B4 currently carries pre-mm1/pre-union values (Task 22, 2026-07-17).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
RUN=$NEW/results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_PACK_SUFFIX=union \
  --env RIBO_ORF_TRACK="$NEW/data/packed_union/orf_track_v2.npy" \
  --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$RUN" --eval_cap 0 --device cuda
echo "done=$(date -Iseconds)"
