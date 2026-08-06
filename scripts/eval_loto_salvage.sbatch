#!/usr/bin/env bash
#SBATCH --job-name=loto_salvage
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_salvage_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_salvage_%A_%a.err
set -euo pipefail
# Salvage: run eval_extra on the partially-trained best.pt from the failed pilot (35180395, died on
# a filesystem-quota write failure ~epoch 5-8, val ~0.47) to get a PRELIMINARY held-out Hepatocytes
# number AND validate the eval_extra LOTO branch on real data before the full re-run's eval step.
TASKS=(baseline orf_v2_attn)
CFG=${TASKS[$SLURM_ARRAY_TASK_ID]}
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track.npy
[ "$CFG" = orf_v2_attn ] && ORFTRACK=$NEW/data/packed/orf_track_v2.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/loto/${CFG}_holdout_Hepatocytes
echo "cfg=$CFG host=$(hostname) start=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 --device cuda
echo "cfg=$CFG done=$(date -Iseconds)"
