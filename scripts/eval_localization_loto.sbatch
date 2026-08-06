#!/usr/bin/env bash
#SBATCH --job-name=ribo_loc_loto
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/locloto_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/locloto_%A_%a.err
set -euo pipefail
# Cross-tissue localization eval: run eval_localization.py on the held-out-Hepatocytes LOTO models,
# scored against Hepatocytes' OWN RiboCode ORF calls. Tests whether the translated-ORF localization
# survives cross-tissue transfer, not just within Fibroblast. Array indexes into the completed LOTO
# runs; add the Orthrus runs to RUNS once their training + eval_extra finish.
RUNS=(
  results/loto/baseline_holdout_Hepatocytes
  results/loto/orf_v2_attn_holdout_Hepatocytes
  results/loto/baseline_orthrus_holdout_Hepatocytes
  results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes
  results/loto/baseline_onehot_holdout_Hepatocytes
  results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
)
RUN=${RUNS[$SLURM_ARRAY_TASK_ID]}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2.npy   # baseline ignores it (use_orf_track=false)
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
mkdir -p "$NEW/logs/eval"

echo "run=$RUN host=$(hostname) start=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" \
    --run "$NEW/$RUN" --min_orf_nt 30 --eval_cap 0 --device cuda
echo "run=$RUN done=$(date -Iseconds)"
