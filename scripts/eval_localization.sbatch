#!/usr/bin/env bash
#SBATCH --job-name=ribo_loc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=110G
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/loc_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/loc_%A_%a.err
set -euo pipefail
# Localization eval (eval_localization.py): does the model concentrate predicted P-site signal on
# RiboCode-translated ORFs vs untranslated candidate AUG ORFs, on held-out transcripts? Runs the two
# core within-tissue fold-0 models to isolate the ORF-track + attention contribution. Inference only,
# so short. LOTO runs are evaluated separately once the Orthrus arm lands.
RUNS=(results/improve/baseline_f0 results/improve/orf_v2_attn_f0)
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
