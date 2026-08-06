#!/usr/bin/env bash
#SBATCH --job-name=dump_peaky
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=01:30:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dump_peaky_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/dump_peaky_%j.err
set -euo pipefail
# O3 track-2 A/B: dump the PK=0.2 anti-smoothing model's predicted profiles on the Hepatocytes holdout
# (dump_pred_profiles.py reads holdout=Hepatocytes + emb_backend=onehot from the run's args.json), using the
# NO-KOZAK ORF track it was trained on. Feeds the drop-in -> orf_features -> fp_filter A/B vs the standard model.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
RUN=results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_peaky0p2_holdout_Hepatocytes
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
mkdir -p "$NEW/logs/dropin"
echo "run=$RUN host=$(hostname) start=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --device cuda
echo "done=$(date -Iseconds)"
ls -la "$NEW/$RUN/dropin/"
