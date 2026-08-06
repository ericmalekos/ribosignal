#!/usr/bin/env bash
#SBATCH --job-name=eval_attn
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/eval_attn_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/eval_attn_%j.err
set -euo pipefail

# Re-run eval_extra (count-head magnitude + posture-B) for the 3 attention runs, whose
# in-job eval_extra OOMed at the old hardcoded budget 24000 (attention is O(L^2); they train
# at 16000). eval_extra now caps its budget at the run's training budget. orf_v2_attn must
# point at orf_track_v2.npy so its ORF features match training.
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model

for r in orf_attn orf_v2_attn attn; do
  OT=""
  [ "$r" = "orf_v2_attn" ] && OT="$NEW/data/packed/orf_track_v2.npy"
  echo "=== eval_extra $r (orf_track=${OT:-default}) $(date -Iseconds) ==="
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    ${OT:+--env RIBO_ORF_TRACK=$OT} \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/eval_extra.py" --run "$NEW/results/improve/${r}_f0" --eval_cap 0 \
    || echo "WARN eval_extra failed for $r"
done
echo "=== eval_attn done $(date -Iseconds) ==="
