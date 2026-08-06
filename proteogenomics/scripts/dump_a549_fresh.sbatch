#!/usr/bin/env bash
#SBATCH --job-name=dump_a549f
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_a549f_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_a549f_%j.err
set -euo pipefail
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
# Defaults reproduce the ORIGINAL run (old deployed attn + heuristic-Kozak track). Override for the
# released models -- see dump_line.sbatch for the same contract:
#   MODEL_RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes \
#   SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img KOZAK=nokozak OUTSUB=pred_fresh_mamba4
SIF=${SIF:-$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img}
RUN=${MODEL_RUN:-results/loto/orf_v2_attn_onehot_holdout_Hepatocytes}
KOZAK=${KOZAK:-heuristic}
OUTSUB=${OUTSUB:-pred_fresh}
PILOT=$NEW/proteogenomics/data/A549_pilot
case "$KOZAK" in
  heuristic) ORFTRACK=$NEW/data/packed_heldout_human_a549_fresh/orf_track_v2.npy ;;
  nokozak)   ORFTRACK=$NEW/data/packed_heldout_human_a549_fresh/orf_track_v2_nokozak.npy ;;
  *) echo "KOZAK must be heuristic|nokozak, got '$KOZAK'" >&2; exit 2 ;;
esac
# The bare filename is NOT authoritative (packed_union carries nokozak content under it) -- check meta.
META=${ORFTRACK%.npy}_meta.json
if [[ -s "$META" ]]; then
  want=$([[ "$KOZAK" == nokozak ]] && echo none || echo heuristic)
  grep -q "\"kozak\": *\"$want\"" "$META" || { echo "ABORT: $META != kozak=$want" >&2; exit 1; }
fi
echo "host=$(hostname) MODEL_RUN=$RUN SIF=$(basename "$SIF") KOZAK=$KOZAK start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$PILOT/A549_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$RUN" \
    --heldout human_a549_fresh --tx_list "$PILOT/A549_universe_tx.txt" \
    --fasta "$PILOT/A549_universe.fa" --out "$PILOT/$OUTSUB" --budget 16000 --device cuda
echo "done=$(date -Iseconds)"
