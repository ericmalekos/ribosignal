#!/usr/bin/env bash
#SBATCH --job-name=dump_line
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_line_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dump_line_%j.err
set -euo pipefail
# Generalized model-predict (deployed orf_v2_attn one-hot, LOTO Hepatocyte). Pass via --export:
# LINE (dir name), HELDOUT (heldout dataset key, e.g. human_dohh2). Prereqs (built by the prep driver):
# <LINE>_universe_tx.txt, <LINE>_universe.fa, packed_heldout_<HELDOUT>/orf_track_v2.npy.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
#
# MODEL SELECTION. Defaults reproduce the ORIGINAL 2026-07-21 runs (old deployed orf_v2_attn +
# heuristic-Kozak track) so those results stay regenerable. Override for the RELEASED models:
#
#   released mamba4 (primary, needs the orthrus image -- rinalmo has no mamba_ssm):
#     MODEL_RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes \
#     SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img KOZAK=nokozak
#   released attn (CPU-capable):
#     MODEL_RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes \
#     KOZAK=nokozak
#
# KOZAK picks the ORF track, and the track MUST match what the model trained on: a nokozak-trained
# model fed the heuristic track is a train/inference mismatch (feedback_kozak_never_set).
SIF=${SIF:-$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img}
RUN=${MODEL_RUN:-results/loto/orf_v2_attn_onehot_holdout_Hepatocytes}
KOZAK=${KOZAK:-heuristic}
PILOT=$NEW/proteogenomics/data/${LINE}_pilot
OUTSUB=${OUTSUB:-pred}
# UNIV / DEVICE / MAMBA_CPU are overridable; defaults reproduce the original behaviour exactly.
#   UNIV=<prefix>  selects <PILOT>/<prefix>_universe{_tx.txt,.fa}. The full-depth all-run rebuild
#                  writes <LINE>_full_universe*, so pass UNIV=${LINE}_full for those packs.
#   DEVICE=cpu     runs on CPU. mamba needs RIBO_MAMBA_CPU=1 as well (mamba_ssm's fast path asserts
#                  x.is_cuda; the pure-PyTorch reference twins are swapped in). CPU mamba is fine for
#                  small tx sets but did NOT finish a 41k-tx dump in 12 h -- prefer a GPU there.
UNIV=${UNIV:-$LINE}
DEVICE=${DEVICE:-cuda}
NVFLAG=""; [ "$DEVICE" = cuda ] && NVFLAG="--nv" 
case "$KOZAK" in
  heuristic) ORFTRACK=$NEW/data/packed_heldout_${HELDOUT}/orf_track_v2.npy ;;
  nokozak)   ORFTRACK=$NEW/data/packed_heldout_${HELDOUT}/orf_track_v2_nokozak.npy ;;
  *) echo "KOZAK must be heuristic|nokozak, got '$KOZAK'" >&2; exit 2 ;;
esac
# Name must match content: verify the track's own metadata agrees with KOZAK before predicting.
META=${ORFTRACK%.npy}_meta.json
if [[ -s "$META" ]]; then
  want=$([[ "$KOZAK" == nokozak ]] && echo none || echo heuristic)
  grep -q "\"kozak\": *\"$want\"" "$META" || {
    echo "ABORT: $META does not say kozak=$want (KOZAK=$KOZAK). The bare filename is NOT" >&2
    echo "       authoritative -- packed_union carries nokozak content under the bare name." >&2
    exit 1; }
fi
echo "MODEL_RUN=$RUN  SIF=$(basename "$SIF")  KOZAK=$KOZAK  track=$(basename "$ORFTRACK")"
echo "host=$(hostname) LINE=$LINE HELDOUT=$HELDOUT start=$(date -Iseconds)"
[[ -f "$ORFTRACK" ]] || { echo "missing orf_track $ORFTRACK" >&2; exit 1; }
[[ -f "$PILOT/${UNIV}_universe.fa" ]] || { echo "missing universe fasta $PILOT/${UNIV}_universe.fa" >&2; exit 1; }
# shellcheck disable=SC2086
singularity exec $NVFLAG --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env MPLCONFIGDIR=/tmp --env HF_HOME=/tmp \
  --env RIBO_MAMBA_CPU="${RIBO_MAMBA_CPU:-}" \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$PILOT/${UNIV}_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$RUN" \
    --heldout "$HELDOUT" --tx_list "$PILOT/${UNIV}_universe_tx.txt" \
    --fasta "$PILOT/${UNIV}_universe.fa" --out "$PILOT/$OUTSUB" --budget 16000 --device "$DEVICE" 
echo "done=$(date -Iseconds)"