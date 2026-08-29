#!/usr/bin/env bash
#SBATCH --job-name=xsp_pred
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/pred_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/pred_%A_%a.err
#
#   sbatch --export=ALL,DS=xsp_yeast,ARCH=attn scripts/xspecies/predict_xspecies.sbatch
#
# Runs the released no-Kozak union checkpoint over one cross-species pack.
#
# THE ORF TRACK IS PASSED EXPLICITLY. Each pack has its own universe, so tx_order differs and a
# track built for another pack is silently misaligned -- there is no shape check that catches it.
# RIBO_ORF_TRACK points at this pack's own orf_track_v2_nokozak.npy.
#
# --kozak none is not optional: the released models trained on the no-Kozak track, and feeding
# the Kozak-heuristic track is a train/inference mismatch. The _nokozak filename is the guard.
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/riboseq-model.sif
DS=${DS:?set DS=<xsp_...>}
ARCH=${ARCH:?set ARCH=attn|mamba4}
UNI=${UNI:?set UNI=<ribo dataset label for the universe fasta>}
PACK=$NEW/data/packed_heldout_$DS
RUN=$NEW/results/loto/orf_v2_${ARCH}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
FASTA=$NEW/data/xspecies_refs/${UNI}_universe.fa
ORF=$PACK/orf_track_v2_nokozak.npy
OUT=$RUN/dropin_$DS

for f in "$PACK/tx_order.txt" "$FASTA" "$ORF" "$RUN"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# The ORF track's first axis must match the pack's tx count, or every downstream index is off.
NTX=$(wc -l < "$PACK/tx_order.txt")
echo "dataset=$DS arch=$ARCH tx=$NTX device=${DEV:-cpu} host=$(hostname) start=$(date -Iseconds)"

GPUFLAG=""; [[ "${DEV:-cpu}" == "cuda" ]] && GPUFLAG="--nv"
# mamba_ssm's fast path asserts x.is_cuda; RIBO_MAMBA_CPU=1 swaps in the pure-PyTorch twins.
MAMBACPU=""
[[ "$ARCH" == "mamba4" && "${DEV:-cpu}" == "cpu" ]] && MAMBACPU="--env RIBO_MAMBA_CPU=1"

singularity exec $GPUFLAG --no-home --cleanenv \
  --env PYTHONNOUSERSITE=1 --env MPLCONFIGDIR=/tmp --env HF_HOME=/tmp \
  --env RIBO_ORF_TRACK="$ORF" --env RIBO_ONEHOT_FASTA="$FASTA" $MAMBACPU \
  -B /private/groups "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" \
    --run "$RUN" --heldout "$DS" --out "$OUT" --fasta "$FASTA" \
    --device "${DEV:-cpu}" ${NSHARDS:+--nshards $NSHARDS} ${SHARD:+--shard $SHARD}

echo "=== $DS/$ARCH DONE=$(date -Iseconds)"
ls -lh "$OUT" 2>/dev/null | tail -3
