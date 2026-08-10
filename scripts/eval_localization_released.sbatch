#!/usr/bin/env bash
#SBATCH --job-name=loc_released
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/loc_released_%A_%a.out
set -euo pipefail
# Localization eval on the RELEASED models, so main Figure 1 stops being built from a stale checkpoint.
#
# WHY. A blast-radius scan of every figure generator (not just the ones under suspicion -- see
# feedback_blast_radius_scan_by_artifact) found that FIVE of the six Fig 1 panels read
# `orf_v2_attn_onehot_holdout_Hepatocytes`, the pre-nokozak / pre-mm1 / pre-union checkpoint:
#   A1 localization ceiling, A2 drop-in, B5 expression independence, C10 saliency, prediction_examples.
# Only arch_attn was on a released model. A2 and prediction_examples only needed repointing (their
# released artifacts already exist). A1 and B5 both read `localization_metrics.json`, which the
# released runs do NOT have -- that is what this job produces.
#
# A3 (LOTO spread) reads results/loto_9fold/onehot_orf_v2_attn_holdout_*, for which there is no union
# equivalent; re-running it means a 9-fold retrain, so it is left labelled rather than silently mixed
# with released-model panels.
#
#   task 0  attn   Hepatocytes    1  mamba4 Hepatocytes
#
#   sbatch --array=0-1 scripts/eval_localization_released.sbatch
#
# CPU FALLBACK, attn ONLY. When the gpu partition is saturated the queue wait dominates the work:
#   sbatch --partition=medium --gres=NONE --export=ALL,DEVICE=cpu --array=0 \
#       scripts/eval_localization_released.sbatch
# mamba4 (task 1) CANNOT: mamba_ssm dispatches to causal_conv1d_cuda and raises
# `Expected x.is_cuda() to be true` at the first block, with no CPU fallback in the installed build.
# The script aborts on that combination rather than failing 15 s in. A1 and B5 read only the ATTN
# localization_metrics.json, so the CPU task alone unblocks both figures.
# Never queue a GPU and a CPU copy of the same task: they write the same path.
DEVICE=${DEVICE:-cuda}
case "$DEVICE" in
  cuda) NVFLAG=(--nv) ;;
  cpu)  NVFLAG=() ;;
  *) echo "bad DEVICE=$DEVICE (want cuda|cpu)" >&2; exit 1 ;;
esac
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
# EXACTLY what train_loto_union.sbatch passed. All three must move together -- pack suffix, ORF track,
# one-hot FASTA -- because they define the same universe.
#
# The first version of this script used data/packed/orf_track_v2_nokozak.npy and its kozak=none check
# PASSED, because that file IS kozak=none. It is also the FIBROBLAST-universe track (980 MB) rather
# than the union one (2,126 MB), so the assertion tested the right property on the wrong axis. Both
# identity dimensions are now checked: kozak, and byte size against the training track.
ORFTRACK=$NEW/data/packed_union/orf_track_v2.npy
FASTA=$NEW/data/union_universe.fa
PACK_SUFFIX=union

MODELS=(attn mamba4)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
          [[ "$DEVICE" == "cpu" ]] && { echo "ABORT: mamba4 cannot run with DEVICE=cpu (mamba_ssm -> \
causal_conv1d_cuda has no CPU fallback). Submit task 1 to the gpu partition." >&2; exit 1; } ;;
esac
for f in "$ORFTRACK" "$FASTA" "$NEW/$RUN/best.pt"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# Identity check 1: no-Kozak, matching how the released models trained (feedback_kozak_never_set).
META=${ORFTRACK%.npy}_meta.json
[[ -s "$META" ]] && { grep -q '"kozak": *"none"' "$META" || \
  { echo "ABORT: $META is not kozak=none" >&2; exit 1; }; }
# Identity check 2: the UNIVERSE. kozak=none is true of both the union track and the Fibroblast one,
# so it cannot distinguish them; their sizes can (2,126,471,848 vs 979,634,968 bytes).
WANT_BYTES=2126471848
GOT_BYTES=$(stat -c%s "$ORFTRACK")
[[ "$GOT_BYTES" -eq "$WANT_BYTES" ]] || { echo "ABORT: $ORFTRACK is $GOT_BYTES bytes, expected \
$WANT_BYTES (the union track). A same-kozak track from a different universe would pass the check \
above and silently mis-align every ORF." >&2; exit 1; }
echo "model=$M run=$RUN device=$DEVICE track=$(basename "$ORFTRACK") ($GOT_BYTES B) host=$(hostname) start=$(date -Iseconds)"

singularity exec "${NVFLAG[@]}" --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_PACK_SUFFIX="$PACK_SUFFIX" \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" \
    --run "$NEW/$RUN" --min_orf_nt 30 --eval_cap 0 --device "$DEVICE"

OUT=$NEW/$RUN/localization_metrics.json
[[ -s "$OUT" ]] || { echo "ABORT: no $OUT written" >&2; exit 1; }
echo "done=$(date -Iseconds) -> $OUT"
