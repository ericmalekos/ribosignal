#!/usr/bin/env bash
#SBATCH --job-name=eval_canon
#SBATCH --partition=gpu
# ANY GPU, not gpu:A5500:1. This is a ~5M-parameter model; pinning the A5500s (copied from
# eval_union.sbatch) left the job queued behind a full A5500 pool while 5 A100s sat idle on
# phoenix-00. Override per-submission with `sbatch --gres=gpu:A5500:1` if a specific card is needed.
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval_canon_%A_%a.out
set -euo pipefail
# Task #86: does training on CANONICAL packs improve non-canonical ORF CALLING?
#
# THE QUESTION PROFILE METRICS CANNOT ANSWER. The final-recipe retrains exist and their profile metrics
# are a wash (attn pearson 0.659 canon vs the deployed run's comparable value), but the standing rule
# is that models are judged on how ORF CALLS change, not on profile correlation
# (feedback_orf_call_eval_standing). None of the three final-recipe retrains has ever been ORF-call
# evaluated -- none has a dropin/ directory. This produces them.
#
# THE COMPARISON IS CONTROLLED BY EVALUATING EVERYTHING ON ONE SUBSTRATE.
# Every task here runs with RIBO_PACK_SUFFIX=canon, so all five arms are scored against the SAME
# canonical Hepatocytes reference and consume the same canonical RNA coverage. The ONLY thing that
# differs between the deployed arms and the final-recipe arms is what the model was TRAINED on.
#
# Scoring the deployed model against the mm1 reference and the canonical model against the canonical
# reference would compare two different reference sets, and a difference could then be reference
# quality rather than model quality. That is the trap this layout exists to avoid, and it is the same
# rule as feedback_orf_call_transcript_space: state the space, hold it fixed.
#
#   task 0  CANON     attn            (loto_canon)
#   task 1  CANON     attn nof012     (frame-channel ablation, its own ORF track)
#   task 2  CANON     mamba4          (loto_canon)
#   task 3  DEPLOYED  attn      -> results/loto/<run>/dropin_canonpack   (NOT its mm1 dropin)
#   task 4  DEPLOYED  mamba4    -> results/loto/<run>/dropin_canonpack
#
# THE `real` ARM IS A BUILT-IN CONTROL. It depends only on the pack, so all five tasks must produce
# an IDENTICAL real_collapsed.txt. compare_canon_orfcalls.py asserts that; if they diverge, something
# about the pack resolution differs between arms and no other number here can be trusted.
#
# ORF TRACK MUST MATCH TRAINING. nof012 was trained on orf_track_v2_nof012.npy and must be dumped
# with it. orf_track_v2_nokozak.npy is a SYMLINK to orf_track_v2.npy (verified: identical md5 and
# size), so the deployed and canonical full-track arms share one file.
#
#   sbatch --array=0-4%3 scripts/eval_canon_orfcalls.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
# The live RiboCode annotation. eval_union.sbatch still points at
# expression_context_human/data/ribocode_annot, which NO LONGER EXISTS -- that script would fail
# today. GENCODE v49 primary, matching data/union_universe.fa (74,704 of 84,472 universe tx carry a
# CDS; the remainder are non-coding, as expected).
ANNOT=$NEW/data/human_ribocode_annot_primary
FASTA=$NEW/data/union_universe.fa
TRACK_FULL=$NEW/data/packed_union/orf_track_v2.npy
TRACK_NOF=$NEW/data/packed_union/orf_track_v2_nof012.npy
RIN=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
ORT=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img

T=${SLURM_ARRAY_TASK_ID:-0}
case "$T" in
  0) RUN=results/loto_canon/orf_v2_attn_onehot_canon_noBrain_nokozak_holdout_Hepatocytes
     TRACK=$TRACK_FULL; SIF=$RIN; OUT=""; LBL="canon/attn" ;;
  1) RUN=results/loto_canon/orf_v2_attn_onehot_canon_noBrain_nokozak_nof012_holdout_Hepatocytes
     TRACK=$TRACK_NOF;  SIF=$RIN; OUT=""; LBL="canon/attn_nof012" ;;
  2) RUN=results/loto_canon/orf_v2_mamba4_onehot_canon_noBrain_nokozak_holdout_Hepatocytes
     TRACK=$TRACK_FULL; SIF=$ORT; OUT=""; LBL="canon/mamba4" ;;
  3) RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
     TRACK=$TRACK_FULL; SIF=$RIN; OUT=dropin_canonpack; LBL="deployed/attn" ;;
  4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
     TRACK=$TRACK_FULL; SIF=$ORT; OUT=dropin_canonpack; LBL="deployed/mamba4" ;;
  *) echo "bad task $T" >&2; exit 1 ;;
esac
D=$NEW/$RUN/${OUT:-dropin}
mkdir -p "$D" "$NEW/logs"

for f in "$NEW/$RUN/best.pt" "$TRACK" "$FASTA" "$ANNOT/transcripts_cds.txt" \
         "$NEW/data/packed_canon/tx_order.txt"; do
  [[ -e "$f" ]] || { echo "ABORT: missing $f" >&2; exit 1; }
done
# Refuse to overwrite a deployed run's ORIGINAL mm1 dropin. Tasks 3/4 must land in dropin_canonpack.
if [[ "$RUN" == results/loto/* && "${OUT}" != "dropin_canonpack" ]]; then
  echo "ABORT: deployed run must write to dropin_canonpack, not $OUT" >&2; exit 1
fi
echo "task=$T $LBL run=$RUN out=$D track=$(basename "$TRACK") host=$(hostname) start=$(date -Iseconds)"

if [[ ! -s "$D/pred_profiles.npz" ]]; then
  N=${SLURM_CPUS_PER_TASK:-8}
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env OMP_NUM_THREADS="$N" --env MKL_NUM_THREADS="$N" \
    --env OPENBLAS_NUM_THREADS="$N" --env NUMEXPR_NUM_THREADS="$N" \
    --env RIBO_PACK_SUFFIX=canon \
    --env RIBO_ORF_TRACK="$TRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --out "$D" \
      --fasta "$FASTA" --device cuda --budget 16000
fi

# A heavily-skipped dump means a FASTA/pack mismatch and silently evaluates a biased subset
# (feedback_dump_silent_skip_onehot_fasta). Fail before any call is made.
"$PY" -c "
import numpy as np,sys
d=np.load('$D/pred_profiles.npz',allow_pickle=True)
nt=len(d['tx_ids']); ns=int(str(d['meta'][3])); frac=ns/(nt+ns) if (nt+ns) else 0
print(f'  dump: {nt:,} tx, {ns:,} skipped ({frac:.1%})')
sys.exit(1 if frac>0.05 else 0)" || { echo "ABORT: >5% skipped -- FASTA/pack mismatch" >&2; exit 1; }

export PYTHONNOUSERSITE=1
for V in real pred_obsdepth pred_preddepth; do
  [[ -s "$D/${V}_collapsed.txt" ]] && { echo "  skip $V (exists)"; continue; }
  echo "  === dropin $V ==="
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$D/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$D" --min_aa 5 --pval 0.05
  echo "  $V: $(($(wc -l < "$D/${V}_collapsed.txt") - 1)) calls"
done
echo "done=$(date -Iseconds) -> $D"
