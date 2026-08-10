#!/usr/bin/env bash
#SBATCH --job-name=rnaq_dump
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaq_dump_%A_%a.out
set -euo pipefail
# Task 68, stage 2: predict arms b and c with both released models, then run the 3 RiboCode drop-in
# density variants in the same job.
#
# Arm a (universe Janich-decon + Janich coverage) is already dumped and scored at
# results/liver_released/{attn,mamba4}_mouse_janich_liver_decon, on the SAME released checkpoints, so
# it is reused as-is rather than re-run. That is the whole comparison: the only thing that differs
# between a, b and c is the RNA-seq.
#
#   task 0  attn   arm b     1  mamba4 arm b
#   task 2  attn   arm c     3  mamba4 arm c
#
#   sbatch --array=0-3%2 scripts/rna_quality_factorial_dump.sbatch
#
# GPU jobs, but the drop-in tail is CPU. Kept in one job anyway: splitting it would need an afterok
# chain, and an afterok chain is what poisoned the original held-out run when one shard timed out.
#
# CPU FALLBACK. The dump itself is only 3-6 min on an A5500 for a pack this size, so when the gpu
# partition is saturated the wait dominates the work. Run the identical job on CPU with:
#
#   sbatch --partition=medium --gres=NONE --export=ALL,DEVICE=cpu --array=0-3%2 \
#       scripts/rna_quality_factorial_dump.sbatch
#
# Same script, same outputs -- only --device and the --nv flag change, so the two paths cannot drift.
# Do NOT have both a GPU and a CPU copy of the same array in the queue: they write the same paths, and
# the usual result of two writers is an output of exactly the right size that is wrong.
DEVICE=${DEVICE:-cuda}
case "$DEVICE" in
  cuda) NVFLAG=(--nv) ;;
  cpu)  NVFLAG=() ;;
  *) echo "bad DEVICE=$DEVICE (want cuda|cpu)" >&2; exit 1 ;;
esac
# ONLY the attn tasks (0, 2) can take DEVICE=cpu. mamba4 cannot: mamba_ssm dispatches to
# causal_conv1d_cuda, which raises `Expected x.is_cuda() to be true` at the first block -- there is no
# CPU fallback in the installed build. Fail here with that explanation rather than 15 s into a job.
# CPU array:  sbatch --partition=medium --gres=NONE --export=ALL,DEVICE=cpu --array=0,2 ...
# GPU array:  sbatch --array=1,3 ...
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"

MODELS=(attn mamba4 attn mamba4)
ARMS=(rnaq_b_uniJanich_cov34766 rnaq_b_uniJanich_cov34766 rnaq_c_uni34766_cov34766 rnaq_c_uni34766_cov34766)
# arm b reuses the Janich-decon universe FASTA verbatim (--ref-pack copied its tx_order), arm c has
# its own. Getting this wrong is the silent-skip failure mode: dump_pred_profiles drops transcripts
# absent from the one-hot FASTA WITHOUT erroring (feedback_dump_silent_skip_onehot_fasta), so the
# dumped tx count is gated against the pack below rather than trusted.
REFS=(mouse_janich_decon mouse_janich_decon mouse_prjeb34766 mouse_prjeb34766)

M=${MODELS[$SLURM_ARRAY_TASK_ID]}
DS=${ARMS[$SLURM_ARRAY_TASK_ID]}
REF=${REFS[$SLURM_ARRAY_TASK_ID]}

case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
          [[ "$DEVICE" == "cpu" ]] && { echo "ABORT: mamba4 cannot run with DEVICE=cpu. mamba_ssm \
dispatches to causal_conv1d_cuda, which raises 'Expected x.is_cuda() to be true' at the first block; \
the installed build has no CPU fallback. Submit the mamba4 tasks (1,3) to the gpu partition." >&2
            exit 1; } ;;
esac
PACK=$NEW/data/packed_heldout_$DS
ORFTRACK=$PACK/orf_track_v2_nokozak.npy
FASTA=$NEW/data/heldout_refs/${REF}_universe.fa
TXLIST=$NEW/data/heldout_refs/${REF}_universe_tx.txt
OUT=$NEW/results/rna_quality_factorial/${M}_${DS}
mkdir -p "$OUT"

for f in "$ORFTRACK" "$FASTA" "$TXLIST" "$NEW/$RUN/best.pt" "$PACK/tx_order.txt"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# Name must match content: the released models trained no-Kozak (feedback_kozak_never_set).
META=$(readlink -f "$ORFTRACK"); META=${META%.npy}_meta.json
if [[ -s "$META" ]]; then
  grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }
fi
NTX=$(wc -l < "$PACK/tx_order.txt")
echo "model=$M arm=$DS pack_tx=$NTX device=$DEVICE host=$(hostname) start=$(date -Iseconds)"

if [[ ! -s "$OUT/pred_profiles.npz" ]]; then
  singularity exec "${NVFLAG[@]}" --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
      --heldout "$DS" --tx_list "$TXLIST" --fasta "$FASTA" \
      --out "$OUT" --budget 16000 --device "$DEVICE"
fi

# Gate on tx actually dumped, not on exit status: a silent skip yields a biased subset, not an error.
GOT=$(/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np;z=np.load('$OUT/pred_profiles.npz',allow_pickle=True);print(len(z['tx_ids']))")
echo "dumped tx: $GOT of $NTX pack tx"
[[ "$GOT" -ge $((NTX * 95 / 100)) ]] || { echo "ABORT: only $GOT of $NTX dumped -- biased subset" >&2; exit 1; }

for V in real pred_obsdepth pred_preddepth; do
  [[ -s "$OUT/${V}_collapsed.txt" ]] && { echo "skip $V"; continue; }
  "$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
  echo "  $V: $(($(wc -l < "$OUT/${V}_collapsed.txt") - 1)) calls"
done
echo "done=$(date -Iseconds) -> $OUT"
