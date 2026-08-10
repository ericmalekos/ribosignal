#!/usr/bin/env bash
#SBATCH --job-name=l3x3_dump
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/l3x3_dump_%A_%a.out
set -euo pipefail
# STAGE 4 of the mouse-liver 3x3: predict every cell with both released models, then run all three
# RiboCode drop-in density variants in the same job.
#
# GRID. 18 tasks = 2 models x 9 cells. Task index t: model = t/9, cell = t%9 (cell = ribo*3 + rna).
#
# WHY 18 AND NOT 6. What each variant actually depends on:
#   real            observed P-sites only            -> 3 distinct results (one per Ribo source)
#   pred_preddepth  model + sequence + RNA coverage  -> 6 distinct (2 models x 3 RNA)
#   pred_obsdepth   predicted shape x OBSERVED depth -> 18 distinct (2 x 3 RNA x 3 Ribo)
# pred_obsdepth is the arm that needs the full grid, so the full grid is what runs. The redundancy in
# the other two arms is not waste: the same Ribo source must yield byte-identical `real` calls no
# matter which RNA sits in the pack, and that is a free internal consistency check on the assembly.
#
# BOTH ARMS ALWAYS (feedback_two_arm_orf_calling): pred_obsdepth and pred_preddepth are both reported,
# never one alone. Poisson calibration is applied at scoring time, not here.
#
# PARTITION. attn runs on CPU, so submit tasks 0-8 to a CPU partition. mamba4 CANNOT: mamba_ssm
# dispatches to causal_conv1d_cuda and raises `Expected x.is_cuda() to be true` at the first block,
# with no CPU fallback in the installed build -- tasks 9-17 need a GPU. Guarded below so that failure
# is immediate and explained rather than 15 s into a job.
#   CPU:  sbatch --partition=medium --export=ALL,DEVICE=cpu --array=0-8%4  scripts/liver3x3_04_dump.sbatch
#   GPU:  sbatch --partition=gpu --gres=gpu:A5500:1 --array=9-17%3        scripts/liver3x3_04_dump.sbatch
# Do NOT queue a CPU and a GPU copy of the same task: they write the same paths, and two writers
# usually produce an output of exactly the right size that is wrong.
DEVICE=${DEVICE:-cuda}
case "$DEVICE" in
  cuda) NVFLAG=(--nv) ;;
  cpu)  NVFLAG=() ;;
  *) echo "bad DEVICE=$DEVICE (want cuda|cpu)" >&2; exit 1 ;;
esac
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
ANNOT=$NEW/data/mouse_ribocode_annot
L=$NEW/data/liver3x3

DS=(janich gse243134 wang)
MI=$((SLURM_ARRAY_TASK_ID / 9))
CELL=$((SLURM_ARRAY_TASK_ID % 9))
RB=${DS[$((CELL / 3))]}
RN=${DS[$((CELL % 3))]}
MODELS=(attn mamba4)
M=${MODELS[$MI]}
NAME=l3x3_ribo-${RB}_rna-${RN}

case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
          [[ "$DEVICE" == "cpu" ]] && { echo "ABORT: mamba4 cannot run with DEVICE=cpu. mamba_ssm \
dispatches to causal_conv1d_cuda ('Expected x.is_cuda() to be true'); the installed build has no CPU \
fallback. Submit tasks 9-17 to the gpu partition." >&2; exit 1; } ;;
esac

PACK=$NEW/data/packed_heldout_$NAME
ORFTRACK=$L/orf_track_union3_nokozak.npy
FASTA=$NEW/data/heldout_refs/mouse_liver_union3_universe.fa
TXLIST=$NEW/data/heldout_refs/mouse_liver_union3_tx.txt
OUT=$NEW/results/mouse_liver_3x3/${M}_ribo-${RB}_rna-${RN}
mkdir -p "$OUT"

for f in "$ORFTRACK" "$FASTA" "$TXLIST" "$NEW/$RUN/best.pt" "$PACK/tx_order.txt"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# Name must match content: the released models trained no-Kozak (feedback_kozak_never_set).
META=${ORFTRACK%.npy}_meta.json
if [[ -s "$META" ]]; then
  grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }
fi
NTX=$(wc -l < "$PACK/tx_order.txt")
echo "task=$SLURM_ARRAY_TASK_ID model=$M ribo=$RB rna=$RN tx=$NTX device=$DEVICE host=$(hostname) start=$(date -Iseconds)"

if [[ ! -s "$OUT/pred_profiles.npz" ]]; then
  singularity exec "${NVFLAG[@]}" --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
      --heldout "$NAME" --tx_list "$TXLIST" --fasta "$FASTA" \
      --out "$OUT" --budget 16000 --device "$DEVICE"
fi

# Gate on transcripts actually dumped, not on exit status: dump_pred_profiles DROPS transcripts absent
# from the one-hot FASTA without erroring, so a silent skip yields a biased subset rather than a crash
# (feedback_dump_silent_skip_onehot_fasta).
GOT=$(/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np;z=np.load('$OUT/pred_profiles.npz',allow_pickle=True);print(len(z['tx_ids']))")
echo "  dumped tx: $GOT of $NTX"
[[ "$GOT" -ge $((NTX * 95 / 100)) ]] || { echo "ABORT: only $GOT of $NTX dumped -- biased subset" >&2; exit 1; }

for V in real pred_obsdepth pred_preddepth; do
  [[ -s "$OUT/${V}_collapsed.txt" ]] && { echo "  skip $V (exists)"; continue; }
  "$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
  echo "  $V: $(($(wc -l < "$OUT/${V}_collapsed.txt") - 1)) calls"
done
echo "done=$(date -Iseconds) -> $OUT"
