#!/usr/bin/env bash
#SBATCH --job-name=hs_dump
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hs_dump_%A_%a.out
set -euo pipefail
# Predict each human dataset with both released models, then run all three RiboCode drop-in variants.
#
#   task 0-2  attn    x {gse208041, gse39561, cart}
#   task 3-5  mamba4  x {gse208041, gse39561, cart}
#
# GPU FOR BOTH ARMS. mamba4 has no choice (mamba_ssm dispatches to causal_conv1d_cuda). attn *can* run
# on CPU, and doing so on the mouse 3x3 was a mistake worth not repeating: mamba4 averaged 9 min per
# cell on an A5500 while attn on CPU ran 2h50m without finishing one. The CPU path is a fallback for a
# saturated GPU partition, not a default.
#
# Each dataset has its OWN universe and therefore its OWN ORF track, unlike the mouse 3x3 where all
# nine cells shared one. RIBO_ORF_TRACK and RIBO_ONEHOT_FASTA are set per task and must move together
# with the pack -- they define one universe between them.
#
# PACKSUF selects the pack VARIANT without forking this script. Empty (default) = the original
# packs; `_canon` = the realigned with the final recipe ones (data/packed_heldout_<name>_canon), which write
# to results/human_orf_calls_canon/. The universe, ORF track and one-hot FASTA are unchanged by a
# final-recipe rebuild (rebuild_heldout_canon.sbatch asserts tx_order and lengths identical and copies
# the track), so only the pack and output paths move. Forking this file into a `human_02_dump_canon`
# is exactly the mistake train_canon.sbatch made -- see docs/CANONICAL_ENTRYPOINTS.md.
#
#   sbatch --array=0-5%3 scripts/human_02_dump.sbatch
#   PACKSUF=_canon sbatch --array=2,5 scripts/human_02_dump.sbatch      # CAR-T only, both models
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
ANNOT=$NEW/data/human_ribocode_annot_primary
REF=$NEW/data/heldout_refs
DEVICE=${DEVICE:-cuda}
case "$DEVICE" in cuda) NVFLAG=(--nv);; cpu) NVFLAG=();; *) echo "bad DEVICE=$DEVICE" >&2; exit 1;; esac

DS=(gse208041 gse39561 cart)
UNIS=(human_thp1_gse208041 human_thp1_gse39561 human_cart_gse304796)
MI=$((SLURM_ARRAY_TASK_ID / 3))
CI=$((SLURM_ARRAY_TASK_ID % 3))
D=${DS[$CI]}
U=${UNIS[$CI]}
MODELS=(attn mamba4)
M=${MODELS[$MI]}

case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
          [[ "$DEVICE" == "cpu" ]] && { echo "ABORT: mamba4 cannot run on CPU (causal_conv1d_cuda)." >&2
                                        exit 1; } ;;
esac

PACKSUF=${PACKSUF:-}
NAME=human_$D
PACK=$NEW/data/packed_heldout_${NAME}${PACKSUF}
ORFTRACK=$PACK/orf_track_v2_nokozak.npy
FASTA=$REF/${U}_universe.fa
TXLIST=$REF/${U}_tx.txt
OUT=$NEW/results/human_orf_calls${PACKSUF}/${M}_${D}
# The pack variant must not silently score against the wrong universe: dump_pred_profiles resolves
# the pack from --heldout, so that name has to carry the suffix too.
NAME=${NAME}${PACKSUF}
mkdir -p "$OUT"

for f in "$ORFTRACK" "$FASTA" "$TXLIST" "$NEW/$RUN/best.pt" "$PACK/tx_order.txt"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
META=${ORFTRACK%.npy}_meta.json
if [[ -s "$META" ]]; then
  grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }
fi
NTX=$(wc -l < "$PACK/tx_order.txt")
echo "task=$SLURM_ARRAY_TASK_ID model=$M dataset=$D pack=${PACKSUF:-<default>} tx=$NTX device=$DEVICE threads=${SLURM_CPUS_PER_TASK:-8} host=$(hostname) start=$(date -Iseconds)"

if [[ ! -s "$OUT/pred_profiles.npz" ]]; then
  # THREAD COUNT MUST BE PASSED EXPLICITLY. --cleanenv drops the host environment, so torch inside
  # the container sees the NODE's core count (256 on a phoenix GPU node), not this job's cgroup
  # allocation, and spawns that many OMP threads. Oversubscribing ~4x is typically SLOWER than
  # running single-threaded, and is the likely reason an earlier CPU attempt ran 2h50m without
  # finishing a comparable dataset. Harmless on the GPU path, decisive on the CPU path.
  NTHREADS=${SLURM_CPUS_PER_TASK:-8}
  singularity exec "${NVFLAG[@]}" --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env OMP_NUM_THREADS="$NTHREADS" --env MKL_NUM_THREADS="$NTHREADS" \
    --env OPENBLAS_NUM_THREADS="$NTHREADS" --env NUMEXPR_NUM_THREADS="$NTHREADS" \
    --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
      --heldout "$NAME" --tx_list "$TXLIST" --fasta "$FASTA" \
      --out "$OUT" --budget 16000 --device "$DEVICE"
fi

# Gate on transcripts dumped, not exit status: absent-from-FASTA transcripts are DROPPED silently
# (feedback_dump_silent_skip_onehot_fasta), turning a bad universe into a biased subset, not a crash.
GOT=$(/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np;z=np.load('$OUT/pred_profiles.npz',allow_pickle=True);print(len(z['tx_ids']))")
echo "  dumped tx: $GOT of $NTX"
[[ "$GOT" -ge $((NTX * 95 / 100)) ]] || { echo "ABORT: only $GOT of $NTX dumped -- biased subset" >&2; exit 1; }

# Both prediction arms, always (feedback_two_arm_orf_calling) -- EXCEPT where the data cannot support
# one. GSE39561 has no detectable 3-nt periodicity (metaplots selects zero read lengths), so its pack
# is coverage-only: there are no observed P-sites, hence no `real` arm and no `pred_obsdepth`, which
# scales predicted shape by observed depth. Only the standalone arm is defined. Running the other two
# would produce empty or meaningless output rather than an error, so they are skipped explicitly.
ARMS=(real pred_obsdepth pred_preddepth)
if ! /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np,sys; sys.exit(0 if np.load('$PACK/target_counts.npy').sum()>0 else 1)" 2>/dev/null; then
  echo "  pack has ZERO P-sites (coverage-only): standalone arm only"
  ARMS=(pred_preddepth)
fi
for V in "${ARMS[@]}"; do
  [[ -s "$OUT/${V}_collapsed.txt" ]] && { echo "  skip $V (exists)"; continue; }
  "$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
  echo "  $V: $(($(wc -l < "$OUT/${V}_collapsed.txt") - 1)) calls"
done
echo "done=$(date -Iseconds) -> $OUT"
