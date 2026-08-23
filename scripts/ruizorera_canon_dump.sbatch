#!/usr/bin/env bash
#SBATCH --job-name=rz_dump
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rz_dump_%A_%a.out
set -euo pipefail
# Drop-in dump for the FINAL-RECIPE iPSC-CM pack. Task #92, the last off-recipe arm of P5.
#
#   task 0  attn     task 1  mamba4
#
# Same shape as scripts/ruizorera_released_dumps.sbatch, differing only in the pack: --heldout
# human_ruizorera_canon resolves to data/packed_heldout_human_ruizorera_canon, and the output goes to
# released_<M>/dropin_canonpack so the ORIGINAL off-recipe dropin is left intact for comparison.
# Keeping both is the point -- the whole question is how much the recipe moved this arm.
#
# THE ORF TRACK AND ONE-HOT FASTA ARE UNCHANGED, and must be: the rebuild reused the universe
# verbatim (asserted byte-identical tx_order/lengths/coverage before the track was copied), so the
# track is still row-aligned. A track that no longer matched would produce plausible numbers rather
# than an error, which is why the rebuild gated on it.
#
#   sbatch --array=0-1 scripts/ruizorera_canon_dump.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PACK=$NEW/data/packed_heldout_human_ruizorera_canon
FASTA=$NEW/data/union_universe.fa                          # one-hot input: chrM-free
TXLIST=$NEW/data/heldout_refs/ruizorera_nochrM_tx.txt      # 36,655 tx
ORFTRACK=$PACK/orf_track_v2_nokozak.npy
NTX=$(wc -l < "$TXLIST")

MODELS=(attn mamba4); M=${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}
case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
esac
for f in "$ORFTRACK" "$FASTA" "$TXLIST" "$NEW/$RUN/best.pt" "$PACK/tx_order.txt"; do
  [[ -e "$f" ]] || { echo "ABORT: missing $f" >&2; exit 1; }
done
META=${ORFTRACK%.npy}_meta.json
[[ -s "$META" ]] && { grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }; }

OUT=$NEW/results/heldout/human_ruizorera/released_${M}/dropin_canonpack
mkdir -p "$OUT"
echo "model=$M pack=$(basename $PACK) pack_tx=$NTX host=$(hostname) start=$(date -Iseconds)"

if [[ ! -s "$OUT/pred_profiles.npz" ]]; then
  N=${SLURM_CPUS_PER_TASK:-8}
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env OMP_NUM_THREADS="$N" --env MKL_NUM_THREADS="$N" \
    --env OPENBLAS_NUM_THREADS="$N" --env NUMEXPR_NUM_THREADS="$N" \
    --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
      --heldout human_ruizorera_canon --fasta "$FASTA" --tx_list "$TXLIST" \
      --out "$OUT" --budget 16000 --device cuda
fi

# Gate on transcripts dumped, not exit status: absent-from-FASTA tx are dropped SILENTLY, turning a
# bad pairing into a biased subset rather than a crash (feedback_dump_silent_skip_onehot_fasta).
GOT=$("$PY" -c "
import numpy as np;z=np.load('$OUT/pred_profiles.npz',allow_pickle=True);print(len(z['tx_ids']))")
echo "  dumped tx: $GOT of $NTX (skipped $((NTX-GOT)))"
[[ "$GOT" -ge $((NTX * 95 / 100)) ]] || { echo "ABORT: only $GOT of $NTX dumped -- biased subset" >&2; exit 1; }

# Both arms, then the metrics P5 reads.
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ANNOT=$NEW/data/human_ribocode_annot_primary
export PYTHONNOUSERSITE=1
for V in real pred_obsdepth pred_preddepth; do
  [[ -s "$OUT/${V}_collapsed.txt" ]] && { echo "  skip $V (exists)"; continue; }
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
  echo "  $V: $(($(wc -l < "$OUT/${V}_collapsed.txt") - 1)) calls"
done
"$PY" "$NEW/scripts/compare_dropin_calls.py" --dropin_dir "$OUT" \
  --profiles "$OUT/pred_profiles.npz" --out "$OUT" 2>&1 | tail -4
echo "done=$(date -Iseconds) -> $OUT"
