#!/usr/bin/env bash
#SBATCH --job-name=orf_abl
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf_abl_%A_%a.out
set -euo pipefail
# Phase 3.2: zero one ORF-track channel at a time, re-predict, re-call, and measure the per-class
# F1 cost. Substrate is mouse GSE243134 (transcripts never seen in training).
#
# The five channels have clean, separable meanings (scripts/build_orf_track.py):
#   0 orf_f0     inside an ATG..in-frame-stop stretch, start pos %3 == 0
#   1 orf_f1     ... %3 == 1
#   2 orf_f2     ... %3 == 2
#   3 start_ext  graded start propensity: AUG + 9 near-cognates, weighted by initiation efficiency
#   4 is_stop    first nt of any stop codon
# so an ablation reads directly as "how much does the model depend on frame occupancy / on knowing
# where starts are / on knowing where stops are".
#
# A sixth arm zeroes channels 0-2 TOGETHER. The three frame channels are mutually redundant -- a
# position excluded from f0 and f1 is likely in f2 -- so zeroing one at a time can understate the
# value of frame information entirely. The joint arm is the one to read for "does the model need
# frame occupancy at all".
#
# WHY THIS DECIDES PHASE 5. Internal ORFs score F1 0.022 against a 0.539-0.620 between-experiment
# ceiling, yet aggregate saliency shows the model localises their start codons reasonably well
# (mamba4 median start_rank 77 on mouse) and attends strongly inside them. If ablating the start
# channel barely moves internal F1, the failure is not a missing sequence feature and ISM around
# start codons -- Phase 5.1 -- would be the wrong experiment.
#
#   sbatch --array=0-11%4 scripts/orf_channel_ablation.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
ANNOT=$NEW/data/mouse_ribocode_annot
L=$NEW/data/liver3x3
SRC=$L/orf_track_union3_nokozak.npy
ABL=$L/ablated_tracks
PACK=$NEW/data/packed_heldout_l3x3canon_ribo-gse243134_rna-gse243134
FASTA=$NEW/data/heldout_refs/mouse_liver_union3_universe.fa
TXLIST=$NEW/data/heldout_refs/mouse_liver_union3_tx.txt
NAME=l3x3canon_ribo-gse243134_rna-gse243134
mkdir -p "$ABL"

ARMS=(ch0 ch1 ch2 ch3 ch4 f012)
MODELS=(attn mamba4)
AI=$((SLURM_ARRAY_TASK_ID % 6))
MI=$((SLURM_ARRAY_TASK_ID / 6))
ARM=${ARMS[$AI]}
M=${MODELS[$MI]}

case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
esac

TRACK=$ABL/orf_track_${ARM}.npy
# Build the ablated track once; concurrent tasks share it, so write to a temp and rename atomically
# rather than letting two array tasks interleave writes into the same path.
if [[ ! -s "$TRACK" ]]; then
  "$PY" - <<PYEOF
import numpy as np, pathlib, os
src = np.load("$SRC", mmap_mode="r")
out = np.array(src)                       # (sum_L, 5) float16, ~563 MB
arm = "$ARM"
cols = [0, 1, 2] if arm == "f012" else [int(arm[-1])]
out[:, cols] = 0
tmp = pathlib.Path("$TRACK").with_suffix(".tmp$SLURM_ARRAY_TASK_ID.npy")
np.save(tmp, out)                          # np.save writes exactly this path (already ends .npy)
os.replace(tmp, "$TRACK")                  # atomic within one filesystem: no torn file for a reader
print(f"  built ablated track {arm}: zeroed columns {cols}, shape {out.shape}")
PYEOF
fi
[[ -s "$TRACK" ]] || { echo "ABORT: no ablated track at $TRACK" >&2; exit 1; }

OUT=$NEW/results/orf_channel_ablation_canon/${M}_${ARM}
mkdir -p "$OUT"
NTX=$(wc -l < "$PACK/tx_order.txt")
echo "task=$SLURM_ARRAY_TASK_ID model=$M arm=$ARM tx=$NTX host=$(hostname) start=$(date -Iseconds)"

if [[ ! -s "$OUT/pred_profiles.npz" ]]; then
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env RIBO_ORF_TRACK="$TRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
    -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
      --heldout "$NAME" --tx_list "$TXLIST" --fasta "$FASTA" \
      --out "$OUT" --budget 16000 --device cuda
fi
GOT=$("$PY" -c "
import numpy as np;z=np.load('$OUT/pred_profiles.npz',allow_pickle=True);print(len(z['tx_ids']))")
echo "  dumped tx: $GOT of $NTX"
[[ "$GOT" -ge $((NTX * 95 / 100)) ]] || { echo "ABORT: only $GOT of $NTX dumped" >&2; exit 1; }

for V in pred_obsdepth pred_preddepth; do
  [[ -s "$OUT/${V}_collapsed.txt" ]] && { echo "  skip $V"; continue; }
  "$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$OUT" --min_aa 5 --pval 0.05
  echo "  $V: $(($(wc -l < "$OUT/${V}_collapsed.txt") - 1)) calls"
done
echo "done=$(date -Iseconds) -> $OUT"
