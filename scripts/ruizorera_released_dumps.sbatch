#!/usr/bin/env bash
#SBATCH --job-name=ruiz_released
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/ruiz_released_%A_%a.out
set -euo pipefail
# Re-dump Ruiz-Orera on the RELEASED models. This is the ONLY genuinely independent human dataset in
# the project -- Hepatocytes is the same Chothani panel with a tissue withheld -- so the cross-study
# human claim rests entirely on it, and its only existing drop-in is on
# `orf_v2_attn_onehot_holdout_Hepatocytes`: pre-nokozak, pre-mm1, pre-union.
#
# The pack ships NO orf_track (it was built reusing another pack's universe), so the track is built
# here, row-aligned to this pack, with --kozak none to match how the released models trained
# (feedback_kozak_never_set).
#
# SILENT-SKIP GUARD. union_universe.fa covers 36,655 of the pack's 36,668 tx (99.96%); 13 have no
# one-hot entry. dump_pred_profiles SKIPS such transcripts without erroring, and a cross-universe eval
# can therefore collapse onto a biased subset with no warning
# (feedback_dump_silent_skip_onehot_fasta). The dumped tx count is checked against a floor below.
#
#   sbatch --array=0-1 scripts/ruizorera_released_dumps.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PACK=$NEW/data/packed_heldout_human_ruizorera
# The pack still contains the 13 chrM protein-coding mRNAs (MT-ND1/2/3/4/4L/5/6, MT-CO1/2/3, MT-CYB,
# MT-ATP6/8); it predates prepare_pack.py and has no provenance.json. union_universe.fa correctly
# EXCLUDES them per feedback_riboseq_exclude_mito_genes, so build_orf_track -- which hard-fails on any
# packed tx missing from the FASTA -- aborted in 3 s. That guard is right; the pack is what is stale.
#
# Rather than rebuild a pack whose inputs are no longer recorded: build the track over the FULL pack
# using a fasta augmented with the 13 chrM sequences (so the row-aligned track is complete), then
# EXCLUDE chrM from what is actually predicted and scored via --tx_list. The chrM track rows exist and
# are never used. Net effect matches the rule: no chrM in the eval.
TRACK_FASTA=$NEW/data/heldout_refs/ruizorera_track_universe.fa   # union + 13 chrM
FASTA=$NEW/data/union_universe.fa                                # one-hot input: chrM-free
TXLIST=$NEW/data/heldout_refs/ruizorera_nochrM_tx.txt            # 36,655 tx
ORFTRACK=$PACK/orf_track_v2_nokozak.npy
NTX=$(wc -l < "$TXLIST")

MODELS=(attn mamba4)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
esac

# Build the pack-aligned ORF track once; guard so the two array tasks cannot race on one path.
if [[ ! -s "$ORFTRACK" ]]; then
  LOCK=$PACK/.building_track
  if mkdir "$LOCK" 2>/dev/null; then
    # Release the lock even if the build fails. The first attempt died under `set -e` BEFORE its
    # rmdir, so the sibling array task spun on a lock that could never clear until walltime.
    trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
    echo "[$(date +%T)] building ORF track over the full pack (kozak none)"
    "$PY" "$NEW/scripts/build_orf_track.py" --pack "$PACK" --fasta "$TRACK_FASTA" \
        --kozak none --mode ext --out "$ORFTRACK"
    rmdir "$LOCK"; trap - EXIT
  else
    while [[ -d "$LOCK" ]]; do sleep 10; done
  fi
fi
[[ -s "$ORFTRACK" ]] || { echo "no ORF track at $ORFTRACK" >&2; exit 1; }

OUT=$NEW/results/heldout/human_ruizorera/released_${M}/dropin
mkdir -p "$OUT"
echo "model=$M pack_tx=$NTX host=$(hostname) start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --heldout human_ruizorera --fasta "$FASTA" --tx_list "$TXLIST" \
    --out "$OUT" --budget 16000 --device cuda

# Gate on how many transcripts actually made it through, not on exit status.
GOT=$("$PY" -c "
import numpy as np;z=np.load('$OUT/pred_profiles.npz',allow_pickle=True);print(len(z['tx_ids']))")
echo "dumped tx: $GOT of $NTX pack tx (skipped $((NTX-GOT)))"
[[ "$GOT" -ge $((NTX * 95 / 100)) ]] || { echo "ABORT: only $GOT of $NTX dumped -- biased subset" >&2; exit 1; }
echo "done=$(date -Iseconds) -> $OUT"
