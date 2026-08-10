#!/usr/bin/env bash
#SBATCH --job-name=dump_liver_rel
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dump_liver_rel_%A_%a.out
set -euo pipefail
# Predict mouse liver with the RELEASED models. 2 models x 2 datasets = 4 dumps.
#
# WHY THIS IS NEEDED. Every mouse-liver model number in results/o2_liver came from
# `orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes` -- the PRE-UNION Fibroblast-universe
# checkpoint (33,918 test tx, Pearson 0.6230), not either released model (union, 70,883 tx, attn
# 0.6585 / mamba4 0.680). There was also no mamba4 liver arm at all, even though mamba4 is the
# primary release (D1b). The directory names encode the DATASET, not the checkpoint, so the staleness
# is invisible from the paths -- it is only in the npz `meta` field.
#
# DATASETS
#   mouse_janich_liver_decon  14,887 tx. The REBUILT Janich universe, from the decontaminated salmon
#                             quants (the shipped ones put 52-70% of TPM into Gm59647 + n-TKctt14,
#                             suppressing ~4,900 real transcripts below TPM >= 1). Janich's own RNA.
#   mouse_gse243134_liver     17,225 tx. GSE243134 Ribo-seq paired with GSE302188 poly(A) RNA
#                             (SRR34449454/471/553/569), because GSE243134's own totalRNA was rejected
#                             on chemistry (21.5% salmon, universe crushed to 4,451). This is a
#                             CROSS-STUDY RNA/Ribo pairing and must be labelled as such wherever the
#                             resulting numbers appear -- it is not a same-sample arm like the others.
#
#   mouse_wang_liver          Wang's own RNA (2 reps, 92.8% salmon, poly(A) -- the cleanest RNA of
#                             the three). Its PACK is the full 221,835-tx annotation, so the dump is
#                             restricted with --tx_list to the 13,041-tx EXPRESSED universe. That
#                             matches how Janich/GSE243134 are scoped (14,887 / 17,225 expressed) and
#                             avoids repeating the O2 mistake, where dumping over the whole
#                             annotation let the model call 122k ORFs against a ~17k reference and
#                             made precision an artifact of the transcript space.
#
# KOZAK: both packs carry orf_track_v2_nokozak.npy and the released models trained no-Kozak
# (feedback_kozak_never_set), so the track matches training. mamba4 needs the orthrus SIF for
# mamba_ssm; attn runs in either.
#
#   sbatch --array=0-5 scripts/dump_liver_released.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache

MODELS=(attn attn mamba4 mamba4 attn mamba4)
DSETS=(mouse_janich_liver_decon mouse_gse243134_liver mouse_janich_liver_decon mouse_gse243134_liver mouse_wang_liver mouse_wang_liver)
REFS=(mouse_janich_decon mouse_gse243134 mouse_janich_decon mouse_gse243134 mouse_wang_expressed mouse_wang_expressed)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
DS=${DSETS[$SLURM_ARRAY_TASK_ID]}
REF=${REFS[$SLURM_ARRAY_TASK_ID]}

case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
esac
PACK=$NEW/data/packed_heldout_$DS
ORFTRACK=$PACK/orf_track_v2_nokozak.npy
FASTA=$NEW/data/heldout_refs/${REF}_universe.fa
TXLIST=$NEW/data/heldout_refs/${REF}_universe_tx.txt
OUT=$NEW/results/liver_released/${M}_${DS}
mkdir -p "$OUT"

for f in "$ORFTRACK" "$FASTA" "$TXLIST" "$NEW/$RUN/best.pt"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# Name must match content: the released models trained no-Kozak, so refuse a heuristic track.
META=${ORFTRACK%.npy}_meta.json
if [[ -s "$META" ]]; then
  grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }
fi
echo "model=$M dataset=$DS pack_tx=$(wc -l < "$PACK/tx_order.txt") host=$(hostname) start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --heldout "$DS" --tx_list "$TXLIST" --fasta "$FASTA" \
    --out "$OUT" --budget 16000 --device cuda
echo "done=$(date -Iseconds) -> $OUT"
ls -la "$OUT" | tail -n +2
