#!/usr/bin/env bash
#SBATCH --job-name=train_canon
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/train_canon_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/train_canon_%j.err
set -euo pipefail
# Retrain the deployed architectures on CANONICAL training packs (data/packed_canon*).
#
# Everything writes under results/loto_canon/ so nothing here can be confused with, or overwrite,
# the deployed models in results/loto/.
#
# THE ONLY INTENDED DIFFERENCE FROM THE DEPLOYED RUN IS THE TRAINING DATA. Every hyperparameter
# below was read back out of the deployed runs' args.json rather than retyped from a script, because
# the comparison is only interpretable if the data is the sole variable. attn and mamba4 differ from
# each other in exactly two args (mixer, n_attn_layers); both match their deployed counterparts on
# all 30+ others (seed 0, val_fold 0, budget 16000, max_tx_per_tissue 6000, min_train_signal 50,
# lr 3e-4, epochs 30, patience 6, warmup 2, channels 256, n_blocks 10, dropout 0.1, count_weight 0.1,
# peakiness_weight 0.0, weight_decay 0.01, cov_norm global_mean, eval_cap 2500).
#
# PACKS: RIBO_PACK_SUFFIX=canon makes dataset.tissue_pack_dir resolve to data/packed_canon (Fibroblast)
# and data/packed_canon_<T> (others). No code change; that hook already existed for the union family.
#
# ORF TRACK: reused from data/packed_union/orf_track_v2.npy. It is a function of sequence and
# annotation, not of reads, and the final-recipe packs share the universe byte-for-byte (asserted in
# build_canonical_chothani_packs.sbatch). Its meta records kozak=none, which is required
# (feedback_kozak_never_set); the file named *_nokozak.npy in that directory is a symlink to this
# same array, so both names denote the no-Kozak track.
#
#   MODEL=attn   sbatch scripts/train_canon.sbatch
#   MODEL=mamba4 sbatch scripts/train_canon.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
MODEL=${MODEL:?set MODEL=attn or MODEL=mamba4}
# The one-hot FASTA is NOT optional here. dataset.py defaults RIBO_ONEHOT_FASTA to
# data/fibroblast_universe.fa (36,668 tx), but the final-recipe packs are on the UNION universe
# (84,472 tx). Leaving it unset cost a run: KeyError 'ENST00000869382.1' 68 s in, on the first
# transcript outside the fibroblast subset. Same value train_loto_union.sbatch uses.
ONEHOT_FASTA=$NEW/data/union_universe.fa
# ORFTRACK/TAG let an architecture arm swap the ORF track WITHOUT editing this script. Default is the
# full union track; TAG must change with it so the run cannot land in the same output directory as
# the full-track run and be mistaken for it.
#   frame ablation:  ORFTRACK=$NEW/data/packed_union/orf_track_v2_nof012.npy TAG=nof012
ORFTRACK=${ORFTRACK:-$NEW/data/packed_union/orf_track_v2.npy}
TAG=${TAG:-}

# SIF FOLLOWS THE MIXER, not the project's default image: mamba_ssm ships only in the orthrus image.
# Hardcoding rinalmo for both cost a run (ModuleNotFoundError: mamba_ssm, 7 s in).
case "$MODEL" in
  attn)   MIXER=transformer; NATTN=2; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) MIXER=mamba;       NATTN=4; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  *) echo "bad MODEL=$MODEL (want attn|mamba4)" >&2; exit 1 ;;
esac
OUT=$NEW/results/loto_canon/orf_v2_${MODEL}_onehot_canon_noBrain_nokozak${TAG:+_$TAG}_holdout_Hepatocytes
# A non-default track with no TAG would overwrite the full-track run's directory. Refuse.
if [[ "$ORFTRACK" != "$NEW/data/packed_union/orf_track_v2.npy" && -z "$TAG" ]]; then
  echo "ABORT: non-default ORFTRACK requires TAG so the run gets its own output dir" >&2; exit 1
fi
mkdir -p "$NEW/logs/train" "$OUT"

# Fail before burning a GPU slot rather than 15 s into the run.
META=${ORFTRACK%.npy}_meta.json
[[ -s "$META" ]] && { grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }; }
for T in Fibroblast VSMC ES Fat HA_EC HCAEC HUVEC Hepatocytes; do
  d=$NEW/data/packed_canon; [[ "$T" != "Fibroblast" ]] && d=$NEW/data/packed_canon_$T
  [[ -s "$d/target_counts.npy" && -s "$d/pack_meta.tsv" ]] || { echo "ABORT: final-recipe pack missing: $d" >&2; exit 1; }
done
# PREFLIGHT: the two failures above were both silent until the dataloader hit a missing key or the
# model built a mamba block. Check the cheap invariants here instead.
[[ -s "$ONEHOT_FASTA" ]] || { echo "ABORT: missing one-hot FASTA $ONEHOT_FASTA" >&2; exit 1; }
[[ -s "$SIF" ]] || { echo "ABORT: missing SIF $SIF" >&2; exit 1; }
"$PY" - <<PYEOF || exit 1
import sys, pathlib
fa = pathlib.Path("$ONEHOT_FASTA")
ids = {l[1:].split()[0].split("|")[0] for l in open(fa) if l.startswith(">")}
order = [l.strip() for l in open("$NEW/data/packed_canon/tx_order.txt") if l.strip()]
miss = [t for t in order if t not in ids]
if miss:
    sys.exit(f"ABORT: one-hot FASTA {fa.name} misses {len(miss):,} of {len(order):,} universe tx "
             f"(e.g. {miss[:3]}). The final-recipe packs are on the UNION universe.")
print(f"  preflight OK: {fa.name} covers all {len(order):,} universe tx")
PYEOF

echo "model=$MODEL mixer=$MIXER n_attn_layers=$NATTN packs=canon holdout=Hepatocytes host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_PACK_SUFFIX=canon \
  --env RIBO_ONEHOT_FASTA="$ONEHOT_FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes \
    --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC \
    --emb_backend onehot --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    --use_orf_track --mixer "$MIXER" --n_attn_layers "$NATTN" \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --n_heads 8 --lr 3e-4 --weight_decay 0.01 \
    --count_weight 0.1 --peakiness_weight 0.0 \
    --budget 16000 --eval_cap 2500 --seed 0 \
    --num_workers 8 --device cuda
echo "done=$(date -Iseconds) -> $OUT"
