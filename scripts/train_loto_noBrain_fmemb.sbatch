#!/usr/bin/env bash
#SBATCH --job-name=loto_fmemb
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_fmemb_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_fmemb_%j.err
set -euo pipefail
# FM-embedding variants of the current FINAL model (mm1 coverage + no-Kozak track, noBrain 7 tissues, holdout
# Hepatocytes) -- ONLY the input embedding differs from the one-hot deployed model, for a controlled comparison:
#   BACKEND=rinalmo  (1280-d)   BACKEND=orthrus (512-d)   BACKEND=hydrarna (1024-d)
#
# FM_TO_MIXER=1 adds the direct embedding->mixer path (dedicated projection + zero-init gate) on top of
# the normal in_proj->CNN route, and writes to a _fmmix run dir. Everything else is byte-identical to the
# baseline arm above, so the two are a controlled A/B on the FUSION POINT alone. Motivation: by default
# the embedding reaches the mixer only after a SHARED 1x1 in_proj squeezes it to `channels` jointly with
# the ORF track and coverage (1280 -> 256 for RiNALMo, a 5x bottleneck the 4-d one-hot backend never
# pays), then 10 dilated convs rewrite it. The zero-init gate means training starts numerically
# identical to the baseline, so the learned |gate| is a direct readout of how much the model wants the
# un-convolved embedding.
# Training just READS the pre-computed per-token .npy (via the backend's tx_index) -> the rinalmo SIF (torch)
# works for every backend; the FM model itself is not loaded here. KOZAK never set (no-Kozak track).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
TAG=${TAG:-mm1}
BACKEND=${BACKEND:?set BACKEND=rinalmo|orthrus|hydrarna}
FM_TO_MIXER=${FM_TO_MIXER:-0}
FMARG=""; SUFFIX=""
[ "$FM_TO_MIXER" = 1 ] && { FMARG="--fm_to_mixer"; SUFFIX="_fmmix"; }
OUT=$NEW/results/loto/orf_v2_attn_${BACKEND}_noBrain_nokozak_${TAG}_holdout_Hepatocytes${SUFFIX}
mkdir -p "$NEW/logs/train" "$OUT"
echo "BACKEND=$BACKEND TAG=$TAG fm_to_mixer=$FM_TO_MIXER holdout=Hepatocytes host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes \
    --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC \
    --emb_backend "$BACKEND" --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    --use_orf_track --n_attn_layers 2 $FMARG \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "train_done=$(date -Iseconds) -> $OUT"
