#!/usr/bin/env bash
#SBATCH --job-name=loto_arch
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_arch_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_arch_%j.err
set -euo pipefail
# Architecture variants of the current FINAL model (one-hot + mm1 coverage + no-Kozak track, noBrain 7 tissues,
# holdout Hepatocytes) -- ONLY the global-context mixer differs, for a controlled comparison:
#   ARCH=attn4  -> expanded transformer, 4 attention layers (vs the deployed 2), 8 heads
#   ARCH=mamba  -> bidirectional Mamba mixer, 2 layers (drop-in swap for the 2 transformer layers)
# Everything else is byte-identical to train_loto_noBrain.sbatch. KOZAK never set (no-Kozak track).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
TAG=${TAG:-mm1}
SEED=${SEED:-0}                 # 0 = the original single-seed runs; 1,2,... for the seed-robustness confirm
ARCH=${ARCH:?set ARCH=attn2|attn4|mamba}
# attn2/attn4 are pure transformer (rinalmo SIF); mamba needs mamba_ssm -> the orthrus SIF (the mamba ablation's image)
case "$ARCH" in
  attn2) ARCHARGS="--n_attn_layers 2 --n_heads 8";   NAME=orf_v2_attn;  SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  attn4) ARCHARGS="--n_attn_layers 4 --n_heads 8";   NAME=orf_v2_attn4; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba) ARCHARGS="--n_attn_layers 2 --mixer mamba"; NAME=orf_v2_mamba;  SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  mamba4) ARCHARGS="--n_attn_layers 4 --mixer mamba"; NAME=orf_v2_mamba4; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  *) echo "bad ARCH=$ARCH (want attn2|attn4|mamba|mamba4)" >&2; exit 1 ;;
esac
SUF=""; [ "$SEED" = 0 ] || SUF="_seed${SEED}"
OUT=$NEW/results/loto/${NAME}_onehot_noBrain_nokozak_${TAG}_holdout_Hepatocytes${SUF}
mkdir -p "$NEW/logs/train" "$OUT"
echo "ARCH=$ARCH ($ARCHARGS) TAG=$TAG SEED=$SEED holdout=Hepatocytes host=$(hostname) start=$(date -Iseconds)"
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
    --emb_backend onehot --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    --use_orf_track $ARCHARGS \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0
echo "train_done=$(date -Iseconds) -> $OUT"
