#!/usr/bin/env bash
#SBATCH --job-name=ribo_loto_1hot
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_1hot_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_1hot_%A_%a.err
set -euo pipefail
# ONE-HOT control LOTO (hold out Hepatocytes, train on the other 8 tissues): identical experiment to
# train_loto_orthrus.sbatch but with --emb_backend onehot -- a (L,4) A/C/G/T one-hot in place of the FM
# per-token embedding. Measures how much lift the FM embeddings actually provide over raw sequence,
# holding the ORF-candidate track + RNA-seq coverage + architecture constant. baseline = onehot+coverage
# (no ORF track); orf_v2_attn = onehot + ORF track v2 + 2 attention layers (the direct FM-lift control,
# comparable to the rinalmo/orthrus orf_v2_attn Task 13/14 numbers). Splits are byte-identical to the FM
# runs (onehot universe = rinalmo intersect orthrus). Output dirs get an _onehot suffix.
HOLDOUT=Hepatocytes
TASKS=(baseline orf_v2_attn)
CFG=${TASKS[$SLURM_ARRAY_TASK_ID]}

BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed/orf_track.npy
case $CFG in
  baseline)     EXTRA="";                                  BUDGET=24000 ;;
  orf_v2_attn)  EXTRA="--use_orf_track --n_attn_layers 2"; BUDGET=16000; ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/loto/${CFG}_onehot_holdout_${HOLDOUT}
mkdir -p "$NEW/logs/train" "$OUT"

echo "cfg=$CFG backend=onehot holdout=$HOLDOUT budget=$BUDGET extra='$EXTRA' orftrack=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout "$HOLDOUT" \
    --emb_backend onehot --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    $EXTRA \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget "$BUDGET" --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0

echo "cfg=$CFG backend=onehot holdout=$HOLDOUT train_done=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 --device cuda || echo "WARN eval_extra failed $CFG onehot"
echo "cfg=$CFG backend=onehot holdout=$HOLDOUT extra_done=$(date -Iseconds)"
