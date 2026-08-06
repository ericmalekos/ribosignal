#!/usr/bin/env bash
#SBATCH --job-name=kozak_arm
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=18:00:00
#SBATCH --array=0-3%4
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/kozak_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/kozak_%A_%a.err
set -euo pipefail
# Kozak start-context ablation: 4 arms, one-hot backend, TRAIN Fibroblast -> TEST Hepatocytes.
# Only ORF-track channel 3 (start propensity) differs (occupancy/stop identical), except v3 which uses
# the no-Kozak track + a learnable in-model start-context module. See KOZAK_PLAN.md.
#   v0_heuristic : orf_track_v2.npy         (0.5*[purine@-3]+0.5*[G@+4])
#   v1_nokozak   : orf_track_v2_nokozak.npy (codon weight only)
#   v2_pwm       : orf_track_v2_pwm.npy     (empirical Kozak PWM)
#   v3_learned   : orf_track_v2_nokozak.npy + --learn_start_context (learned conv over one-hot)
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
TASKS=(v0_heuristic v1_nokozak v2_pwm v3_learned)
CFG=${TASKS[$SLURM_ARRAY_TASK_ID]}
EXTRA="--use_orf_track --n_attn_layers 2"
case $CFG in
  v0_heuristic) ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;;
  v1_nokozak)   ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy ;;
  v2_pwm)       ORFTRACK=$NEW/data/packed/orf_track_v2_pwm.npy ;;
  v3_learned)   ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy; EXTRA="$EXTRA --learn_start_context" ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac

SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
OUT=$NEW/results/kozak/${CFG}_onehot_fib2hep
mkdir -p "$NEW/logs/train" "$OUT"
echo "cfg=$CFG orftrack=$ORFTRACK extra='$EXTRA' host=$(hostname) start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes --train_tissues Fibroblast \
    --emb_backend onehot --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 --resume \
    $EXTRA \
    --epochs 40 --patience 8 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed 0

echo "cfg=$CFG train_done=$(date -Iseconds)"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0 || echo "WARN eval_extra failed $CFG"
echo "cfg=$CFG extra_done=$(date -Iseconds)"
