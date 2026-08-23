#!/usr/bin/env bash
#SBATCH --job-name=pB_xeval
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/pB_xeval_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/pB_xeval_%A_%a.err
set -euo pipefail
# Posture-B CROSS-EVAL: the cell Task 21 flagged as missing and never ran.
#
# Task 21 scored A on A's full fold-0 test set (6,962 tx) and B on B's representative subset
# (2,313 tx), so training-composition and eval-composition were confounded. This scores the
# B-TRAINED checkpoints on the FULL test set, giving A and B a COMMON denominator.
#
# It is also the direct test of "train restricted, predict broad": restrict training to the
# highest-expressed isoform per gene, then predict on everything passing the RNA threshold. No
# Ribo-seq is needed at inference, so the rule is deployable -- the open question is whether a
# model that only ever saw dominant isoforms smears predicted density onto the low-expressed
# siblings it was never taught to zero out.
#
# THE MECHANISM: dataset.load_split() calls representative_tx() at CALL time, which reads
# $RIBO_REPRESENTATIVE_TX. Leaving it UNSET yields the full posture-A test set from the same
# checkpoint. Nothing about the model changes; only the transcripts it is scored on.
#
# Writes to *_crossevalA/ copies so the Task 21 extra_metrics.json are NOT clobbered.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
RUNS=(baseline_B_onehot_f0_crossevalA orf_v2_attn_B_onehot_f0_crossevalA)
OUT=$NEW/results/postureB/${RUNS[$SLURM_ARRAY_TASK_ID]}
[ -f "$OUT/best.pt" ] || { echo "FATAL: no checkpoint at $OUT" >&2; exit 2; }
echo "run=$(basename "$OUT") host=$(hostname) start=$(date -Iseconds)"
echo "RIBO_REPRESENTATIVE_TX is deliberately UNSET -> full posture-A test set"
# NOTE: no --env RIBO_REPRESENTATIVE_TX, and --cleanenv guarantees it cannot leak in from the shell.
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$ORFTRACK" -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$OUT" --eval_cap 0
echo "run=$(basename "$OUT") done=$(date -Iseconds)"
