#!/usr/bin/env bash
#SBATCH --job-name=kozak_eval
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --array=0-3%4
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/kozak_eval_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/kozak_eval_%A_%a.err
set -euo pipefail
# Kozak arms eval (GPU): per arm, dump the predicted per-nt profiles (-> dropin/pred_profiles.npz, for
# the ATG+CTG drop-in) AND run the localization eval (uORF / non-canonical AUROC) on the held-out
# Hepatocytes test set. CRITICAL: each arm re-runs inference, so RIBO_ORF_TRACK MUST be the arm's OWN
# training track (v0 heuristic / v1 nokozak / v2 pwm / v3 nokozak). The model config (incl. v3's
# learn_start_context) is read from the run's args.json. See KOZAK_PLAN.md.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
TASKS=(v0_heuristic v1_nokozak v2_pwm v3_learned)
CFG=${TASKS[$SLURM_ARRAY_TASK_ID]}
case $CFG in
  v0_heuristic) ORFTRACK=$NEW/data/packed/orf_track_v2.npy ;;
  v1_nokozak)   ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy ;;
  v2_pwm)       ORFTRACK=$NEW/data/packed/orf_track_v2_pwm.npy ;;
  v3_learned)   ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy ;;
  *) echo "bad CFG $CFG"; exit 1 ;;
esac
RUN=$NEW/results/kozak/${CFG}_onehot_fib2hep
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
mkdir -p "$NEW/logs/eval"
[ -f "$RUN/best.pt" ] || { echo "ERROR: no best.pt in $RUN (training not finished?)"; exit 1; }
echo "cfg=$CFG orftrack=$ORFTRACK host=$(hostname) start=$(date -Iseconds)"

echo "[$(date +%T)] $CFG: dump pred_profiles"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$RUN" --device cuda

echo "[$(date +%T)] $CFG: localization eval"
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$RUN" --min_orf_nt 30 --eval_cap 0 --device cuda
echo "cfg=$CFG done=$(date -Iseconds)"
