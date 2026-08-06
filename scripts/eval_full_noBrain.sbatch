#!/usr/bin/env bash
#SBATCH --job-name=eval_full_noBrain
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/full_noBrain_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/full_noBrain_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
mkdir -p "$NEW/logs/eval"
# Parameterized: RUN = run dir (relative to $NEW); ORFTRACK MUST match the track the model was trained on.
# Defaults eval the no-kozak mm20 Brain-drop model with the no-kozak ORF track. For mm1: RUN=...nokozak_mm1...
RUN=${RUN:-results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm20_holdout_Hepatocytes}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2_nokozak.npy}
echo "eval RUN=$RUN ORFTRACK=$(basename "$ORFTRACK") start=$(date -Iseconds)"
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta "/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img" bash -c '
    set -e
    python3 "'"$NEW"'/scripts/eval_localization.py" --run "'"$NEW/$RUN"'" --min_orf_nt 30 --eval_cap 0 --device cuda
    python3 "'"$NEW"'/scripts/dump_pred_profiles.py" --run "'"$NEW/$RUN"'" --device cuda
  '
echo "localization + dump done=$(date -Iseconds)"
