#!/usr/bin/env bash
#SBATCH --job-name=loto_pack
#SBATCH --partition=medium
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --array=0-7%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pack_loto_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pack_loto_%A_%a.err
set -euo pipefail
# Pack each non-Fibroblast tissue's target + pooled coverage onto the shared 36,668-tx
# universe axis -> data/packed_<Tissue>/ (+ per-tissue coverage_norm.json). Runs after both
# build_target_loto and pool_coverage_loto (afterok): needs the tissue's target hd5 + pooled
# coverage hd5. offsets/lengths/tx_order come out identical across tissues (same universe),
# so the shared orf_track + embeddings align to every pack.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
TISSUES=(Brain ES Fat HA_EC HCAEC HUVEC Hepatocytes VSMC)
T=${TISSUES[$SLURM_ARRAY_TASK_ID]}
echo "tissue=$T host=$(hostname) start=$(date -Iseconds)"
$PY "$NEW/scripts/pack_target_coverage_tissue.py" "$T"
echo "tissue=$T done=$(date -Iseconds)"
