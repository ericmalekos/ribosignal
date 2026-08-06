#!/usr/bin/env bash
#SBATCH --job-name=loto_target
#SBATCH --partition=medium
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --array=0-7%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/build_target_loto_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/build_target_loto_%A_%a.err
set -euo pipefail
# Per-tissue pooled per-nt P-site target for the 8 non-Fibroblast Chothani tissues (Fibroblast
# already built). Pools the tissue's existing RiboCode per-nt P-site hd5 (posture A). Cheap CPU.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
TISSUES=(Brain ES Fat HA_EC HCAEC HUVEC Hepatocytes VSMC)
T=${TISSUES[$SLURM_ARRAY_TASK_ID]}
echo "tissue=$T host=$(hostname) start=$(date -Iseconds)"
$PY "$NEW/scripts/build_psite_target.py" "$T"
echo "tissue=$T done=$(date -Iseconds)"
