#!/usr/bin/env bash
#SBATCH --job-name=ho_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/pack_%x_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/pack_%x_%j.err
set -euo pipefail
# Pool an external held-out dataset's per-SRR P-sites + RNAseq coverage into a model-input pack.
# Pass the dataset name via --export=ALL,DS=<dataset>. h5py -> cas12a env (not a SIF).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/logs/heldout"
DS=${DS:?set DS=<dataset> via --export}
echo "dataset=$DS host=$(hostname) start=$(date -Iseconds)"
"$PY" "$NEW/scripts/heldout/build_heldout_pack.py" "$DS"
echo "done=$(date -Iseconds)"
