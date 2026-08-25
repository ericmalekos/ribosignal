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
# Pass DS=<dataset> via --export. Optional RIBO_COVERAGE_DIR / RIBO_PACK_OUT_SUFFIX select a
# variant coverage posture (e.g. mm10) without editing DATASETS. h5py -> riboseq env.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
# cas12a has torch but NO h5py -- this script imports h5py, so it must use riboseq.
PY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python
mkdir -p "$NEW/logs/heldout"
DS=${DS:?set DS=<dataset> via --export}
echo "dataset=$DS host=$(hostname) start=$(date -Iseconds)"
"$PY" "$NEW/scripts/heldout/build_heldout_pack.py" "$DS"
echo "done=$(date -Iseconds)"
