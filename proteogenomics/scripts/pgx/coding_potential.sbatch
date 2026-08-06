#!/usr/bin/env bash
# CPAT + CPC2 coding-potential DB arms for one dataset's expressed universe (locked decision D3).
#   UNIVERSE_FA=<fa> SPECIES=human LABEL=<name> OUT=<dir> sbatch coding_potential.sbatch
#SBATCH --job-name=pgx_cp
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=medium
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/pgx/%x_%j.out
set -euo pipefail
: "${UNIVERSE_FA:?}" "${OUT:?}"
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
echo "label=${LABEL:-?} host=$(hostname) start=$(date -Iseconds)"
cd "$NEW/proteogenomics/scripts"
"$PY" -m pgx.coding_potential --universe-fa "$UNIVERSE_FA" --species "${SPECIES:-human}" \
  --out "$OUT" --tools "${TOOLS:-cpat,cpc2}" --min-aa "${MIN_AA:-7}" --starts "${STARTS:-ATG}"
echo "done=$(date -Iseconds)"
