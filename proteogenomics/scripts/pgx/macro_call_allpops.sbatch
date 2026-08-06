#!/usr/bin/env bash
# Stage 1 of the macrophage cross-subtype experiment: calibrate theta, call ORFs (both arms) and
# scan N-terminal extensions for ONE macrophage population per array task. No databases, no
# searches -- those are built once from the pooled calls in stage 2, because a per-population
# database manufactures false tissue specificity (a peptide can only be found where its database
# contains it; a previous pass found 96% of "population-unique" peptides were that artifact).
#
#   MODEL=mamba4 sbatch --array=0-11 proteogenomics/scripts/pgx/macro_call_allpops.sbatch
#SBATCH --job-name=pgx_macro_call
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=medium
#SBATCH --output=%x_%A_%a.out
set -euo pipefail

MODEL=${MODEL:-mamba4}
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
PY=${PGX_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3}

mapfile -t POPS < <(ls "$PD"/universes/*_universe_tx.txt | xargs -n1 basename | sed 's/_universe_tx.txt//' | sort)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"

PROF=$PD/pred_${MODEL}_union/$POP/pred_profiles.npz
OUT=$PD/pgx_macro_${MODEL}/$POP
[ -s "$PROF" ] || { echo "ABORT: no dump at $PROF" >&2; exit 1; }
echo "pop=$POP model=$MODEL host=$(hostname) start=$(date -Iseconds)"

cd "$NEW/proteogenomics/scripts"
"$PY" -m pgx.run \
  --species mouse --label "$POP" \
  --profiles "$PROF" \
  --universe-fa "$PD/universes/${POP}_universe.fa" \
  --out "$OUT" \
  --dbs model \
  --skip-searches
echo "pop=$POP done=$(date -Iseconds)"
