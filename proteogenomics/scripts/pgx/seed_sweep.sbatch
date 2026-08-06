#!/usr/bin/env bash
# Poisson seed-sensitivity check: same theta, different RNG seeds. One seed per array task.
#   PROFILES=<npz> ANNOT=<dir> OUTDIR=<dir> THETA=0.05 sbatch pgx/seed_sweep.sbatch
#SBATCH --job-name=pgx_seed
#SBATCH --array=1-4
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --output=%x_%A_%a.out
set -euo pipefail
: "${PROFILES:?}" "${ANNOT:?}" "${OUTDIR:?}"
S=$SLURM_ARRAY_TASK_ID
RPY=${PGX_RIBOCODE_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python}
_cmd=$(scontrol show job "$SLURM_JOB_ID" | tr ' ' '\n' | sed -n 's/^Command=//p' | head -1)
ROOT=$(cd "$(dirname "$(dirname "$(dirname "$(dirname "$_cmd")")")")" && pwd)
mkdir -p "$OUTDIR/seed_$S"
"$RPY" "$ROOT/scripts/ribocode_dropin.py" --profiles "$PROFILES" --annot "$ANNOT" \
  --variant pred_preddepth --out "$OUTDIR/seed_$S" --min_aa 5 --pval 0.05 \
  --pred_scale "${THETA:-0.05}" --pred_poisson --seed "$S"
echo "seed=$S done"
