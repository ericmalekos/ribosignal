#!/usr/bin/env bash
# Poisson/theta calibration sweep for the pgx pipeline: one theta per array task.
# Nothing is hardcoded; every input arrives by environment variable.
#
#   PROFILES=<pred_profiles.npz> SPECIES=mouse OUTDIR=<sweep dir> \
#     sbatch proteogenomics/scripts/pgx/sweep_theta.sbatch
#
# Optional: THETAS (space separated, default below), MIN_AA, PVAL, SEED, ANNOT, NO_POISSON=1.
# After the array finishes, run the pick step:
#   python -m pgx.calibrate --profiles $PROFILES --species $SPECIES --out $OUTDIR --pick-only
#SBATCH --job-name=pgx_theta
#SBATCH --array=0-6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --output=%x_%A_%a.out
set -euo pipefail

: "${PROFILES:?set PROFILES=<pred_profiles.npz>}"
: "${SPECIES:?set SPECIES=human|mouse}"
: "${OUTDIR:?set OUTDIR=<sweep dir>}"

read -r -a TH <<< "${THETAS:-0.02 0.05 0.10 0.20 0.35 0.60 1.00}"
T=${TH[$SLURM_ARRAY_TASK_ID]}

# SLURM executes a COPY of this script from its spool dir, so BASH_SOURCE/$0 do NOT point at the
# repo. Resolve the package dir from PGX_SCRIPTS if set, else ask slurm for the submitted path.
if [[ -z "${PGX_SCRIPTS:-}" ]]; then
  _cmd=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | tr ' ' '\n' | sed -n 's/^Command=//p' | head -1)
  [[ -n "$_cmd" ]] || { echo "cannot resolve script path; set PGX_SCRIPTS=<.../proteogenomics/scripts>" >&2; exit 2; }
  PGX_SCRIPTS=$(cd "$(dirname "$(dirname "$_cmd")")" && pwd)
fi
PKG=$PGX_SCRIPTS                                     # .../proteogenomics/scripts (holds pgx/)
[[ -d "$PKG/pgx" ]] || { echo "no pgx package under $PKG" >&2; exit 2; }
PY=${PGX_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3}

EXTRA=()
[[ -n "${ANNOT:-}" ]] && EXTRA+=(--annot "$ANNOT")
[[ "${NO_POISSON:-0}" == "1" ]] && EXTRA+=(--no-poisson)

echo "theta=$T host=$(hostname) start=$(date -Iseconds)"
cd "$PKG"
"$PY" -m pgx.calibrate \
  --profiles "$PROFILES" --species "$SPECIES" --out "$OUTDIR" \
  --only-theta "$T" --min-aa "${MIN_AA:-5}" --pval "${PVAL:-0.05}" --seed "${SEED:-0}" \
  "${EXTRA[@]}"
echo "theta=$T done=$(date -Iseconds)"
