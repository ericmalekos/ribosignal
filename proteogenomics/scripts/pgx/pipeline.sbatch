#!/usr/bin/env bash
# Run the whole pgx pipeline for one dataset on a compute node, then submit the MSFragger searches.
# The light steps (universe, calibration sweep, calling, extensions, database build) run here;
# MSFragger is submitted separately because it needs hours and ~160 GB.
#
#   SPECIES=mouse LABEL=BMDM \
#   PROFILES=<...>/pred_profiles.npz \
#   SALMON="<...>/r1/quant.sf <...>/r2/quant.sf" \
#   MZML='<...>/mzml/BMDM_Rep*_F*.mzML' \
#   OUT=<...>/pgx/BMDM \
#     sbatch proteogenomics/scripts/pgx/pipeline.sbatch
#
# Optional: ENZYME (default tryptic), CDS_RECALL, PVAL, ARM, START_CODONS, DBS, MIN_AA,
#           BIOTYPE_FIELD, SHARED_SEARCH_ROOT, UNIVERSE_FA, EXTRA (raw extra args), PGX_SCRIPTS.
# Then, once the searches land:
#   ... REPORT_ONLY=1 sbatch proteogenomics/scripts/pgx/pipeline.sbatch
#SBATCH --job-name=pgx_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=110G
#SBATCH --time=4:00:00
#SBATCH --partition=medium
#SBATCH --output=%x_%j.out
set -euo pipefail

: "${SPECIES:?set SPECIES=human|mouse}"
: "${LABEL:?set LABEL=<dataset name>}"
: "${PROFILES:?set PROFILES=<pred_profiles.npz>}"
: "${OUT:?set OUT=<output dir>}"

if [[ -z "${PGX_SCRIPTS:-}" ]]; then
  _cmd=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | tr ' ' '\n' | sed -n 's/^Command=//p' | head -1)
  [[ -n "$_cmd" ]] || { echo "set PGX_SCRIPTS=<.../proteogenomics/scripts>" >&2; exit 2; }
  PGX_SCRIPTS=$(cd "$(dirname "$(dirname "$_cmd")")" && pwd)
fi
[[ -d "$PGX_SCRIPTS/pgx" ]] || { echo "no pgx package under $PGX_SCRIPTS" >&2; exit 2; }
PY=${PGX_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3}

ARGS=(--species "$SPECIES" --label "$LABEL" --profiles "$PROFILES" --out "$OUT"
      --enzyme "${ENZYME:-tryptic}")
[[ -n "${SALMON:-}" ]]             && { read -r -a _s <<< "$SALMON"; ARGS+=(--salmon "${_s[@]}"); }
[[ -n "${UNIVERSE_FA:-}" ]]        && ARGS+=(--universe-fa "$UNIVERSE_FA")
[[ -n "${MZML:-}" ]]               && ARGS+=(--mzml "$MZML")
[[ -n "${CDS_RECALL:-}" ]]         && ARGS+=(--cds-recall "$CDS_RECALL")
[[ -n "${PVAL:-}" ]]               && ARGS+=(--pval "$PVAL")
[[ -n "${ARM:-}" ]]                && ARGS+=(--arm "$ARM")
[[ -n "${START_CODONS:-}" ]]       && ARGS+=(--start-codons "$START_CODONS")
[[ -n "${DBS:-}" ]]                && ARGS+=(--dbs "$DBS")
[[ -n "${MIN_AA:-}" ]]             && ARGS+=(--min-aa "$MIN_AA")
[[ -n "${BIOTYPE_FIELD:-}" ]]      && ARGS+=(--biotype-field "$BIOTYPE_FIELD")
[[ -n "${SHARED_SEARCH_ROOT:-}" ]] && ARGS+=(--shared-search-root "$SHARED_SEARCH_ROOT")
[[ "${REPORT_ONLY:-0}" == "1" ]]   && ARGS+=(--report-only)
[[ "${SKIP_SEARCHES:-0}" == "1" ]] && ARGS+=(--skip-searches)
[[ "${SUBMIT:-1}" == "1" && "${REPORT_ONLY:-0}" != "1" ]] && ARGS+=(--submit-searches)
[[ -n "${EXTRA:-}" ]]              && { read -r -a _e <<< "$EXTRA"; ARGS+=("${_e[@]}"); }

echo "host=$(hostname) start=$(date -Iseconds)"
cd "$PGX_SCRIPTS"
"$PY" -m pgx.run "${ARGS[@]}"
echo "done=$(date -Iseconds)"
