#!/usr/bin/env bash
# Comet search for one pgx database, as an engine cross-check on MSFragger.
#
# WHY. MSFragger's `calibrate_mass = 2` derives a mass calibration from a FIRST-PASS search, so the
# calibration depends on the database. Measured on Microglia: swapping in 699 extra sequences moved
# 28,195 spectra in or out of the matched set and shifted canonical identifications by +11,017,
# which is not target-decoy competition at all. Turning calibration off removes the artefact
# (+16 spectra, dPSM 0) but costs 3x the canonical IDs (205,348 -> 68,411), so it is not usable.
#
# Comet is single-pass with no data-derived recalibration, so the only thing that differs between
# two arms is the database. If dPSM behaves sanely here, the anomaly is an MSFragger property and
# Comet is the right engine for database comparisons.
#
#   DB=<db.fasta> MZML="<glob>" OUT=<dir> sbatch proteogenomics/scripts/pgx/comet_search.sbatch
#
# Optional: PARAMS, KEEP_RAW=1 (skip the digest), FORCE=1, PGX_SCRIPTS.
#SBATCH --job-name=pgx_comet
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --output=%x_%j.out
set -euo pipefail

: "${DB:?set DB=<db.fasta>}"
: "${MZML:?set MZML=<glob of .mzML>}"
: "${OUT:?set OUT=<output dir>}"

if [[ -z "${PGX_SCRIPTS:-}" ]]; then
  _cmd=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | tr ' ' '\n' | sed -n 's/^Command=//p' | head -1)
  [[ -n "$_cmd" ]] || { echo "set PGX_SCRIPTS=<.../proteogenomics/scripts>" >&2; exit 2; }
  PGX_SCRIPTS=$(cd "$(dirname "$(dirname "$_cmd")")" && pwd)
fi
NEW=$(cd "$PGX_SCRIPTS/../.." && pwd)
PY=${PGX_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3}
COMET=${PGX_COMET:-/private/groups/carpenterlab/emalekos/conda_envs/proteomics/bin/comet}
PARAMS=${PARAMS:-$NEW/proteogenomics/msfragger/comet_macro_lfq.params}
[[ -s "$DB" ]] || { echo "ABORT: missing DB $DB" >&2; exit 1; }
[[ -x "$COMET" ]] || { echo "ABORT: no comet at $COMET" >&2; exit 1; }
mkdir -p "$OUT"

SCR=/data/tmp/emalekos/pgx_comet_$(basename "$OUT")_${SLURM_JOB_ID:-$$}
mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT

if [[ "${FORCE:-0}" != "1" && -s "$OUT/rank1.tsv.gz" ]]; then
  echo "SKIP $OUT (digest present)"; exit 0
fi

# Same staging gate as the MSFragger path: skip vendor stubs and truncated files.
n=0
shopt -s nullglob
FILES=($MZML)
for m in "${FILES[@]}"; do
  [[ $(stat -c%s "$m" 2>/dev/null || echo 0) -gt 1000000 ]] || { echo "skip stub $(basename "$m")"; continue; }
  tail -c 400 "$m" 2>/dev/null | grep -qE "</indexedmzML>|</mzML>" || { echo "skip truncated $(basename "$m")"; continue; }
  cp "$m" "$SCR/"; n=$((n + 1))
done
[[ $n -gt 0 ]] || { echo "ABORT: no mzML staged from $MZML" >&2; exit 1; }
sed "s|^database_name.*|database_name = $DB|" "$PARAMS" > "$SCR/comet.params"
echo "staged $n mzML  db=$(basename "$DB")  engine=comet  host=$(hostname)  start=$(date -Iseconds)"

cd "$SCR"
"$COMET" -P"$SCR/comet.params" "$SCR"/*.mzML

PYTHONPATH="$PGX_SCRIPTS" "$PY" -m pgx.digest --search-dir "$SCR" --engine comet \
  --out "$OUT/rank1.tsv.gz" || { echo "ABORT: digest failed" >&2; exit 1; }
[[ -s "$OUT/rank1.tsv.gz" ]] || { echo "ABORT: empty digest" >&2; exit 1; }
cp "$SCR/comet.params" "$OUT/comet.params"
echo "DONE $(date -Iseconds): $(du -h "$OUT/rank1.tsv.gz" | cut -f1) digest in $OUT"
