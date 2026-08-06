#!/usr/bin/env bash
# MSFragger search for one pgx database. Nothing hardcoded; inputs arrive by environment variable.
#
#   DB=<db.fasta> MZML="<glob or dir>" OUT=<dir> ENZYME=tryptic \
#     sbatch proteogenomics/scripts/pgx/search.sbatch
#
# Optional: TEMPLATE (params template), NUM_SLICES, XMX (default 140g), PGX_SCRIPTS, FORCE=1.
#
# Skips automatically when <OUT>/.pgx_hash already matches sha256(db)+sha256(params) and outputs
# are present, which is what stops the GENCODE baseline being re-searched once per arm.
#SBATCH --job-name=pgx_msf
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --output=%x_%j.out
set -euo pipefail

: "${DB:?set DB=<db.fasta>}"
: "${MZML:?set MZML=<dir or glob of .mzML>}"
: "${OUT:?set OUT=<output dir>}"
: "${ENZYME:?set ENZYME=tryptic|nonspecific}"

if [[ -z "${PGX_SCRIPTS:-}" ]]; then
  _cmd=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | tr ' ' '\n' | sed -n 's/^Command=//p' | head -1)
  [[ -n "$_cmd" ]] || { echo "set PGX_SCRIPTS=<.../proteogenomics/scripts>" >&2; exit 2; }
  PGX_SCRIPTS=$(cd "$(dirname "$(dirname "$_cmd")")" && pwd)
fi
[[ -d "$PGX_SCRIPTS/pgx" ]] || { echo "no pgx package under $PGX_SCRIPTS" >&2; exit 2; }

PY=${PGX_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3}
JAR=${PGX_MSFRAGGER_JAR:-/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar}
JAVA=${PGX_JAVA:-/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java}
[[ -s "$DB" ]] || { echo "ABORT: missing DB $DB" >&2; exit 1; }
mkdir -p "$OUT"

SCR=/data/tmp/emalekos/pgx_msf_$(basename "$OUT")_${SLURM_JOB_ID:-$$}
mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT
cd "$PGX_SCRIPTS"

EXTRA=()
[[ -n "${TEMPLATE:-}" ]]   && EXTRA+=(--template "$TEMPLATE")
[[ -n "${NUM_SLICES:-}" ]] && EXTRA+=(--num-slices "$NUM_SLICES")
# SET="key=value key2=value2" overrides individual MSFragger parameters. Used to test whether
# calibrate_mass makes the search database-dependent: with calibration on, MSFragger derives the
# mass calibration FROM a first-pass search, so a different database changes which spectra match
# at all in the second pass, and canonical counts move for reasons unrelated to the added sequences.
if [[ -n "${SET:-}" ]]; then read -r -a _sets <<< "$SET"; for kv in "${_sets[@]}"; do EXTRA+=(--set "$kv"); done; fi
"$PY" -m pgx.search --db "$DB" --enzyme "$ENZYME" --out-params "$SCR/fragger.params" "${EXTRA[@]}"
HASH=$("$PY" -m pgx.search --db "$DB" --enzyme "$ENZYME" --print-hash "${EXTRA[@]}")

# A cache hit requires the digest AND non-empty results. Testing only that *.tsv files EXIST is not
# enough: a copy that fails partway (a full filesystem, most likely) leaves zero-byte .tsv behind,
# which would then be honoured as a cached search and silently report no identifications.
_n_ok() { local n=0; for f in "$1"/*.tsv "$1"/rank1.tsv.gz; do [[ -s "$f" ]] && n=$((n + 1)); done; echo "$n"; }
if [[ "${FORCE:-0}" != "1" && -f "$OUT/.pgx_hash" && "$(cat "$OUT/.pgx_hash")" == "$HASH" \
      && $(_n_ok "$OUT") -gt 0 ]]; then
  echo "SKIP $OUT (cache hit $HASH, $(_n_ok "$OUT") non-empty tsv)"; exit 0
fi

# Stage COMPLETE mzML only: > 1MB (drops the sub-1MB vendor-bug stubs) AND a proper closing tag
# (drops OOM-truncated files). MSFragger aborts the whole batch on any corrupt input, so a size
# gate alone is not sufficient -- a truncated large file must be caught here.
n=0
shopt -s nullglob
if [[ -d "$MZML" ]]; then FILES=("$MZML"/*.mzML); else FILES=($MZML); fi
for m in "${FILES[@]}"; do
  [[ $(stat -c%s "$m" 2>/dev/null || echo 0) -gt 1000000 ]] || { echo "skip stub $(basename "$m")"; continue; }
  tail -c 400 "$m" 2>/dev/null | grep -qE "</indexedmzML>|</mzML>" || { echo "skip truncated $(basename "$m")"; continue; }
  cp "$m" "$SCR/"; n=$((n + 1))
done
echo "staged $n mzML  db=$(basename "$DB")  enzyme=$ENZYME  hash=$HASH  host=$(hostname)"
[[ $n -gt 0 ]] || { echo "ABORT: no mzML staged from $MZML" >&2; exit 1; }

cd "$SCR"
"$JAVA" "-Xmx${XMX:-140g}" -jar "$JAR" "$SCR/fragger.params" "$SCR"/*.mzML
# DIGEST=1 (default): reduce to rank-1 rows here, on the node, and copy only the ~1%-sized digest
# to the shared filesystem. The raw .tsv/.pin are ~1-5 GB per arm and nothing downstream reads more
# than rank-1 rows, so keeping them is pure quota cost. /data/tmp is node-local, so this reduction
# MUST happen in the search job; a later job on another node cannot reach the scratch.
if [[ "${DIGEST:-1}" == "1" ]]; then
  # PYTHONPATH is explicit because the working directory is $SCR by this point, not $PGX_SCRIPTS,
  # so `-m pgx.digest` would not resolve the package.
  PYTHONPATH="$PGX_SCRIPTS" "$PY" -m pgx.digest --search-dir "$SCR" --out "$OUT/rank1.tsv.gz" \
    || { echo "ABORT: digest failed" >&2; exit 1; }
  [[ -s "$OUT/rank1.tsv.gz" ]] || { echo "ABORT: empty digest" >&2; exit 1; }
else
  cp "$SCR"/*.pin "$SCR"/*.tsv "$OUT"/ || { echo "ABORT: copy of results failed (disk full?)" >&2; exit 1; }
fi
# pepXML is only an MS2Rescore input, and rescoring proved untrainable on this novel class
# (too few real novel peptides among the decoys), so it is optional and skipped when space is tight.
if [[ "${KEEP_PEPXML:-1}" == "1" ]]; then
  cp "$SCR"/*.pepXML "$OUT"/ 2>/dev/null || echo "WARN: pepXML not copied" >&2
fi
cp "$SCR/fragger.params" "$OUT/fragger.params"

# Only claim the search is cached once the results are actually on disk and non-empty. The earlier
# version copied with `|| true` and then wrote the digest unconditionally, so a filesystem that
# filled during copy-back produced zero-byte .tsv marked as a completed search.
if [[ "${DIGEST:-1}" == "1" ]]; then
  # One digest file covers every fraction, so the per-fraction count does not apply. Instead check
  # the digest actually contains rows from all staged fractions.
  n_frac=$(PYTHONPATH="$PGX_SCRIPTS" "$PY" -c "
import gzip,sys
seen=set()
with gzip.open('$OUT/rank1.tsv.gz','rt') as fh:
    next(fh)
    for ln in fh: seen.add(ln.split(':',1)[0])
print(len(seen))")
  if [[ "$n_frac" -lt "$n" ]]; then
    echo "ABORT: digest covers $n_frac of $n staged fractions; not writing .pgx_hash" >&2
    exit 1
  fi
  echo "$HASH" > "$OUT/.pgx_hash"
  echo "DONE $(date -Iseconds): digest covers $n_frac fractions ($(du -h "$OUT/rank1.tsv.gz" | cut -f1)) in $OUT"
else
  n_ok=$(_n_ok "$OUT")
  if [[ "$n_ok" -lt "$n" ]]; then
    echo "ABORT: only $n_ok of $n fractions have non-empty tsv; not writing .pgx_hash" >&2
    exit 1
  fi
  echo "$HASH" > "$OUT/.pgx_hash"
  echo "DONE $(date -Iseconds): $n_ok non-empty fraction tsv in $OUT"
fi
