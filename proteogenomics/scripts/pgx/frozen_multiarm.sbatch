#!/usr/bin/env bash
# Search every database arm of one population under IDENTICAL, frozen search parameters.
#
# THE PROBLEM THIS SOLVES. With calibrate_mass = 2, MSFragger re-derives six search parameters
# from a first-pass search against whichever database it is handed. Arms of the same population
# are then scored under different rules, so canonical counts move for reasons unrelated to the
# sequences added and dPSM stops measuring target-decoy competition. Freezing the parameters
# (fragger_macro_lfq_frozen.params, calibrate_mass = 0) removes every data-derived step, so all
# arms share one scoring regime -- verified on one Microglia fraction as a shared score cutoff of
# 18.598 across arms, dPSM +5, and no sensitivity lost versus calibrate_mass = 2.
#
# TWO EARLIER ROUTES THAT DO NOT WORK, recorded so they are not retried:
#   write_calibrated_mzml = 1 -> writes MS1 spectra ONLY (6,051 of 77,289 on a test fraction).
#                                Not searchable; it is an input for precursor-level quant.
#   write_mzbin_all = 1       -> writes UNCALIBRATED spectra (matches calibrate_mass = 0 exactly,
#                                2,336 vs 2,333 targets, identical cutoff). A format cache only.
# The .mzBIN_calibrated MSFragger builds internally is deleted before the run ends and no
# parameter retains it. With frozen parameters none of this matters -- there is nothing to share.
#
# WHY mzBIN IS STILL USED. Converting mzML is the expensive step (minutes per fraction); searching
# the cache took 10 s. So convert once per population, then search every arm against the cache.
# This also guarantees all arms see byte-identical spectra. MSFragger rejects the .mzBIN_all
# extension it writes ("Unrecognized file extension: mzbin_all"), so the files are renamed .mzBIN.
#
#   POPS_FILE=<pop_mzml_map.tsv> DBDIR=<dir> MZMLDIR=<dir> OUT=<dir> ARMS="a b c" \
#     sbatch --array=0-11 proteogenomics/scripts/pgx/frozen_multiarm.sbatch
#
# Optional: FORCE=1 (re-run arms with a digest present), XMX, PGX_SCRIPTS.
#SBATCH --job-name=pgx_frozen
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=long
#SBATCH --output=%x_%A_%a.out
set -euo pipefail

: "${POPS_FILE:?}" "${OUT:?}" "${MZMLDIR:?}"
# Databases are per-population: pass DBROOT (resolved to $DBROOT/$POP) for a whole array, or
# DBDIR to pin one directory for every task.
: "${DBROOT:=}" "${DBDIR:=}"
[[ -n "$DBROOT" || -n "$DBDIR" ]] || { echo "ABORT: set DBROOT=<dir> or DBDIR=<dir>" >&2; exit 2; }
ARMS=${ARMS:-"gencode model_predicted model_ribocode_bmdm_nt model_ribocode_bmdm null_atg null_nc"}

if [[ -z "${PGX_SCRIPTS:-}" ]]; then
  _cmd=$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | tr ' ' '\n' | sed -n 's/^Command=//p' | head -1)
  PGX_SCRIPTS=$(cd "$(dirname "$(dirname "$_cmd")")" && pwd)
fi
NEW=$(cd "$PGX_SCRIPTS/../.." && pwd)
PY=${PGX_PYTHON:-/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3}
JAR=${PGX_MSFRAGGER_JAR:-/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar}
JAVA=${PGX_JAVA:-/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java}
TEMPLATE=${PGX_TEMPLATE:-$NEW/proteogenomics/msfragger/fragger_macro_lfq_frozen.params}
[[ -s "$TEMPLATE" ]] || { echo "ABORT: no frozen template at $TEMPLATE" >&2; exit 1; }

# POPS_FILE must live on the GROUP filesystem: /data/tmp is node-local, so a file written on the
# head node is invisible to the compute node that runs the array task.
mapfile -t ROWS < <(sort "$POPS_FILE")
IFS=$'\t' read -r POP PREF <<< "${ROWS[$SLURM_ARRAY_TASK_ID]:?no row}"
[[ -n "$DBROOT" ]] && DBDIR="$DBROOT/$POP"
[[ -d "$DBDIR" ]] || { echo "ABORT: no database dir $DBDIR" >&2; exit 1; }
SCR=/data/tmp/emalekos/pgx_frozen_${POP}_${SLURM_JOB_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "pop=$POP host=$(hostname) template=$(basename "$TEMPLATE") start=$(date -Iseconds)"

# ---- stage mzML (skip vendor stubs and truncated files) ----
n=0; shopt -s nullglob
for m in "$MZMLDIR"/${PREF}_Rep*_F*.mzML; do
  [[ $(stat -c%s "$m" 2>/dev/null || echo 0) -gt 1000000 ]] || continue
  tail -c 400 "$m" | grep -qE "</indexedmzML>|</mzML>" || continue
  cp "$m" "$SCR/"; n=$((n + 1))
done
[[ $n -gt 0 ]] || { echo "ABORT: no mzML staged for $POP" >&2; exit 1; }
echo "staged $n mzML"

# ---- convert once: mzML -> mzBIN (results discarded; only the cache is wanted) ----
cd "$PGX_SCRIPTS"
"$PY" -m pgx.search --db "$DBDIR/db_gencode.fasta" --enzyme tryptic --template "$TEMPLATE" \
  --set write_mzbin_all=1 --out-params "$SCR/convert.params"
cd "$SCR"
echo "--- mzBIN conversion $(date -Iseconds) ---"
"$JAVA" "-Xmx${XMX:-180g}" -jar "$JAR" "$SCR/convert.params" "$SCR"/*.mzML
for b in "$SCR"/*.mzBIN_all; do mv -f "$b" "${b%.mzBIN_all}.mzBIN"; done
BIN=("$SCR"/*.mzBIN)
# Gate on CONTENT, not just presence: a file can exist and still be unusable.
[[ ${#BIN[@]} -eq $n ]] || { echo "ABORT: $n staged but ${#BIN[@]} mzBIN written" >&2; exit 1; }
for b in "${BIN[@]}"; do
  [[ $(stat -c%s "$b") -gt 1000000 ]] || { echo "ABORT: $b is truncated" >&2; exit 1; }
done
echo "mzBIN cache: ${#BIN[@]} files, $(du -sh "$SCR" | cut -f1) scratch"
rm -f "$SCR"/*.mzML "$SCR"/*.tsv "$SCR"/*.pin "$SCR"/*.pepXML

# ---- every arm against the SAME cache, identical frozen parameters ----
for arm in $ARMS; do
  db="$DBDIR/db_${arm}.fasta"
  [[ -s "$db" ]] || { echo "SKIP $arm (no $db)"; continue; }
  o="$OUT/$POP/$arm"; mkdir -p "$o"
  if [[ "${FORCE:-0}" != "1" && -s "$o/rank1.tsv.gz" ]]; then echo "SKIP $arm (digest present)"; continue; fi
  cd "$PGX_SCRIPTS"
  "$PY" -m pgx.search --db "$db" --enzyme tryptic --template "$TEMPLATE" \
    --out-params "$SCR/${arm}.params"
  cd "$SCR"
  echo "--- arm $arm $(date -Iseconds) ---"
  rm -f "$SCR"/*.tsv "$SCR"/*.pin "$SCR"/*.pepXML
  "$JAVA" "-Xmx${XMX:-180g}" -jar "$JAR" "$SCR/${arm}.params" "${BIN[@]}"
  ls "$SCR"/*.tsv >/dev/null 2>&1 || { echo "ABORT: $arm produced no tsv" >&2; exit 1; }
  PYTHONPATH="$PGX_SCRIPTS" "$PY" -m pgx.digest --search-dir "$SCR" --out "$o/rank1.tsv.gz" \
    || { echo "ABORT: digest failed for $arm" >&2; exit 1; }
  cp "$SCR/${arm}.params" "$o/fragger.params"
  echo "$arm done $(du -h "$o/rank1.tsv.gz" | cut -f1)"
done
echo "pop=$POP ALL ARMS DONE $(date -Iseconds)"
