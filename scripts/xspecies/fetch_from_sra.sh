#!/usr/bin/env bash
# Fetch runs from NCBI SRA instead of ENA, for the files ENA cannot serve.
#
# WHY: on 2026-08-25/26 ENA's FTP was down globally and its HTTPS layer answered some
# file paths with an ~800 byte HTML directory listing instead of the file. Four
# zebrafish runs could not be retrieved from ENA at all. NCBI SRA holds the same runs
# and is an independent source.
#
# THE VALIDATION PROBLEM, AND HOW IT IS SOLVED.
# ENA publishes an md5 per fastq.gz, and that is what every other file in this project
# was verified against. An SRA-derived FASTQ CANNOT match it: fasterq-dump re-generates
# the file, so gzip framing, read order and header formatting all differ even when the
# underlying reads are identical. Comparing to the ENA md5 here would fail on a
# perfectly good file, and silently accepting "md5 differs, probably fine" is exactly
# the habit this project is trying not to have. So two other checks are used instead:
#
#   1. READ COUNT must equal ENA's published `read_count` for the run. That is an
#      independent number from the run table, not from the bytes just downloaded.
#   2. If a mate of the same run was already retrieved from ENA and md5-verified, the
#      SRA-derived mate must carry the SAME read IDs in the same order. That is a much
#      stronger check than a count, and it is available for exactly the case that
#      motivated this script (SRR2087765, whose _1 verified from ENA and whose _2 did not).
#
# Provenance is recorded per file in data/external/xspecies/sra_provenance.tsv, because
# these files are NOT ENA-md5-verified like everything else and that difference must be
# visible later.
#
# Usage: bash scripts/xspecies/fetch_from_sra.sh <label> <run> [<run> ...]
set -uo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SRA=/private/groups/carpenterlab/emalekos/conda_envs/sratools/bin
XDIR=$NEW/data/external/xspecies
FQROOT=$XDIR/fastq
PROV=$XDIR/sra_provenance.tsv
LOG=$NEW/logs/xspecies/fetch_sra.log
# Group fs, not /data/tmp: /data/tmp is node-local, and anything a later step reads
# must be visible from other nodes.
CACHE=$XDIR/.sra_cache
mkdir -p "$CACHE" "$(dirname "$LOG")"

LABEL=${1:?usage: fetch_from_sra.sh <label> <run> [run ...]}; shift
RUNS=("$@"); [ "${#RUNS[@]}" -gt 0 ] || { echo "no runs given" >&2; exit 2; }
OUT="$FQROOT/$LABEL"; mkdir -p "$OUT"
[ -s "$PROV" ] || printf 'utc\tlabel\trun\tfile\tbytes\tmd5\tsource\tvalidation\tread_count\n' > "$PROV"

# ENA's published read_count, for check 1.
expected_reads() {
  local run=$1 f
  for f in "$XDIR"/PRJ*_runs.tsv; do
    awk -F'\t' -v r="$run" 'NR==1{for(i=1;i<=NF;i++){if($i=="run_accession")ra=i; if($i=="read_count")rc=i}}
                            NR>1 && $ra==r {print $rc; exit}' "$f"
  done | head -1
}

echo "=== SRA fetch label=$LABEL runs=${RUNS[*]} start=$(date -Iseconds)" | tee -a "$LOG"
fail=0
for run in "${RUNS[@]}"; do
  exp=$(expected_reads "$run")
  echo "  [$run] ENA read_count=${exp:-unknown}" | tee -a "$LOG"

  SCR=$(mktemp -d "$CACHE/${run}.XXXXXX")
  echo "  [$run] prefetch $(date +%H:%M:%S)" | tee -a "$LOG"
  if ! "$SRA/prefetch" --max-size u --output-directory "$SCR" "$run" >>"$LOG" 2>&1; then
    echo "  [$run] PREFETCH FAILED" | tee -a "$LOG"; rm -rf "$SCR"; fail=$((fail+1)); continue
  fi
  echo "  [$run] fasterq-dump $(date +%H:%M:%S)" | tee -a "$LOG"
  if ! "$SRA/fasterq-dump" --split-files --threads 6 --temp "$SCR" \
        --outdir "$SCR" "$SCR/$run" >>"$LOG" 2>&1; then
    echo "  [$run] FASTERQ-DUMP FAILED" | tee -a "$LOG"; rm -rf "$SCR"; fail=$((fail+1)); continue
  fi

  shopt -s nullglob
  produced=("$SCR/${run}"*.fastq)
  shopt -u nullglob
  [ "${#produced[@]}" -gt 0 ] || { echo "  [$run] no fastq produced" | tee -a "$LOG"; rm -rf "$SCR"; fail=$((fail+1)); continue; }

  ok_run=1
  for fq in "${produced[@]}"; do
    base=$(basename "$fq"); n=$(( $(wc -l < "$fq") / 4 ))
    echo "    $base: $n reads" | tee -a "$LOG"
    if [ -n "$exp" ] && [ "$n" != "$exp" ]; then
      echo "    MISMATCH: $base has $n reads, ENA says $exp" | tee -a "$LOG"; ok_run=0
    fi
  done
  if [ "$ok_run" != 1 ]; then
    echo "  [$run] read-count check FAILED, nothing installed" | tee -a "$LOG"
    rm -rf "$SCR"; fail=$((fail+1)); continue
  fi

  # Cross-check BEFORE installing anything, against the mate as it exists on disk NOW.
  # The first version installed _1 first and then compared _2 against it -- by which
  # point _1 was itself the SRA copy, so the check was circular AND it had already
  # overwritten an ENA-md5-verified file. Never clobber a verified file with an
  # unverifiable one.
  declare -A XCHECK=()
  for fq in "${produced[@]}"; do
    base=$(basename "$fq"); mate="${base%_*}"; idx="${base##*_}"; idx="${idx%.fastq}"
    other=""
    [ "$idx" = "2" ] && other="$OUT/${mate}_1.fastq.gz"
    [ "$idx" = "1" ] && other="$OUT/${mate}_2.fastq.gz"
    if [ -n "$other" ] && [ -s "$other" ]; then
      a=$(awk 'NR%4==1{sub(/ .*/,"");sub(/\/[12]$/,"");print}' "$fq" | head -200000 | md5sum | cut -d' ' -f1)
      b=$(zcat "$other" | awk 'NR%4==1{sub(/ .*/,"");sub(/\/[12]$/,"");print}' | head -200000 | md5sum | cut -d' ' -f1)
      if [ "$a" = "$b" ]; then
        echo "    $base: read IDs match the on-disk mate over the first 200k reads" | tee -a "$LOG"
        XCHECK[$base]="read_ids==on_disk_mate(200k)"
      else
        echo "    $base: READ IDs DO NOT MATCH the on-disk mate -- refusing" | tee -a "$LOG"; ok_run=0
      fi
    fi
  done
  [ "$ok_run" = 1 ] || { echo "  [$run] cross-check FAILED, nothing installed" | tee -a "$LOG"; rm -rf "$SCR"; fail=$((fail+1)); continue; }

  for fq in "${produced[@]}"; do
    base=$(basename "$fq"); dest="$OUT/${base}.gz"
    # Refuse to replace a file that already verifies against its ENA md5. An SRA-derived
    # file can never reproduce that md5, so overwriting would DOWNGRADE the provenance.
    if [ -s "$dest" ]; then
      want=$(awk -F"\t" -v r="$run" -v b="${base}.gz" '$6==r{n=split($18,m,";"); split($17,f,";");
              for(i=1;i<=n;i++){fn=f[i]; sub(/.*\//,"",fn); if(fn==b) print m[i]}}' \
              "$XDIR/download_manifest.tsv" | head -1)
      if [ -n "$want" ] && [ "$(md5sum "$dest" | cut -d" " -f1)" = "$want" ]; then
        echo "    $base: already on disk and ENA-md5-verified, KEEPING the ENA copy" | tee -a "$LOG"
        continue
      fi
    fi
    validation="read_count==ENA($exp)${XCHECK[$base]:+; ${XCHECK[$base]}}"
    "$SRA/../bin/gzip" -c "$fq" > "$dest" 2>/dev/null || gzip -c "$fq" > "$dest"
    printf '%s\t%s\t%s\t%s\t%s\t%s\tNCBI_SRA\t%s\t%s\n' \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$LABEL" "$run" "${base}.gz" \
      "$(stat -c%s "$dest")" "$(md5sum "$dest" | cut -d' ' -f1)" "$validation" "$exp" >> "$PROV"
    echo "    installed $dest ($(du -h "$dest" | cut -f1))" | tee -a "$LOG"
  done
  [ "$ok_run" = 1 ] || fail=$((fail+1))
  rm -rf "$SCR"
done

echo "=== SRA fetch done=$(date -Iseconds) failures=$fail" | tee -a "$LOG"
echo "provenance: $PROV"
[ "$fail" -eq 0 ] || exit 1
