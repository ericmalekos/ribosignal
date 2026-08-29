#!/usr/bin/env bash
# Serial, md5-verified ENA fetch driven by data/external/xspecies/download_manifest.tsv.
#
# Encodes the same rules as scripts/fetch_heldout_rna.sh, which it is derived from:
#   HEAD NODE, SERIAL, ONE CONNECTION. Compute nodes throttle a single stream to
#   ~0.4 MB/s and parallel array downloads trip ENA's per-IP limit. Downloads are
#   I/O bound so they are a light head-node task even when they run for days.
#   `wget -nv -c`, never aria2c -x4: multi-connection silently corrupts large
#   files at exactly the right size. `-nv` not `-q`: -q suppresses wget's ERRORS
#   as well as progress, and a truncated transfer once left an empty log.
#   VERIFY BY md5, NEVER BY SIZE. On md5 failure DELETE and refetch from scratch;
#   `wget -c` would resume onto the bad bytes and append to them.
#
# Three differences from fetch_heldout_rna.sh, each fixing a real failure mode:
#
#   1. No per-accession ENA round trip. The URL and md5 come from the committed
#      manifest. fetch_heldout_rna.sh calls the ENA API inside `set -euo
#      pipefail`, so one transient curl failure aborts the whole chain and loses
#      every remaining accession. It also parses only the FIRST line of the
#      response, so a study accession silently fetches one run.
#   2. Redundant paired file skipped. ENA often lists a merged <ACC>.fastq.gz
#      alongside _1/_2; fetching all three wastes ~50% of the bytes. (Verified
#      absent from the projects in this manifest, but the guard stays: it costs
#      nothing and the next project will not be checked by hand.)
#   3. A row whose fastq_ftp is empty is a loud error, not a silent skip.
#
# Usage:
#   bash scripts/xspecies/fetch_from_manifest.sh [--dry-run] [label ...]
# With no labels, fetches every label in manifest order (smallest species first).
set -uo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MANIFEST="$NEW/data/external/xspecies/download_manifest.tsv"
DESTROOT="$NEW/data/external/xspecies/fastq"
LOGDIR="$NEW/logs/xspecies"
mkdir -p "$DESTROOT" "$LOGDIR"

DRY=0
[ "${1:-}" = "--dry-run" ] && { DRY=1; shift; }
WANT=("$@")

[ -s "$MANIFEST" ] || { echo "missing $MANIFEST" >&2; exit 2; }

# Column indices, read from the header so a manifest schema change fails loudly
# rather than fetching the wrong field.
# `exit` inside $( ) only leaves the subshell, so a missing column has to be
# checked in the caller or the fetch proceeds with an empty index and silently
# compares every md5 against the empty string. Observed once; do not simplify.
col() {
  awk -F'\t' -v want="$1" \
    'NR==1{for(i=1;i<=NF;i++){gsub(/\r$/,"",$i); if($i==want){print i; exit}}}' "$MANIFEST"
}
for want in label run_accession fastq_ftp fastq_md5 gib; do
  if [ -z "$(col "$want")" ]; then
    echo "FATAL: manifest $MANIFEST has no column '$want'" >&2; exit 2
  fi
done
C_LABEL=$(col label); C_RUN=$(col run_accession)
C_FTP=$(col fastq_ftp); C_MD5=$(col fastq_md5); C_GIB=$(col gib)

# Refuse a manifest with CRLF endings rather than fail every checksum on it.
if head -2 "$MANIFEST" | grep -q $'\r'; then
  echo "FATAL: $MANIFEST has CRLF line endings; md5 values would carry a trailing" >&2
  echo "       carriage return and every verification would fail. Regenerate with" >&2
  echo "       scripts/xspecies/select_runs.py (lineterminator must be \\n)." >&2
  exit 2
fi

labels=$(awk -F'\t' -v c="$C_LABEL" 'NR>1{if(!seen[$c]++) print $c}' "$MANIFEST")
if [ "${#WANT[@]}" -gt 0 ]; then
  labels=$(for l in $labels; do for w in "${WANT[@]}"; do [ "$l" = "$w" ] && echo "$l"; done; done)
  [ -n "$labels" ] || { echo "none of the requested labels are in the manifest" >&2; exit 2; }
fi

DRIVER="$LOGDIR/fetch_driver.log"
echo "=== xspecies fetch start=$(date -Iseconds) host=$(hostname) dry=$DRY" | tee -a "$DRIVER"
echo "    labels: $(echo $labels | tr '\n' ' ')" | tee -a "$DRIVER"

TOTAL_FAIL=0
for label in $labels; do
  OUT="$DESTROOT/$label"; mkdir -p "$OUT"
  LOG="$LOGDIR/fetch_${label}.log"
  nrun=$(awk -F'\t' -v c="$C_LABEL" -v l="$label" 'NR>1 && $c==l' "$MANIFEST" | wc -l)
  gib=$(awk -F'\t' -v c="$C_LABEL" -v l="$label" -v g="$C_GIB" 'NR>1 && $c==l{s+=$g} END{printf "%.1f",s}' "$MANIFEST")
  echo "=== $label: $nrun runs, $gib GiB  start=$(date -Iseconds)" | tee -a "$LOG" "$DRIVER"
  FAILED=()

  while IFS=$'\t' read -r run ftp md5s; do
    if [ -z "$ftp" ]; then
      echo "  [ERROR] $run has an EMPTY fastq_ftp -- submitted-format-only run, not fetched" | tee -a "$LOG"
      FAILED+=("$run:no-ftp"); continue
    fi
    IFS=';' read -ra F <<< "$ftp"
    IFS=';' read -ra M <<< "$md5s"
    if [ "${#F[@]}" -ne "${#M[@]}" ]; then
      echo "  [ERROR] $run: ${#F[@]} files but ${#M[@]} md5s -- refusing to guess the pairing" | tee -a "$LOG"
      FAILED+=("$run:md5-count"); continue
    fi
    # Redundant merged file: if any _1/_2 is listed, the bare <ACC>.fastq.gz is a
    # duplicate of both concatenated and must not be downloaded.
    has_pair=0
    for p in "${F[@]}"; do case "$(basename "$p")" in *_1.fastq.gz|*_2.fastq.gz) has_pair=1;; esac; done

    for i in "${!F[@]}"; do
      base=$(basename "${F[$i]}")
      if [ "$has_pair" = 1 ]; then
        case "$base" in
          *_1.fastq.gz|*_2.fastq.gz) : ;;
          *) echo "  [skip] $base is the redundant merged file (_1/_2 present)" | tee -a "$LOG"; continue ;;
        esac
      fi
      url="ftp://${F[$i]}"; want="${M[$i]}"; f="$OUT/$base"

      if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ]; then
        echo "  [ok] $base already verified" | tee -a "$LOG"; continue
      fi
      if [ -s "$f" ]; then
        rm -f "$f"; echo "  [refetch] $base failed md5, deleted" | tee -a "$LOG"
      fi
      if [ "$DRY" = 1 ]; then
        echo "  [dry] would fetch $base from $url" | tee -a "$LOG"; continue
      fi

      ok=0
      for attempt in 1 2 3; do
        [ "$attempt" -gt 1 ] && { echo "  [retry $attempt] $base" | tee -a "$LOG"; rm -f "$f"; sleep 15; }
        echo "  [get] $base attempt $attempt $(date +%H:%M:%S)" | tee -a "$LOG"
        wget -nv -c --tries=5 --timeout=60 --waitretry=10 -O "$f" "$url" 2>>"$LOG" || true
        got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)
        if [ "$got" = "$want" ]; then
          echo "  [md5 OK] $base" | tee -a "$LOG"; ok=1; break
        fi
        sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
        echo "  [md5 FAIL] $base attempt $attempt size=$sz want=$want got=${got:-none}" | tee -a "$LOG"
      done
      [ "$ok" = 1 ] || { echo "  [GIVE UP] $base after 3 attempts" | tee -a "$LOG"; FAILED+=("$base"); }
    done
  done < <(awk -F'\t' -v c="$C_LABEL" -v l="$label" -v r="$C_RUN" -v ftp="$C_FTP" -v m="$C_MD5" \
             'NR>1 && $c==l {print $r"\t"$ftp"\t"$m}' "$MANIFEST")

  if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "=== $label FAILED (${#FAILED[@]}): ${FAILED[*]}" | tee -a "$LOG" "$DRIVER"
    TOTAL_FAIL=$(( TOTAL_FAIL + ${#FAILED[@]} ))
  else
    echo "=== $label done=$(date -Iseconds)  all md5 verified" | tee -a "$LOG" "$DRIVER"
  fi
  du -sh "$OUT" 2>/dev/null | tee -a "$DRIVER"
  # Disk trend, so a quota problem is visible before EDQUOT silently kills a job.
  getfattr -n ceph.dir.rbytes --only-values /private/groups/carpenterlab 2>/dev/null \
    | awk -v l="$label" '{printf "%s\tafter %s\t%.3f TB used\n", strftime("%Y-%m-%dT%H:%M:%S"), l, $1/1e12}' \
    | tee -a "$NEW/logs/disk_watch.log"
done

echo "=== xspecies fetch ALL DONE=$(date -Iseconds) total_failures=$TOTAL_FAIL" | tee -a "$DRIVER"
[ "$TOTAL_FAIL" -eq 0 ] || exit 1
