#!/usr/bin/env bash
# Verify every downloaded FASTQ against its ENA md5, and repair failures over HTTPS.
#
# WHY THIS EXISTS. The main fetch chain uses ENA's FTP endpoint and retries 3x. On
# 2026-08-24 two files exhausted all three attempts, and the failure mode was
# TRUNCATION rather than corruption: each attempt produced a DIFFERENT, SHORTER
# size, none matching.
#
#   ERR12549929  got 1,816,666,112 of 3,227,905,296  (56%)
#   SRR5667268   got 2,016,301,056 of 2,148,789,443  (94%)
#
# The same files' HTTPS endpoint returns a correct Content-Length, so FTP is the
# flaky leg. This script refetches over HTTPS instead.
#
# TWO RULES IT KEEPS:
#   * Never resume onto bad bytes. A failed file is DELETED and refetched whole;
#     `wget -c` would append to the truncated remainder and can never recover.
#   * Verify by md5, never by size. Size is checked FIRST only because it is a
#     cheap early diagnostic that names truncation explicitly; md5 is still the gate.
#
# ONE WRITER PER PATH: the label the main chain is currently downloading is SKIPPED,
# detected from logs/xspecies/fetch_driver.log. Repair finished labels only, or stop
# the chain first (`kill "$(cat /tmp/fetch_xspecies_pid.txt)"`).
#
# Usage:
#   bash scripts/xspecies/verify_and_repair.sh [--verify-only] [--fetch-missing]
#                                              [--include-active] [label ...]
#
# --fetch-missing makes this a COMPLETER as well as a verifier: it pulls anything the
# main chain never reached or gave up on. Safe to run in a loop until it exits 0.
set -uo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MANIFEST="$NEW/data/external/xspecies/download_manifest.tsv"
FQROOT="$NEW/data/external/xspecies/fastq"
DRIVER="$NEW/logs/xspecies/fetch_driver.log"
LOG="$NEW/logs/xspecies/verify_repair.log"
SRAPROV="$NEW/data/external/xspecies/sra_provenance.tsv"
ATTEMPTS=4

VERIFY_ONLY=0; INCLUDE_ACTIVE=0; FETCH_MISSING=0
while [ $# -gt 0 ]; do
  case "$1" in
    --verify-only) VERIFY_ONLY=1; shift;;
    --include-active) INCLUDE_ACTIVE=1; shift;;
    --fetch-missing) FETCH_MISSING=1; shift;;
    *) break;;
  esac
done
WANT=("$@")
[ -s "$MANIFEST" ] || { echo "missing $MANIFEST" >&2; exit 2; }

# The label the chain is inside right now: last "=== <label>: N runs" with no later
# "=== <label> done=" or "=== <label> FAILED".
ACTIVE=""
if [ -s "$DRIVER" ] && [ "$INCLUDE_ACTIVE" = 0 ]; then
  ACTIVE=$(awk '
    /^=== [a-z0-9_]+: [0-9]+ runs/ { lab=$2; sub(/:$/,"",lab); started[lab]=NR; order[NR]=lab }
    /^=== [a-z0-9_]+ (done=|FAILED)/ { lab=$2; finished[lab]=NR }
    END { best=""; bn=0
          for (l in started) if (started[l] > (l in finished ? finished[l] : -1) && started[l] > bn) { bn=started[l]; best=l }
          print best }' "$DRIVER")
  [ -n "$ACTIVE" ] && echo "chain is currently inside '$ACTIVE'; skipping it (one writer per path)"
fi

want_label() {
  [ "$1" = "$ACTIVE" ] && return 1
  [ "${#WANT[@]}" -eq 0 ] && return 0
  for w in "${WANT[@]}"; do [ "$1" = "$w" ] && return 0; done
  return 1
}

echo "=== verify_and_repair start=$(date -Iseconds) verify_only=$VERIFY_ONLY" | tee -a "$LOG"
n_ok=0; n_missing=0; n_bad=0; n_fixed=0; n_unfixed=0
declare -a UNFIXED=()

while IFS=$'\t' read -r label run ftp md5s; do
  want_label "$label" || continue
  IFS=';' read -ra F <<< "$ftp"
  IFS=';' read -ra M <<< "$md5s"
  has_pair=0
  for p in "${F[@]}"; do case "$(basename "$p")" in *_1.fastq.gz|*_2.fastq.gz) has_pair=1;; esac; done
  for i in "${!F[@]}"; do
    base=$(basename "${F[$i]}")
    if [ "$has_pair" = 1 ]; then
      case "$base" in *_1.fastq.gz|*_2.fastq.gz) : ;; *) continue ;; esac
    fi
    f="$FQROOT/$label/$base"; want="${M[$i]}"
    # A file sourced from NCBI SRA can never reproduce ENA's md5 (fasterq-dump
    # regenerates the file, so gzip framing and read order differ). Without this guard
    # this script sees the mismatch, DELETES the good SRA file and tries to refetch it
    # from the ENA endpoint that could not serve it in the first place -- which is
    # exactly what happened to SRR2087765_2 on 2026-08-26. SRA-sourced files carry
    # their own validation in sra_provenance.tsv.
    if [ -s "$SRAPROV" ] && awk -F'\t' -v b="$base" 'NR>1 && $4==b{found=1} END{exit !found}' "$SRAPROV"; then
      echo "  [skip] $base is SRA-sourced (see sra_provenance.tsv), not ENA-md5-checkable" | tee -a "$LOG"
      n_ok=$((n_ok+1)); continue
    fi
    if [ ! -e "$f" ]; then
      n_missing=$((n_missing+1))
      # --fetch-missing turns this from a verifier into a completer: it will pull
      # anything the main chain never got to, or gave up on and deleted.
      [ "$FETCH_MISSING" = 1 ] || continue
      echo "  [MISSING] $base" | tee -a "$LOG"
    else
      got=$(md5sum "$f" | cut -d' ' -f1)
      if [ "$got" = "$want" ]; then n_ok=$((n_ok+1)); continue; fi
    fi

    n_bad=$((n_bad+1))
    url="https://${F[$i]}"
    exp=$(curl -sI --max-time 120 "$url" | awk 'BEGIN{IGNORECASE=1}/^content-length:/{gsub(/\r/,"");print $2}')
    # A degraded ENA answers a file path with a ~800 byte HTML directory listing, so
    # its Content-Length is meaningless. No real FASTQ is under 100 kB; ignore it and
    # fall back to md5 as the only gate.
    if [ -n "$exp" ] && [ "$exp" -lt 100000 ] 2>/dev/null; then exp=""; fi
    cur=$(stat -c%s "$f" 2>/dev/null || echo 0)
    if [ -n "$exp" ] && [ "$cur" -lt "$exp" ]; then
      echo "  [BAD] $base TRUNCATED: $cur of $exp bytes ($(awk -v a="$cur" -v b="$exp" 'BEGIN{printf "%.0f",100*a/b}')%)" | tee -a "$LOG"
    else
      echo "  [BAD] $base md5 mismatch at size $cur (expected size ${exp:-unknown})" | tee -a "$LOG"
    fi
    [ "$VERIFY_ONLY" = 1 ] && { UNFIXED+=("$base"); n_unfixed=$((n_unfixed+1)); continue; }

    ok=0
    for a in $(seq 1 "$ATTEMPTS"); do
      # Delete first, ALWAYS: never resume onto the truncated remainder.
      rm -f "$f"
      # Alternate protocols. ENA degrades one at a time: on 2026-08-24 FTP truncated
      # while HTTPS was clean, and on 2026-08-25 HTTPS began 301-ing file paths to an
      # HTML directory listing while FTP refused connections outright.
      if [ $((a % 2)) -eq 1 ]; then proto=https; try="$url"; else proto=ftp; try="ftp://${F[$i]}"; fi
      echo "  [$proto $a/$ATTEMPTS] $base $(date +%H:%M:%S)" | tee -a "$LOG"
      wget -nv --tries=3 --timeout=120 --waitretry=20 -O "$f" "$try" 2>>"$LOG" || true
      sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
      # A 301 to a directory listing lands ~800 bytes of HTML with a .gz name.
      if [ "$sz" -gt 0 ] && [ "$sz" -lt 4096 ] && head -c 200 "$f" | grep -qi "<html\|<!DOCTYPE"; then
        echo "  [html] $base got an HTML page, not the file ($sz bytes) -- endpoint is degraded" | tee -a "$LOG"
        rm -f "$f"; sleep 30; continue
      fi
      if [ -n "$exp" ] && [ "$sz" -ne "$exp" ]; then
        echo "  [short] $base got $sz of $exp" | tee -a "$LOG"; sleep 20; continue
      fi
      got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)
      if [ "$got" = "$want" ]; then
        echo "  [md5 OK] $base (repaired over $proto)" | tee -a "$LOG"; ok=1; n_fixed=$((n_fixed+1)); break
      fi
      echo "  [md5 FAIL] $base attempt $a size=$sz want=$want got=${got:-none}" | tee -a "$LOG"; sleep 20
    done
    [ "$ok" = 1 ] || { echo "  [STILL BAD] $base after $ATTEMPTS attempts across https and ftp" | tee -a "$LOG"
                       UNFIXED+=("$base"); n_unfixed=$((n_unfixed+1)); }
  done
done < <(awk -F'\t' 'NR>1{print $1"\t"$6"\t"$17"\t"$18}' "$MANIFEST")

echo "=== verify_and_repair done=$(date -Iseconds)" | tee -a "$LOG"
echo "    verified_ok=$n_ok  not_yet_fetched=$n_missing  bad=$n_bad  repaired=$n_fixed  still_bad=$n_unfixed" | tee -a "$LOG"
if [ "$n_unfixed" -gt 0 ]; then
  echo "    STILL BAD: ${UNFIXED[*]}" | tee -a "$LOG"
  exit 1
fi
