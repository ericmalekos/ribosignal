#!/usr/bin/env bash
# Fetch held-out RNA-seq FASTQs from ENA so they can be RE-ALIGNED at mm10.
#
# Why re-download at all: mm1 BAMs have already DISCARDED the multimapping reads, so mm10 coverage
# cannot be recovered from them. Only the raw FASTQ works. (And reverting a filtered BAM to FASTQ
# is not raw FASTQ.)
#
# HEAD NODE, SERIAL, SINGLE CONNECTION, md5-verified -- per the standing rules. Compute nodes
# throttle a single stream to ~0.4 MB/s and parallel array downloads trip ENA's per-IP limit;
# multi-connection (aria2c -x4) silently corrupts large files at exactly the right size.
#
# Resume policy: a file already matching its ENA md5 is skipped. A file that FAILS md5 is DELETED
# and refetched from scratch -- `wget -c` would keep the bad bytes and append to them.
#
#   fetch_heldout_rna.sh <dataset> <acc> [acc ...]
# Writes to data/external/heldout_rna_mm10/<dataset>/ and logs to logs/fetch_<dataset>.log
set -euo pipefail
R=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
DS=${1:?dataset name}; shift
ACCS=("$@")
[ "${#ACCS[@]}" -gt 0 ] || { echo "no accessions given" >&2; exit 1; }
OUT=$R/data/external/heldout_rna_mm10/$DS; mkdir -p "$OUT"
LOG=$R/logs/fetch_${DS}.log; mkdir -p "$R/logs"
FAILED=()
echo "=== $DS: ${#ACCS[@]} accessions  start=$(date -Iseconds)" | tee -a "$LOG"
for acc in "${ACCS[@]}"; do
  meta=$(curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" | tail -n +2)
  [ -n "$meta" ] || { echo "  $acc: NO ENA RECORD" | tee -a "$LOG"; continue; }
  IFS=$'\t' read -r run ftp md5s bytes <<< "$meta"
  IFS=';' read -ra F <<< "$ftp"; IFS=';' read -ra M <<< "$md5s"
  for i in "${!F[@]}"; do
    # HTTPS, not FTP. ENA serves both from the same path, but this cluster's FTP data channel
    # fails (wget reports the file then transfers 0 bytes at 0.00 B/s -- passive-mode data ports are
    # filtered; the same filtering blocked MassIVE port 21). HTTPS returns 200 OK on the identical URL.
    url="https://${F[$i]}"; want="${M[$i]}"; f="$OUT/$(basename "${F[$i]}")"
    if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ]; then
      echo "  [ok] $(basename "$f") already verified" | tee -a "$LOG"; continue
    fi
    if [ -s "$f" ]; then
      rm -f "$f"; echo "  [refetch] $(basename "$f") failed md5, deleted" | tee -a "$LOG"
    fi
    # -nv NOT -q: -q suppresses wget's ERRORS too, so a truncated transfer left no diagnostic at
    # all (observed: SRR5262874 stopped at 13 MB of 1,594 MB with an empty log).
    # Retry from scratch on md5 failure -- `wget -c` would resume onto the bad bytes.
    ok=0
    for attempt in 1 2 3; do
      [ "$attempt" -gt 1 ] && { echo "  [retry $attempt] $(basename "$f")" | tee -a "$LOG"; rm -f "$f"; sleep 15; }
      echo "  [get] $(basename "$f") attempt $attempt" | tee -a "$LOG"
      wget -nv -c --tries=5 --timeout=60 --waitretry=10 -O "$f" "$url" 2>>"$LOG" || true
      got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)
      if [ "$got" = "$want" ]; then
        echo "  [md5 OK] $(basename "$f")" | tee -a "$LOG"; ok=1; break
      fi
      sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
      echo "  [md5 FAIL] $(basename "$f") attempt $attempt size=$sz want=$want got=${got:-none}" | tee -a "$LOG"
    done
    [ "$ok" = 1 ] || { echo "  [GIVE UP] $(basename "$f") after 3 attempts" | tee -a "$LOG"; FAILED+=("$(basename "$f")"); }
  done
done
if [ "${#FAILED[@]}" -gt 0 ]; then
  echo "=== $DS FAILED (${#FAILED[@]}): ${FAILED[*]}" | tee -a "$LOG"
  echo "=== $DS done-with-failures=$(date -Iseconds)" | tee -a "$LOG"; exit 1
fi
echo "=== $DS done=$(date -Iseconds)  all md5 verified" | tee -a "$LOG"
