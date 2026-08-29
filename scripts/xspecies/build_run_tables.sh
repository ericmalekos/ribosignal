#!/usr/bin/env bash
# Pull the FULL ENA run table for every project in the cross-species expansion.
#
# One TSV per project under data/external/xspecies/. These are the authoritative
# metadata for everything downstream: the download manifest, the adapter policy,
# and the dataset registry all derive from them.
#
# Columns are deliberately wider than the old 4-column PRJEB65856_runs.tsv, which
# had no fastq_md5 (so nothing could be verified from it) and no sample_alias for
# the non-primate projects.
#
# Usage: bash scripts/xspecies/build_run_tables.sh
# Head node, tiny requests. Safe to re-run; overwrites each TSV.
set -euo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
OUT="$NEW/data/external/xspecies"
LOG="$NEW/logs/xspecies/build_run_tables.log"
mkdir -p "$OUT" "$(dirname "$LOG")"

FIELDS="run_accession,study_accession,secondary_study_accession,sample_accession,sample_alias,sample_description,sample_title,scientific_name,library_strategy,library_selection,library_source,library_layout,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"

# project           label
PROJECTS="
PRJEB65856        ruizorera_primates
PRJNA200706       zebrafish_ribo_gse46512
PRJNA288987       zebrafish_rna_gse70549
PRJNA230441       worm_ribo_gse52905
PRJNA230374       worm_rna_gse52861
PRJNA726548       yeast_gse173654
PRJNA390134       fly_gse99920
"

echo "=== build_run_tables start=$(date -Iseconds) host=$(hostname)" | tee -a "$LOG"

fail=0
while read -r acc label; do
  [ -z "${acc:-}" ] && continue
  dst="$OUT/${acc}_runs.tsv"
  tmp="$dst.tmp"
  printf '%-14s %-28s ' "$acc" "$label" | tee -a "$LOG"
  if ! curl -sS --max-time 180 --retry 3 --retry-delay 5 -G \
        "https://www.ebi.ac.uk/ena/portal/api/filereport" \
        --data-urlencode "accession=$acc" \
        --data-urlencode "result=read_run" \
        --data-urlencode "fields=$FIELDS" \
        --data-urlencode "format=tsv" -o "$tmp"; then
    echo "FETCH FAILED" | tee -a "$LOG"; fail=1; rm -f "$tmp"; continue
  fi
  n=$(( $(wc -l < "$tmp") - 1 ))
  # A filereport error comes back as a one-line message, not a TSV with a header.
  if [ "$n" -lt 1 ] || ! head -1 "$tmp" | grep -q '^run_accession'; then
    echo "BAD RESPONSE: $(head -c 200 "$tmp")" | tee -a "$LOG"; fail=1; rm -f "$tmp"; continue
  fi
  # Every run must carry an md5 per fastq file, or nothing downstream can verify.
  nomd5=$(awk -F'\t' 'NR>1 && $17==""' "$tmp" | wc -l)
  mv "$tmp" "$dst"
  echo "${n} runs  (${nomd5} without md5)  -> ${dst##*/}" | tee -a "$LOG"
  [ "$nomd5" -gt 0 ] && { echo "  WARNING: $nomd5 runs have no fastq_md5" | tee -a "$LOG"; }
  sleep 1
done <<< "$PROJECTS"

echo "=== build_run_tables done=$(date -Iseconds) fail=$fail" | tee -a "$LOG"
exit "$fail"
