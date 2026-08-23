#!/usr/bin/env bash
set -uo pipefail
# Re-fetch the iPSC-CM (Ruiz-Orera, PRJEB65856) FASTQs so the dataset can be rebuilt on the FINAL
# RECIPE. Task #92.
#
# WHY A RE-DOWNLOAD IS THE ONLY OPTION. The 2026-08-15 recipe audit found this dataset off-recipe:
# zero filter_tx_heldout runs, pack built 2026-07-14 (before the 2026-08-12 adoption). It cannot be
# rebuilt from anything on disk -- its pack provenance records `inputs_known: false`, the source
# psites hd5 were deleted, and the 5 aligned BAMs are gone from data/heldout_bam/ (mouse only
# remains). The FASTQs are the last recoverable point.
#
# HEAD NODE, SERIAL, ONE CONNECTION. Per the cluster rules: the head node's network is unthrottled
# while compute nodes throttle a single stream to ~0.4 MB/s and concurrent sbatch-array downloads
# trip EBI's per-IP rate limit. `wget -c` with ONE connection -- multi-connection (aria2c -x4/-x8)
# silently corrupts large files, producing the CORRECT SIZE with a bad md5.
#
# CHECKSUM, NEVER SIZE. Every file is verified against the manifest's published md5. A raced or
# truncated FASTQ routinely lands at exactly the expected length; only the checksum catches it.
# A file that fails is DELETED, not resumed -- `wget -c` onto bad bytes keeps them and appends.
#
#   nohup bash scripts/fetch_ruizorera_ipsc.sh > logs/fetch_ipsc.log 2>&1 &
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MAN=$NEW/data/external/ruizorera2024/manifest_human.tsv
OUT=$NEW/data/external/ruizorera2024/fastq          # GROUP fs: SLURM must read this later
mkdir -p "$OUT" "$NEW/logs"

echo "host=$(hostname) start=$(date -Iseconds)"
[[ "$(hostname)" == mustard* ]] || echo "WARNING: not on mustard; compute nodes throttle downloads"

# The 5 Ribo runs + the 5 matching RNA clones (BIHi242-C1..C5). Selected by ACCESSION and LABEL from
# the manifest, never by a directory glob.
mapfile -t ROWS < <(awk -F'\t' 'NR>1 && ($1 ~ /ERR125499(2[6-9]|30)/ ||
                                 ($4=="RNA-Seq" && $5 ~ /BIHi242-C[1-5]_mR/)) {print $1"\t"$4"\t"$6"\t"$7}' "$MAN")
echo "  ${#ROWS[@]} files to fetch"
(( ${#ROWS[@]} == 15 )) || { echo "ABORT: expected 15 files, got ${#ROWS[@]}" >&2; exit 1; }

ok=0; fail=0; skip=0
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r run strategy url md5 <<< "$row"
  f=$OUT/$(basename "$url")
  if [[ -s "$f" ]] && [[ "$(md5sum "$f" | cut -d' ' -f1)" == "$md5" ]]; then
    echo "  skip $(basename "$f") (present, md5 OK)"; skip=$((skip+1)); continue
  fi
  echo "  --- $run [$strategy] $(basename "$f") $(date -Iseconds)"
  # -c resumes a PARTIAL file; safe here only because any file that fails verification below is
  # deleted outright rather than resumed onto.
  if ! wget -c -q --show-progress --progress=dot:giga -O "$f" "https://$url"; then
    echo "     DOWNLOAD FAILED"; fail=$((fail+1)); continue
  fi
  got=$(md5sum "$f" | cut -d' ' -f1)
  if [[ "$got" == "$md5" ]]; then
    echo "     md5 OK $got  ($(stat -c %s "$f" | numfmt --to=iec))"; ok=$((ok+1))
  else
    echo "     md5 MISMATCH: got $got want $md5 -- DELETING (never resume onto bad bytes)"
    rm -f "$f"; fail=$((fail+1))
  fi
done

echo
echo "done=$(date -Iseconds)  verified=$ok  skipped=$skip  failed=$fail"
echo "  -> $OUT"
[[ "$fail" -eq 0 ]] || { echo "RERUN to retry the $fail failed file(s); verified ones are skipped." >&2; exit 1; }
