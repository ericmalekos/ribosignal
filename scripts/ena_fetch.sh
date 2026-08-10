#!/usr/bin/env bash
# Resilient serial ENA fetcher. One place, so the next dataset does not get a fourth ad-hoc copy.
#
# WHAT WENT WRONG THE NAIVE WAY. Fetching 74 Chothani runs back-to-back over FTP gave 13 successes
# then 61 straight failures -- ENA's per-IP limit, not corruption (all 61 were connection failures
# with zero-byte files, no md5 mismatches). A tight serial loop is enough to trip it.
#
# WHAT THIS DOES DIFFERENTLY
#   * HTTPS (ftp.sra.ebi.ac.uk over https) instead of FTP -- far less aggressively throttled.
#   * A pause between files, and a longer back-off after a failure.
#   * SEVERAL PASSES: each pass retries only what is still missing, so a transient block costs a
#     retry rather than the run. Stops when a pass makes no progress.
#   * md5 against the manifest after every file. A file that fails is DELETED, never resumed --
#     `wget -c` onto corrupt bytes keeps them and appends (feedback: an ENA file once arrived at
#     exactly the right size with the wrong md5).
#
# Manifest: TSV with a header, and columns named `ftp` and `md5` (';'-separated for paired runs).
# Any other columns are ignored, so the same manifests used elsewhere work unchanged.
#
# Usage: ena_fetch.sh <dir-with-manifest.tsv> [passes] [sleep_seconds]
set -uo pipefail
D=${1:?usage: ena_fetch.sh <dir> [passes] [sleep]}
PASSES=${2:-4}
NAP=${3:-3}
M=$D/manifest.tsv
[[ -s "$M" ]] || { echo "no manifest at $M" >&2; exit 1; }
cd "$D" || exit 1

col() { head -1 "$M" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
FC=$(col ftp); MC=$(col md5)
[[ -n "$FC" && -n "$MC" ]] || { echo "manifest needs 'ftp' and 'md5' columns" >&2; exit 1; }

verify() {  # file expected_md5 -> 0 if good
  [[ -s "$1" ]] || return 1
  [[ "$(md5sum "$1" | cut -d' ' -f1)" == "$2" ]]
}

for pass in $(seq 1 "$PASSES"); do
  got=0; left=0
  while IFS=$'\t' read -r -a F; do
    [[ "${F[0]}" == "run" || -z "${F[0]:-}" ]] && continue
    IFS=';' read -ra URLS <<< "${F[$((FC-1))]}"
    IFS=';' read -ra MD5S <<< "${F[$((MC-1))]}"
    for i in "${!URLS[@]}"; do
      raw=${URLS[$i]}; want=${MD5S[$i]}; f=$(basename "$raw")
      verify "$f" "$want" && continue
      rm -f "$f"                                  # never resume onto unverified bytes
      url="https://${raw#ftp://}"                 # HTTPS, not FTP
      if wget -q --tries=3 --timeout=90 -O "$f" "$url" && verify "$f" "$want"; then
        echo "OK   pass$pass $f"; got=$((got+1)); sleep "$NAP"
      else
        rm -f "$f"; echo "MISS pass$pass $f"; left=$((left+1)); sleep $((NAP*4))
      fi
    done
  done < "$M"
  echo "=== pass $pass: $got fetched, $left still missing ==="
  [[ "$left" -eq 0 ]] && { echo "=== complete ==="; exit 0; }
  [[ "$got" -eq 0 ]] && { echo "=== no progress this pass; stopping ==="; exit 1; }
done
echo "=== exhausted $PASSES passes ==="; exit 1
