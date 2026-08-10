#!/usr/bin/env bash
# Adapter + insert-length audit for the three new human datasets, before any alignment.
#
# WHY THIS RUNS AT ALL. Many ENA/SRA submissions ship raw reads with the full Illumina TruSeq adapter
# still attached. Without `cutadapt -a` you get ~0% unique mapping and a 5 KB BAM stub. The check is
# cheap and the failure is expensive.
#
# ONLY ON COMPLETE FILES. Every input here is md5-verified against ENA. A byte-range slice of a
# partially-downloaded gz once reported 0% adapter where the finished file was 93%
# (feedback_adapter_check_full_file), so this decompresses from the start of a verified file rather
# than seeking into it.
#
# SANITY-CHECK AGAINST BIOLOGY, not just the percentage. Ribosome footprints are ~28-31 nt. So:
#   * a 72 nt "RPF" read with no adapter is a contradiction -- the insert cannot fill the read
#   * a read whose modal length is already ~30 nt was trimmed by the submitter
#   * an RNA-seq library has no such constraint and stays long
# That last point is the discriminator for GSE304796, where GEO labels all six samples "Ribo-seq" but
# three are said to be RNA: after adapter removal the true RPFs must collapse to ~30 nt and the RNA
# must not.
#
# Usage: adapter_check_human_2026_08.sh [dataset_dir ...]   (default: all three)
set -uo pipefail
R=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
EXT=$R/data/external
ADAPTER=AGATCGGAAGAGC          # Illumina TruSeq / small-RNA 3' adapter prefix
N_READS=${N_READS:-2000000}    # reads sampled per file
DSETS=("$@"); [[ ${#DSETS[@]} -eq 0 ]] && DSETS=(GSE208041_thp1 GSE39561_thp1 GSE304796_cart)

printf "%-20s %-14s %-6s %6s %8s %8s %9s %s\n" \
  dataset run arm rawlen adapter% modeIns medIns verdict
printf '%.0s-' {1..108}; echo
for ds in "${DSETS[@]}"; do
  m=$EXT/$ds/manifest.tsv
  [[ -s "$m" ]] || { echo "no manifest for $ds" >&2; continue; }
  while IFS=$'\t' read -r run gsm srx arm layout reads byts ftp md5s instr; do
    [[ "$run" == "run" || -z "$run" ]] && continue
    # mate 1 only: adapter read-through is a property of the insert, identical on both mates
    f=$EXT/$ds/${run}.fastq.gz; [[ -s "$f" ]] || f=$EXT/$ds/${run}_1.fastq.gz
    [[ -s "$f" ]] || { printf "%-20s %-14s %-6s  MISSING\n" "$ds" "$run" "$arm"; continue; }
    read -r rawlen pct mode med < <(
      zcat "$f" 2>/dev/null | head -n $((N_READS * 4)) | awk -v ad="$ADAPTER" '
        NR%4==2 {
          n++; L=length($0); if (L>rawmax) rawmax=L
          p = index($0, ad)
          if (p > 0) { hit++; ins = p - 1 } else { ins = L }
          c[ins]++; tot++
        }
        END {
          best=0; bl=0; for (k in c) if (c[k]+0 > best) { best=c[k]+0; bl=k }
          # median insert
          nn=0; split("", ks); i=0; for (k in c) ks[i++]=k+0
          asort(ks); run=0; medv=0
          for (j=1; j<=i; j++) { run += c[ks[j]]; if (run >= tot/2) { medv=ks[j]; break } }
          printf "%d %.1f %d %d\n", rawmax, 100*hit/n, bl, medv
        }' 2>/dev/null)
    [[ -z "${rawlen:-}" ]] && { printf "%-20s %-14s %-6s  AWK FAILED\n" "$ds" "$run" "$arm"; continue; }
    verdict="-"
    if   awk "BEGIN{exit !($pct >= 50)}"; then verdict="TRIM (-a $ADAPTER)"
    elif awk "BEGIN{exit !($pct < 5 && $mode <= 40)}"; then verdict="pre-trimmed by submitter"
    elif awk "BEGIN{exit !($pct < 5 && $mode > 40)}"; then verdict="no adapter, LONG insert -> not RPF-like"
    else verdict="partial ($pct%) -- inspect"; fi
    printf "%-20s %-14s %-6s %6s %7s%% %8s %9s %s\n" "$ds" "$run" "$arm" "$rawlen" "$pct" "$mode" "$med" "$verdict"
  done < "$m"
done
