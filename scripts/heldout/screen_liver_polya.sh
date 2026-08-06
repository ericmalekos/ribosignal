#!/usr/bin/env bash
set -uo pipefail
# Screen candidate poly(A) mouse-liver RNA-seq arms to replace the rejected GSE243134 totalRNA.
#
# The GSE243134 totalRNA arm failed on library chemistry, not on processing: 17.7% salmon rate,
# 93.6% of TPM in ten tRNA/7SL structural transcripts, universe 3,321 tx. Metadata alone cannot
# be trusted to avoid a repeat (SRA files Ribo-Zero libraries as library_selection=cDNA), so this
# runs the failing test directly on one run per candidate before committing to a full download.
#
# PASS BAR (vs Wang P42_Liver_mRNA poly(A): 91-93% mapping, top-10 share 32.9%, 14,331 tx TPM>=1):
#   salmon mapping rate >= 70%   AND   top-10 TPM share <= 50%   AND   TPM>=1 count >= 10,000
#
# Download runs SERIALLY on the head node with a single stream and md5 verification (project rule:
# multi-connection downloads silently corrupt large files at the correct size). salmon is sent to
# SLURM. Head-node cost here is I/O only.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
S=/data/tmp/emalekos/liver_polya_screen
mkdir -p "$S/fastq" "$NEW/logs/liver_polya"
cd "$S"

# candidate  run          label
CANDS="GSE302188:SRR34449471:C57BL/6J_3mo_PE101_NovaSeq
ENCSR000BYS:SRR5047934:C57BL/6J_8wk_PE76_ENCODE"

echo "===== DOWNLOAD (serial, single-stream, md5-verified) ====="
for C in $CANDS; do
  ACC=${C%%:*}; REST=${C#*:}; RUN=${REST%%:*}
  # NB: ENA ALWAYS prepends run_accession as column 1 regardless of the `fields` list, so positional
  # cut -f1/-f2 silently yields the accession instead of the URL (this wrote 0-byte files named after
  # the run). Resolve columns by HEADER NAME.
  TSV=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv")
  URLS=$(echo "$TSV" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="fastq_ftp")c=i;next}{print $c}' | tr ';' ' ')
  MD5S=$(echo "$TSV" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="fastq_md5")c=i;next}{print $c}' | tr ';' ' ')
  [ -n "$URLS" ] || { echo "  NO ENA fastq_ftp for $RUN ($ACC) -- skipping"; continue; }
  i=1
  for U in $URLS; do
    F="$S/fastq/$(basename "$U")"
    WANT=$(echo "$MD5S" | cut -d' ' -f$i); i=$((i+1))
    if [ -s "$F" ]; then
      GOT=$(md5sum "$F" | cut -d' ' -f1)
      [ "$GOT" = "$WANT" ] && { echo "  OK (cached)  $(basename "$F")"; continue; }
      echo "  md5 MISMATCH on cached $(basename "$F") -- removing (wget -c cannot repair a bad prefix)"
      rm -f "$F"
    fi
    # ENA resets connections under load (this project previously recovered 29/30 runs only after
    # adding backoff). Retry with exponential backoff; on an md5 mismatch DELETE and refetch from
    # scratch -- `wget -c` cannot repair a corrupted prefix, it just appends to bad bytes.
    OK=0
    for TRY in 1 2 3 4 5; do
      echo "  fetching $(basename "$F")  try $TRY  ($(date +%H:%M:%S)) ..."
      if wget -q --tries=3 --timeout=60 -O "$F" "https://$U"; then
        GOT=$(md5sum "$F" | cut -d' ' -f1)
        if [ "$GOT" = "$WANT" ]; then
          echo "  OK  $(basename "$F")  $(du -h "$F" | cut -f1)  md5 verified"; OK=1; break
        fi
        echo "  md5 mismatch (want=$WANT got=$GOT) -- refetching from scratch"
      else
        echo "  transfer failed"
      fi
      rm -f "$F"; SLP=$((TRY*TRY*15)); echo "  backing off ${SLP}s"; sleep $SLP
    done
    [ "$OK" = 1 ] || echo "  GIVING UP on $(basename "$F") after 5 tries"
  done
done

echo
echo "===== ADAPTER CHECK (full local file, not a stream slice) ====="
for F in "$S"/fastq/*_1.fastq.gz; do
  [ -s "$F" ] || continue
  N=$(zcat "$F" 2>/dev/null | head -400000 | awk 'NR%4==2' | wc -l)
  A=$(zcat "$F" 2>/dev/null | head -400000 | awk 'NR%4==2' | grep -c AGATCGGAAGAGC)
  L=$(zcat "$F" 2>/dev/null | head -4000   | awk 'NR%4==2{print length($0)}' | sort -n | awk '{v[NR]=$1}END{print v[int(NR/2)]}')
  printf "  %-28s read_len=%s  adapter=%.1f%% (%s/%s)\n" "$(basename "$F")" "$L" "$(echo "$A $N" | awk '{print ($2?100*$1/$2:0)}')" "$A" "$N"
done

echo
echo "===== submitting salmon ====="
sbatch "$NEW/scripts/heldout/screen_liver_polya_salmon.sbatch"
