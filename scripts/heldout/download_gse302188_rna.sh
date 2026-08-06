#!/usr/bin/env bash
set -uo pipefail
# GSE302188 mouse liver poly(A) RNA-seq -- replacement RNA arm for the third liver dataset.
# Replaces GSE243134's rejected Ribo-Zero totalRNA (17.7% salmon, 93.6% of TPM in 10 structural
# transcripts, universe 3,321 tx). See ../../data/external/mouse_liver_polya_rnaseq/README.md.
#
# Four "Liver, 3 Months" biological replicates (Rep 9/14/16/10), 219.4M pairs, 22.4 GB, PE 2x101.
# poly(A) evidence is GEO protocol text (`!Sample_molecule_ch1 = polyA RNA`, "TruSeq RNA stranded kit
# ... with polyA enrichment"), NOT SRA library_selection -- that field reads `cDNA` here and also reads
# `cDNA` for Ribo-Zero libraries, so it proves nothing.
#
# Downloads land on the GROUP FS, not /data/tmp: head-node /data/tmp is node-local and invisible to
# SLURM compute nodes. Serial, single-stream, md5-verified (multi-connection silently corrupts large
# files at the correct size); ENA resets connections under load, so each file retries with backoff.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
D=$NEW/data/external/mouse_liver_polya_rnaseq
FQ=$D/fastq
mkdir -p "$FQ" "$NEW/logs/liver_polya"

# SRR34449471 first: it is the screen run that gates the other three.
RUNS="SRR34449471 SRR34449553 SRR34449454 SRR34449569"

# Adopt any already-verified file sitting in the old screen scratch rather than refetching it.
for f in /data/tmp/emalekos/liver_polya_screen/fastq/SRR*.fastq.gz; do
  [ -s "$f" ] && [ ! -s "$FQ/$(basename "$f")" ] && { echo "adopting $(basename "$f") from scratch"; mv "$f" "$FQ/"; }
done

MAN=$D/manifest.tsv
[ -s "$MAN" ] || printf "run\tsample_title\tread_count\tfile\tbytes\tmd5\tstatus\n" > "$MAN"

for RUN in $RUNS; do
  TSV=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=sample_title,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv")
  # ENA ALWAYS prepends run_accession as column 1 regardless of the `fields` list -- resolve by header
  # name, never by position (positional cut wrote 0-byte files named after the accession).
  col() { echo "$TSV" | awk -F'\t' -v k="$1" 'NR==1{for(i=1;i<=NF;i++)if($i==k)c=i;next}{print $c}'; }
  TITLE=$(col sample_title); NREAD=$(col read_count)
  URLS=$(col fastq_ftp | tr ';' ' '); MD5S=$(col fastq_md5 | tr ';' ' '); SIZES=$(col fastq_bytes | tr ';' ' ')
  [ -n "$URLS" ] || { echo "NO ENA fastq_ftp for $RUN -- skipping"; continue; }
  echo "=== $RUN  $TITLE  ${NREAD} reads"

  i=1
  for U in $URLS; do
    F="$FQ/$(basename "$U")"
    WANT=$(echo "$MD5S" | cut -d' ' -f$i); SZ=$(echo "$SIZES" | cut -d' ' -f$i); i=$((i+1))
    if [ -s "$F" ]; then
      GOT=$(md5sum "$F" | cut -d' ' -f1)
      if [ "$GOT" = "$WANT" ]; then
        echo "  OK (cached)  $(basename "$F")  md5 verified"
        grep -q "	$(basename "$F")	" "$MAN" || printf "%s\t%s\t%s\t%s\t%s\t%s\tverified\n" \
          "$RUN" "$TITLE" "$NREAD" "$(basename "$F")" "$SZ" "$WANT" >> "$MAN"
        continue
      fi
      echo "  cached $(basename "$F") FAILS md5 -- deleting (wget -c appends to bad bytes, it cannot repair a corrupt prefix)"
      rm -f "$F"
    fi
    OK=0
    for TRY in 1 2 3 4 5; do
      echo "  fetching $(basename "$F")  try $TRY  ($(date +%H:%M:%S)) ..."
      if wget -q --tries=3 --timeout=60 -O "$F" "https://$U"; then
        GOT=$(md5sum "$F" | cut -d' ' -f1)
        if [ "$GOT" = "$WANT" ]; then
          echo "  OK  $(basename "$F")  $(du -h "$F" | cut -f1)  md5 verified"
          printf "%s\t%s\t%s\t%s\t%s\t%s\tverified\n" "$RUN" "$TITLE" "$NREAD" "$(basename "$F")" "$SZ" "$WANT" >> "$MAN"
          OK=1; break
        fi
        echo "  md5 mismatch (want=$WANT got=$GOT) -- refetching from scratch"
      else
        echo "  transfer failed"
      fi
      rm -f "$F"; SLP=$((TRY*TRY*15)); echo "  backing off ${SLP}s"; sleep $SLP
    done
    [ "$OK" = 1 ] || { echo "  GIVING UP on $(basename "$F") after 5 tries"
      printf "%s\t%s\t%s\t%s\t%s\t%s\tFAILED\n" "$RUN" "$TITLE" "$NREAD" "$(basename "$F")" "$SZ" "$WANT" >> "$MAN"; }
  done

  # Gate: as soon as the screen run is complete, test it before spending the other ~16 GB.
  if [ "$RUN" = SRR34449471 ] && [ -s "$FQ/SRR34449471_1.fastq.gz" ] && [ -s "$FQ/SRR34449471_2.fastq.gz" ]; then
    echo "  screen run complete -> submitting salmon gate (continuing downloads in parallel)"
    sbatch "$NEW/scripts/heldout/screen_liver_polya_salmon.sbatch"
  fi
done

echo
echo "===== MANIFEST ====="
column -t -s$'\t' "$MAN"
echo
echo "verified: $(grep -c verified "$MAN")  failed: $(grep -c FAILED "$MAN" || true)"
du -sh "$FQ"
