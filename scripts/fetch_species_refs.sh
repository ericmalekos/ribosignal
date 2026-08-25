#!/usr/bin/env bash
# SUPERSEDED 2026-08-24. The cross-species reference tree is now built and owned by
#   scripts/xspecies/build_species_refs.sbatch  ->  genomes/xspecies_refs/<species>/
# (7 species incl. fly). This script wrote to genomes/primates_t2t/, which no longer exists.
# Kept for provenance only -- do NOT run it; it would recreate a stale parallel tree.
# Fetch genome FASTA + GTF for the cross-species expansion, into the SHARED reference tree
# (standing rule: one canonical copy per (tool, species, annotation version), reused across projects).
#
# HEAD NODE, SERIAL, ONE CONNECTION. NCBI datasets zips carry an md5 manifest; we verify it.
# Primates use the T2T assemblies (NCBI), NOT Ensembl 116, which still ships gorGor4 / Pan_tro_3.0.
#
#   fetch_species_refs.sh <label> <GCF accession>
set -euo pipefail
G=/private/groups/carpenterlab/emalekos/genomes/xspecies_refs
LAB=${1:?label}; ACC=${2:?GCF accession}
OUT=$G/$LAB; mkdir -p "$OUT"
LOG=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/fetch_ref_${LAB}.log
Z=$OUT/${ACC}.zip
echo "=== $LAB ($ACC) start=$(date -Iseconds)" | tee -a "$LOG"
if [ -s "$OUT/genome.fna" ] && [ -s "$OUT/genomic.gtf" ]; then
  echo "  already extracted, skipping" | tee -a "$LOG"; exit 0
fi
URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${ACC}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GTF&filename=${ACC}.zip"
wget -nv -c --tries=5 --timeout=120 --waitretry=15 -O "$Z" "$URL" 2>>"$LOG"
unzip -o -q "$Z" -d "$OUT/_x" 2>>"$LOG"
# NCBI zips carry md5checksums.txt -- verify before trusting
if [ -s "$OUT/_x/md5sum.txt" ]; then
  (cd "$OUT/_x" && md5sum -c md5sum.txt --quiet 2>>"$LOG") && echo "  md5 OK" | tee -a "$LOG" \
    || { echo "  MD5 FAIL" | tee -a "$LOG"; exit 1; }
fi
fna=$(find "$OUT/_x" -name "*_genomic.fna" | head -1)
gtf=$(find "$OUT/_x" -name "genomic.gtf" -o -name "*.gtf" | head -1)
[ -s "$fna" ] || { echo "  no genome fasta in zip" | tee -a "$LOG"; exit 1; }
[ -s "$gtf" ] || { echo "  no gtf in zip" | tee -a "$LOG"; exit 1; }
mv "$fna" "$OUT/genome.fna"; mv "$gtf" "$OUT/genomic.gtf"
rm -rf "$OUT/_x" "$Z"
printf "  genome %s\n  gtf    %s\n" "$(du -h "$OUT/genome.fna"|cut -f1)" "$(du -h "$OUT/genomic.gtf"|cut -f1)" | tee -a "$LOG"
echo "{\"species\":\"$LAB\",\"accession\":\"$ACC\",\"fetched\":\"$(date -Iseconds)\",\"source\":\"NCBI datasets v2alpha\"}" > "$OUT/PROVENANCE.json"
echo "=== $LAB done=$(date -Iseconds)" | tee -a "$LOG"
