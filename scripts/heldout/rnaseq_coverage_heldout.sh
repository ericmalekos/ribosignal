#!/usr/bin/env bash
# Align one held-out RNAseq sample (PE or SE, autodetected) to the transcriptome and write per-nt
# coverage hd5, then delete the BAM. Mirrors the training-data method (align_rnaseq_sample.sh):
# STAR --quantMode TranscriptomeSAM (no genomic BAM) -> node-local scratch -> rnaseq_coverage.py
# (vlen coverage hd5, libtype ISR sense-counted, per-nt length == tx length, same axis as the
# P-site target). Parameterized by DATASET (species -> STAR index) via heldout_config.sh. Resumable.
set -euo pipefail
DATASET="$1"; SRR="$2"
HELDOUT_DIR=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/heldout
export DATASET
# shellcheck source=/dev/null
source "$HELDOUT_DIR/heldout_config.sh"

COVDIR="$NEW/data/heldout_rnaseq_coverage/$DATASET"
OUTCOV="$COVDIR/${SRR}_coverage.hd5"
mkdir -p "$COVDIR" "$LOG_DIR"
if [ -s "$OUTCOV" ]; then echo "[$SRR] coverage exists, skipping: $OUTCOV"; exit 0; fi

# autodetect PE (_1/_2) vs SE (single) in the dataset's fastq dir
FQ1="$FASTQ_DIR/${SRR}_1.fastq.gz"; FQ2="$FASTQ_DIR/${SRR}_2.fastq.gz"; FQ0="$FASTQ_DIR/${SRR}.fastq.gz"
if [ -s "$FQ1" ] && [ -s "$FQ2" ]; then READS="$FQ1 $FQ2"; LAYOUT=PE
elif [ -s "$FQ0" ]; then READS="$FQ0"; LAYOUT=SE
else echo "MISSING fastq for $SRR in $FASTQ_DIR" >&2; exit 1; fi

SCRATCH=/data/tmp/emalekos/heldout_rnaseq/$DATASET/$SRR
rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"
NTHREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[$SRR] dataset=$DATASET layout=$LAYOUT idx=$(basename "$STAR_INDEX") host=$(hostname) threads=$NTHREADS start=$(date -Iseconds)"

# shellcheck disable=SC2086
"$RIBOSEQ_ENV/bin/STAR" --genomeDir "$STAR_INDEX" \
  --readFilesIn $READS --readFilesCommand zcat \
  --runThreadN "$NTHREADS" \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 \
  --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCRATCH/${SRR}."

BAM="$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam"
[ -s "$BAM" ] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCRATCH/${SRR}.Log.final.out" "$LOG_DIR/rnaseq_${SRR}.Log.final.out" 2>/dev/null || true

"$RIBOSEQ_ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$BAM" "$OUTCOV" "$SRR"
# ARCHIVE: keep both BAMs (user request -> long-term storage); do not delete.
BAMDIR="$NEW/data/heldout_bam/$DATASET"; mkdir -p "$BAMDIR"
GBAM="$SCRATCH/${SRR}.Aligned.sortedByCoord.out.bam"
mv "$BAM" "$BAMDIR/${SRR}.rnaseq.Aligned.toTranscriptome.out.bam"
[ -s "$GBAM" ] && { mv "$GBAM" "$BAMDIR/${SRR}.rnaseq.Aligned.sortedByCoord.out.bam"; "$RIBOSEQ_ENV/bin/samtools" index -@ "${SLURM_CPUS_PER_TASK:-8}" "$BAMDIR/${SRR}.rnaseq.Aligned.sortedByCoord.out.bam" 2>/dev/null || true; }
rm -rf "$SCRATCH"
echo "[$SRR] done=$(date -Iseconds)  -> cov=$OUTCOV  bams=$BAMDIR/${SRR}.rnaseq.*"
