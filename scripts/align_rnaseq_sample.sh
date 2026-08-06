#!/usr/bin/env bash
# Align one Fibroblast RNAseq PE sample to the transcriptome and write per-nt coverage.
# STAR (--quantMode TranscriptomeSAM, no genomic BAM) -> node-local scratch, then
# rnaseq_coverage.py -> per-sample coverage hd5 on the group fs, then delete the BAM.
# Resumable: skips a sample whose coverage hd5 already exists.
set -euo pipefail
SRR="$1"

ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49   # relocated out of genomes/ (2026-07)
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
FQ1=$ECH/data/rnaseq_input/${SRR}_1.fastq.gz
FQ2=$ECH/data/rnaseq_input/${SRR}_2.fastq.gz
COVDIR=$NEW/data/rnaseq_coverage/per_sample
OUTCOV=$COVDIR/${SRR}_coverage.hd5
BAMDIR=$NEW/data/rnaseq_bam_mm1          # KEEP the alignments here (user archives to long-term storage; do NOT delete)
LOGDIR=$NEW/logs/rnaseq_star
NTHREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$COVDIR" "$LOGDIR" "$BAMDIR"

if [ -s "$OUTCOV" ]; then
  echo "[$SRR] coverage exists, skipping: $OUTCOV"; exit 0
fi
for f in "$FQ1" "$FQ2"; do [ -s "$f" ] || { echo "MISSING $f" >&2; exit 1; }; done

SCRATCH=/data/tmp/emalekos/rnaseq_star/$SRR
rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"
echo "[$SRR] host=$(hostname) threads=$NTHREADS start=$(date -Iseconds)"

# mm1 (unique genomic mappers). Produce BOTH the coord-sorted GENOME BAM and the transcriptome BAM so the
# full alignments can be archived (do not delete -- user request).
"$ENV/bin/STAR" --genomeDir "$IDX" \
  --readFilesIn "$FQ1" "$FQ2" --readFilesCommand zcat \
  --runThreadN "$NTHREADS" \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 \
  --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCRATCH/${SRR}."

BAM="$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam"
GBAM="$SCRATCH/${SRR}.Aligned.sortedByCoord.out.bam"
[ -s "$BAM" ] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCRATCH/${SRR}.Log.final.out" "$LOGDIR/${SRR}.Log.final.out" 2>/dev/null || true

"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$BAM" "$OUTCOV" "$SRR"

# ARCHIVE: move both BAMs to the group fs (persistent), index the genome BAM; then clean node-local scratch.
mv "$BAM" "$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam"
[ -s "$GBAM" ] && { mv "$GBAM" "$BAMDIR/${SRR}.Aligned.sortedByCoord.out.bam"; "$ENV/bin/samtools" index -@ "$NTHREADS" "$BAMDIR/${SRR}.Aligned.sortedByCoord.out.bam" 2>/dev/null || true; }
rm -rf "$SCRATCH"
echo "[$SRR] done=$(date -Iseconds)  -> cov=$OUTCOV  bams=$BAMDIR/${SRR}.*"
