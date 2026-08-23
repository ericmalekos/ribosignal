#!/usr/bin/env bash
# SUPERSEDED by scripts/riboseq_align.sbatch (docs/PIPELINE_POLICY.md rule 1).
# Kept for provenance only -- do NOT use for new work. Historically important: this script was
# ON-RECIPE (EndToEnd + mm1 + filter on query-grouped input), which is why the mouse leukocyte
# packs it built regenerate to 0.13-3.16%. It differs from canonical only in two tightened STAR
# filters (--outFilterMismatchNoverLmax 0.05, --outFilterMatchNminOverLread 0.7).
# Align one held-out Ribo-seq footprint FASTQ (cutadapt + STAR TranscriptomeSAM),
# then clean the transcriptome BAM in place (posture A). Usage: align_and_filter.sh <DATASET> <SRR>
set -euo pipefail
DATASET="${1:?DATASET required}"; SRR="${2:?SRR required}"
HELDOUT_DIR=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/heldout
export DATASET
source "$HELDOUT_DIR/heldout_config.sh"
export PATH="$RIBOSEQ_ENV/bin:$PATH"
mkdir -p "$BAM_DIR" "$LOG_DIR"

FASTQ_IN="$FASTQ_DIR/${SRR}.fastq.gz"
[[ -f "$FASTQ_IN" ]] || { echo "missing $FASTQ_IN" >&2; exit 1; }
SCRATCH="/data/tmp/emalekos/heldout/${DATASET}_${SRR}_$$"
mkdir -p "$SCRATCH"; trap 'rm -rf "$SCRATCH"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-16}"
TRIMMED="$SCRATCH/${SRR}.trimmed.fastq.gz"

echo "[$(date +%T)] $DATASET $SRR: cutadapt ($ADAPTER_MODE)"
if [[ "$ADAPTER_MODE" == "truseq" ]]; then
  # raw footprint reads: clip TruSeq 3' adapter, keep ribo-footprint size window
  cutadapt -j "$THREADS" -a AGATCGGAAGAGC --trim-n --minimum-length 20 --maximum-length 40 \
    -o "$TRIMMED" "$FASTQ_IN" > "$LOG_DIR/${DATASET}_${SRR}.cutadapt.log" 2>&1
else
  # trimn: already adapter-clipped + N-padded to 35 nt (Chothani / Ruiz-Orera)
  cutadapt -j "$THREADS" --trim-n --minimum-length 20 \
    -o "$TRIMMED" "$FASTQ_IN" > "$LOG_DIR/${DATASET}_${SRR}.cutadapt.log" 2>&1
fi

echo "[$(date +%T)] $DATASET $SRR: STAR (index $(basename "$STAR_INDEX"))"
STAR --runThreadN "$THREADS" \
     --genomeDir "$STAR_INDEX" \
     --readFilesIn "$TRIMMED" --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outFilterMismatchNoverLmax 0.05 \
     --outFilterMatchNminOverLread 0.7 \
     --outFilterMultimapNmax 1 \
     --alignEndsType EndToEnd \
     --outSAMattributes NH HI AS nM MD \
     --outFileNamePrefix "$SCRATCH/${SRR}." \
     --outTmpDir "$SCRATCH/_STAR_tmp" \
     --limitBAMsortRAM 30000000000 \
     > "$LOG_DIR/${DATASET}_${SRR}.star.log" 2>&1

# only the transcriptome BAM is needed downstream (RiboCode); genome BAM stays in scratch
mv "$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam" "$BAM_DIR/"
cp "$SCRATCH/${SRR}.Log.final.out" "$LOG_DIR/${DATASET}_${SRR}.Log.final.out"

echo "[$(date +%T)] $DATASET $SRR: filter transcriptome BAM (ncRNA + cross-gene)"
"$RIBOSEQ_ENV/bin/python" "$HELDOUT_DIR/filter_tx_heldout.py" \
  "$BAM_DIR/${SRR}.Aligned.toTranscriptome.out.bam" "$NCRNA_TX" "$TX_TO_GENE"

echo "[$(date +%T)] $DATASET $SRR: align+filter DONE"
grep -E "Uniquely mapped|too many loci|Number of input reads" "$LOG_DIR/${DATASET}_${SRR}.Log.final.out" || true
