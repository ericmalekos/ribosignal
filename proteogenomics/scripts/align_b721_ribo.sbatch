#!/usr/bin/env bash
# B721.221 Ribo-seq alignment (Ouspenskaia 2022, GSE143263 / PRJNA599422). 7 runs, 4 HLA alleles.
#
# LIBRARY PREP IS DIFFERENT FROM HBL-1 -- do not copy align_hbl1_ribo.sbatch.
# GEO protocol (GSM4255693) states: Illumina TruSeq Ribo Profile kit (RPHMR12126), Ribo-Zero rRNA
# depletion, gel-cut 28-32 nt, NextSeq 500. That kit ligates the 3' adapter DIRECTLY and uses NO
# UMIs. HBL-1 (Ruiz-Cuevas) uses the McGlincy-Ingolia prep with a 2 nt 5' UMI and 5 nt 3' UMI, so
# its `-u 2` and `-a NNNNNAGATCGGAAGAGC` would strip 2 real bases off every B721 read and shift
# every P-site by 2 nt. Verified from the protocol, not assumed.
#
# Adapter check on the full downloaded files: 96-97% of reads carry AGATCGGAAGAGC, median raw read
# 50-51 nt. Trimming is mandatory; without it STAR gives ~0% unique and a stub BAM.
#
# Footprint window 26-34 nt brackets the 28-32 nt gel cut with a little slack.
#
# STANDING RIBO-SEQ RULES (CLAUDE.md), both applied here:
#   --outFilterMultimapNmax 1  -- RiboCode's process_bam.py does NOT drop multi-mappers, so they
#                                 silently inflate P-site counts and break periodicity.
#   ncRNA-locus drop           -- rRNA/tRNA/miRNA/Mt_rRNA/Mt_tRNA removed before the caller.
#
#   sbatch --array=0-6 proteogenomics/scripts/align_b721_ribo.sbatch
#SBATCH --job-name=align_b721_ribo
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_b721_ribo_%A_%a.out
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
export PATH="${ENV}/bin:${PATH}"
STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
NCRNA=$ECH/data/phase30/ncrna_filter_tx.txt
TX2GENE=$ECH/data/phase30/tx_to_gene.tsv

SRRS=(SRR10846525 SRR10846526 SRR10846527 SRR10846528 SRR10846529 SRR10846530 SRR10846531)
ALLELE=(A0101 A0101 A3303 B1501 B1501 B4402 B4402)
SRR=${SRRS[$SLURM_ARRAY_TASK_ID]}; AL=${ALLELE[$SLURM_ARRAY_TASK_ID]}
RB=$NEW/proteogenomics/data/B721_pilot/riboseq
FASTQ_IN=$RB/$SRR.fastq.gz
OUT_TX=$RB/$SRR.Aligned.toTranscriptome.out.bam
[[ -f "$FASTQ_IN" ]] || { echo "missing $FASTQ_IN" >&2; exit 1; }
if [[ "${FORCE:-0}" != "1" && -s "$OUT_TX" ]]; then echo "already done: $SRR"; exit 0; fi

THREADS=${SLURM_CPUS_PER_TASK:-16}
# Node-local scratch is correct here: THIS job both writes and reads it, and the deliverable BAM is
# copied out to the group fs before the trap fires.
SCR=/data/tmp/emalekos/b721_align/${SRR}_${SLURM_JOB_ID}; mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT
echo "[$(date +%T)] $SRR ($AL) host=$(hostname) TruSeq Ribo Profile, NO UMI"

cutadapt -j "$THREADS" -a AGATCGGAAGAGC --trim-n --minimum-length 26 --maximum-length 34 \
    -o "$SCR/$SRR.trim.fastq.gz" "$FASTQ_IN" > "$RB/$SRR.cutadapt.log" 2>&1
grep -E "Reads written|too short|too long" "$RB/$SRR.cutadapt.log" | head -3

STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" \
    --readFilesIn "$SCR/$SRR.trim.fastq.gz" --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
    --alignEndsType EndToEnd --seedSearchStartLmax 15 \
    --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 25 \
    --outFilterMultimapNmax 1 \
    --outSAMattributes NH HI AS nM MD --limitBAMsortRAM 45000000000 \
    --outFileNamePrefix "$SCR/$SRR." --outTmpDir "$SCR/_STARtmp" \
    > "$RB/$SRR.star.log" 2>&1
cp "$SCR/$SRR.Log.final.out" "$RB/$SRR.Log.final.out"
UNIQ=$(grep "Uniquely mapped reads %" "$RB/$SRR.Log.final.out" | grep -oE "[0-9.]+%")
echo "  unique mapping: $UNIQ"

# Move the STAR BAM to the GROUP FS **before** filtering. cutadapt + STAR is the expensive part
# (~8 min/run); the filter is seconds. Filtering inside node-local scratch meant a filter error
# tripped `set -e`, the EXIT trap wiped $SCR, and the alignment was lost -- which is exactly what
# happened when this script was first run with the wrong filter CLI.
mv "$SCR/$SRR.Aligned.toTranscriptome.out.bam" "$OUT_TX"

# Drop ncRNA loci before the caller sees the BAM.
# filter_tx_heldout.py takes POSITIONAL args and rewrites the BAM IN PLACE (sorts + indexes):
#   filter_tx_heldout.py <bam> <ncrna_tx.txt> <tx_to_gene.tsv>
if [[ -s "$NCRNA" && -s "$TX2GENE" ]]; then
  python3 "$NEW/scripts/heldout/filter_tx_heldout.py" "$OUT_TX" "$NCRNA" "$TX2GENE" 2>&1 | tail -3
else
  echo "  WARNING: missing $NCRNA or $TX2GENE -- BAM NOT ncRNA-filtered" >&2
  samtools index "$OUT_TX"
fi
echo "[$(date +%T)] $SRR done: $(du -h "$OUT_TX" | cut -f1)  unique=$UNIQ"
