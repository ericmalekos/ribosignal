#!/usr/bin/env bash
#SBATCH --job-name=hbl1_ribo_align
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_ribo_align_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_ribo_align_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
export PATH="${ENV}/bin:${PATH}"
STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
NCRNA=$ECH/data/phase30/ncrna_filter_tx.txt
TX2GENE=$ECH/data/phase30/tx_to_gene.tsv

# PAPER-MATCHED PROCESSING (Ruiz-Cuevas 2021, validated on the u2 periodicity pilot 35607559: offset 12,
# frame-0 ~80%). McGlincy-Ingolia UMI prep: [2nt 5' UMI][footprint][5nt 3' UMI/linker][TruSeq adapter].
#   trim  = cutadapt -u 2 (5' UMI) -a NNNNNAGATCGGAAGAGC (wildcard 5nt 3' UMI + adapter), footprints 26-34
#   STAR  = EndToEnd (NO soft-clip) + --seedSearchStartLmax 15 + mm0.05 + --outFilterMatchNmin 25 (paper flags)
# Outputs kept: unique genome BAM (MAPQ255) for PRICE + cleaned transcriptome BAM for RiboCode.
SRRS=(SRR12285192 SRR12285191)          # D1, D2 (DMSO elongating Ribo-seq)
SRR=${SRRS[$SLURM_ARRAY_TASK_ID]}
RB=$NEW/proteogenomics/data/HBL1_pilot/riboseq
FASTQ_IN=$RB/$SRR.fastq.gz
OUT_TX=$RB/$SRR.Aligned.toTranscriptome.out.bam     # RiboCode (cleaned)
OUT_GN=$RB/$SRR.genome.uniq.bam                      # PRICE (unique, MAPQ255)
[[ -f "$FASTQ_IN" ]] || { echo "missing $FASTQ_IN" >&2; exit 1; }
if [[ -f "$OUT_TX" && -f "$OUT_TX.bai" && -f "$OUT_GN" && -f "$OUT_GN.bai" ]]; then echo "already done: $SRR"; exit 0; fi
THREADS=${SLURM_CPUS_PER_TASK:-16}
SCR=/data/tmp/emalekos/hbl1_align/${SRR}_${SLURM_JOB_ID}; mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT

echo "[$(date +%T)] $SRR host=$(hostname): cutadapt (-u 2 5'UMI; -a NNNNNAGATCGGAAGAGC 3'UMI+adapter; 26-34)"
cutadapt -j "$THREADS" -u 2 -a NNNNNAGATCGGAAGAGC --trim-n --minimum-length 26 --maximum-length 34 \
    -o "$SCR/$SRR.trim.fastq.gz" "$FASTQ_IN" > "$RB/$SRR.cutadapt.log" 2>&1
grep -E "passing filters|too short|too long" "$RB/$SRR.cutadapt.log" | head

echo "[$(date +%T)] $SRR: STAR (paper params: EndToEnd, seedSearchStartLmax 15, mm0.05, matchNmin 25)"
STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" \
    --readFilesIn "$SCR/$SRR.trim.fastq.gz" --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
    --alignEndsType EndToEnd --seedSearchStartLmax 15 \
    --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 25 --outFilterMultimapNmax 20 \
    --outSAMattributes NH HI AS nM MD --limitBAMsortRAM 45000000000 \
    --outFileNamePrefix "$SCR/$SRR." --outTmpDir "$SCR/_STARtmp" \
    > "$RB/$SRR.star.log" 2>&1
grep -E "Uniquely mapped reads number|Uniquely mapped reads %|multiple loci" "$SCR/$SRR.Log.final.out" | head
mv "$SCR/$SRR.Log.final.out" "$RB/$SRR.Log.final.out"

echo "[$(date +%T)] $SRR: unique genome BAM (MAPQ255) for PRICE"
samtools view -@ "$THREADS" -q 255 -b "$SCR/$SRR.Aligned.sortedByCoord.out.bam" > "$OUT_GN"
samtools index -@ "$THREADS" "$OUT_GN"
echo "  unique genome reads: $(samtools view -c -@ $THREADS $OUT_GN)"

echo "[$(date +%T)] $SRR: clean transcriptome BAM for RiboCode (ncRNA + cross-gene MM drop)"
mv "$SCR/$SRR.Aligned.toTranscriptome.out.bam" "$OUT_TX"
$ENV/bin/python3 $NEW/proteogenomics/scripts/filter_tx_bam_generic.py "$OUT_TX" "$NCRNA" "$TX2GENE"
echo "[$(date +%T)] $SRR DONE: genome=$(ls -lh $OUT_GN|awk '{print $5}')  tx=$(ls -lh $OUT_TX|awk '{print $5}')"