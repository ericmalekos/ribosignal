#!/usr/bin/env bash
#SBATCH --job-name=janich_rna_mm10
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --array=0-6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout_mm10/janich_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout_mm10/janich_%A_%a.err
set -euo pipefail
# mm10 RNA coverage for the Janich mouse-liver HELD-OUT eval set.
#
# Derived from scripts/heldout/align_janich_rna.sbatch by changing ONLY:
#   --outFilterMultimapNmax 1 -> 10, and the coverage/BAM output dirs (_mm10 suffix).
# Everything else -- especially the cutadapt recipe -- is copied verbatim on purpose. Janich RNA is
# 95% TruSeq adapter with a 34 nt median insert; without the exact poly-A + adapter trim STAR
# rejects the reads via --outFilterMatchNminOverLread and the library looks empty. Per-dataset
# trimming is the part that is easy to get wrong, so it is not re-derived.
#
# mm1 BAMs CANNOT be reused for this: an mm1 alignment has already discarded the multimapping
# reads, so mm10 coverage is only obtainable by re-aligning the FASTQs.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
SRRS=(SRR1930216 SRR1930217 SRR1930218 SRR1930221 SRR1930222 SRR1930225 SRR1930226)
SRR=${SRRS[$SLURM_ARRAY_TASK_ID]}
COVDIR=$NEW/data/heldout_rnaseq_coverage_mm10/mouse_janich_liver
OUTCOV=$COVDIR/${SRR}_coverage.hd5
mkdir -p "$COVDIR"
[ -s "$OUTCOV" ] && { echo "$SRR mm10 coverage exists"; exit 0; }
FQ=$NEW/data/external/janich_liver_rnaseq/${SRR}.fastq.gz
[ -s "$FQ" ] || { echo "MISSING fastq $FQ" >&2; exit 1; }
SCRATCH=/data/tmp/emalekos/janich_rna_mm10/$SRR; rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"
TRIMMED=$SCRATCH/${SRR}.trimmed.fq.gz
"$ENV/bin/cutadapt" -j 8 -a AGATCGGAAGAGC -a "A{18}" --times 2 --minimum-length 20 \
  -o "$TRIMMED" "$FQ" > "$SCRATCH/${SRR}.cutadapt.log" 2>&1
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$TRIMMED" --readFilesCommand zcat --runThreadN 8 \
  --outSAMtype None --quantMode TranscriptomeSAM --outFilterMultimapNmax 10 \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCRATCH/${SRR}."
TXBAM=$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam
"$ENV/bin/samtools" quickcheck "$TXBAM" || { echo "TRUNCATED $TXBAM" >&2; exit 2; }
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TXBAM" "$OUTCOV" "$SRR"
# keep the BAM (do NOT delete intermediates until the result is verified)
BAMDIR=$NEW/data/heldout_bam_mm10/mouse_janich_liver; mkdir -p "$BAMDIR"
mv "$TXBAM" "$BAMDIR/" && rm -rf "$SCRATCH"
echo "$SRR mm10 done -> $OUTCOV"
