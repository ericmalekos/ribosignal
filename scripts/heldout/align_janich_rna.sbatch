#!/usr/bin/env bash
#SBATCH --job-name=janich_rna
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=0-6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/align_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/align_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
SRRS=(SRR1930216 SRR1930217 SRR1930218 SRR1930221 SRR1930222 SRR1930225 SRR1930226)
SRR=${SRRS[$SLURM_ARRAY_TASK_ID]}
COVDIR=$NEW/data/heldout_rnaseq_coverage/mouse_janich_liver
OUTCOV=$COVDIR/${SRR}_coverage.hd5
mkdir -p $COVDIR
[ -s "$OUTCOV" ] && { echo "$SRR coverage exists"; exit 0; }
FQ=$NEW/data/external/janich_liver_rnaseq/${SRR}.fastq.gz   # CLEAN re-download (originals were corrupt)
[ -s "$FQ" ] || { echo "MISSING clean fastq $FQ (re-download not done?)" >&2; exit 1; }
SCRATCH=/data/tmp/emalekos/janich_rna/$SRR; rm -rf $SCRATCH; mkdir -p $SCRATCH
# Janich RNA is 95% TruSeq adapter, median insert 34nt (short-insert poly-A library): STAR rejects the
# untrimmed reads via --outFilterMatchNminOverLread (0.66). Trim 3' adapter + poly-A, keep reads >=20nt.
TRIMMED=$SCRATCH/${SRR}.trimmed.fq.gz
"$ENV/bin/cutadapt" -j 8 -a AGATCGGAAGAGC -a "A{18}" --times 2 --minimum-length 20 \
  -o "$TRIMMED" "$FQ" > "$SCRATCH/${SRR}.cutadapt.log" 2>&1
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$TRIMMED"   --readFilesCommand zcat --runThreadN 8 --outSAMtype BAM SortedByCoordinate   --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM   --outFileNamePrefix "$SCRATCH/${SRR}."
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam" "$OUTCOV" "$SRR"
BAMDIR=$NEW/data/heldout_bam/mouse_janich_liver; mkdir -p $BAMDIR
mv "$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam" "$BAMDIR/" 2>/dev/null || true
rm -rf $SCRATCH
echo "$SRR done -> $OUTCOV"
