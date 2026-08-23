#!/usr/bin/env bash
#SBATCH --job-name=jw_mm
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-13%7
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/jw_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/jw_%A_%a.err
set -euo pipefail
# DIAGNOSTIC (see data/mm25_diagnostic/README.md). Janich (5 RPF) + Wang (2 RPF) at BOTH mm1 and
# mm25, under ONE recipe.
#
# WHY BOTH SETTINGS. gse243134's historical mm1 BAMs came from a script whose recipe I matched
# exactly, so its mm1-vs-mm25 delta is clean. Janich's and Wang's historical scripts are both marked
# "provenance only -- LACKED --alignEndsType EndToEnd", so their archived BAMs are NOT recipe-matched
# to anything I can reproduce. Aligning each at mm1 AND mm25 here makes the within-dataset delta
# internally controlled and removes that confound entirely, at the cost of 7 extra alignments.
#
# Adapter measured, not assumed (project rule): janich 96.8% TruSeq at 101 nt, wang 92.4% at 50 nt.
# Both raw -> the same cutadapt bounds as gse243134.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
# Accessions are READ FROM DISK, not hardcoded -- a hardcoded janich list was wrong on the first
# attempt (assumed a contiguous SRR range; the real set skips 1930190-92).
declare -a SPEC
for mm in 1 25; do
  for f in $NEW/data/external/GSE72064_janich_ribo/*.fastq.gz; do
    SPEC+=("janich:$(basename "$f" .fastq.gz):$mm"); done
  for f in $NEW/data/external/wang2021_ribo_refetch/*.fastq.gz; do
    SPEC+=("wang:$(basename "$f" .fastq.gz):$mm"); done
done
IFS=: read -r DS SRR MM <<< "${SPEC[$SLURM_ARRAY_TASK_ID]}"
case $DS in
  janich) FQ=$NEW/data/external/GSE72064_janich_ribo/${SRR}.fastq.gz ;;
  wang)   FQ=$NEW/data/external/wang2021_ribo_refetch/${SRR}.fastq.gz ;;
esac
[ -s "$FQ" ] || { echo "MISSING $FQ" >&2; exit 1; }
OUTDIR=$NEW/data/mm25_diagnostic/${DS}_ribo_mm${MM}; mkdir -p "$OUTDIR"
OUTBAM=$OUTDIR/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$OUTBAM" ] && { echo "exists, skip $DS $SRR mm$MM"; exit 0; }
SCRATCH=/data/tmp/emalekos/jw_${DS}_${SRR}_mm${MM}.$$; mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT
TRIM=$SCRATCH/t.fq.gz
"$ENV/bin/cutadapt" -j 8 -a AGATCGGAAGAGC --minimum-length 20 --maximum-length 40 \
  -o "$TRIM" "$FQ" > "$OUTDIR/${SRR}.cutadapt.log" 2>&1
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$TRIM" --readFilesCommand zcat \
  --runThreadN 8 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax "$MM" --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCRATCH/${SRR}."
TBAM=$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$TBAM" ] || { echo "no transcriptome BAM for $DS $SRR mm$MM" >&2; exit 1; }
mv "$TBAM" "$OUTBAM"
cp "$SCRATCH/${SRR}.Log.final.out" "$OUTDIR/${SRR}.Log.final.out" 2>/dev/null || true
u=$(grep "Uniquely mapped reads %" "$OUTDIR/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
m=$(grep "% of reads mapped to multiple loci" "$OUTDIR/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
echo "DONE $DS $SRR mm$MM  unique=$u multi=$m"
