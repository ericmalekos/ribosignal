#!/usr/bin/env bash
#SBATCH --job-name=mm25_aln
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-20%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/aln_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/aln_%A_%a.err
set -euo pipefail
# ============================================================================================
# DIAGNOSTIC ONLY -- DELIBERATELY VIOLATES THE PROJECT'S RIBO-SEQ RULE.
#
# The standing rule is --outFilterMultimapNmax 1 for every Ribo-seq alignment, because RiboCode
# does NOT drop multimappers itself and they inflate P-site counts and corrupt periodicity.
# This job uses **25** on purpose, to test one question:
#
#   are the ORFs the model calls that real Ribo-seq misses actually HIDDEN by multimap filtering?
#
# If those calls appear when multimappers are kept, the model is recovering real translation the
# standard pipeline discards. If they do not, they are model false positives.
#
# Output goes to data/mm25_diagnostic/ and MUST NOT be fed into any pack, training run, or
# headline ORF-call table. Everything else -- cutadapt bounds, index, EndToEnd-equivalent
# settings, --quantMode -- is byte-identical to scripts/heldout/align_gse243134.sbatch so the
# ONLY variable is the multimap cap.
# ============================================================================================
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
D=$NEW/data/external/gse243134_mouse_liver
OUTDIR=$NEW/data/mm25_diagnostic/gse243134_ribo_mm25; mkdir -p "$OUTDIR"
mapfile -t RUNS < $NEW/data/liver3x3/gse243134_ribo_runs.txt
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
OUTBAM=$OUTDIR/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$OUTBAM" ] && { echo "exists, skip $SRR"; exit 0; }
FQ=$D/fastq/${SRR}.fastq.gz
[ -s "$FQ" ] || { echo "MISSING $FQ" >&2; exit 1; }
SCRATCH=/data/tmp/emalekos/mm25_$SRR.$$; mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT
TRIM=$SCRATCH/${SRR}.trim.fq.gz
# RPF bounds, identical to the mm1 recipe: 93.1% adapter, median insert 29 nt
"$ENV/bin/cutadapt" -j 8 -a AGATCGGAAGAGC --minimum-length 20 --maximum-length 40 \
  -o "$TRIM" "$FQ" > "$OUTDIR/${SRR}.cutadapt.log" 2>&1
grep -E "Reads written|Total reads processed" "$OUTDIR/${SRR}.cutadapt.log" | sed 's/^/  /'
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$TRIM" --readFilesCommand zcat \
  --runThreadN 8 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 25 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCRATCH/${SRR}."
TBAM=$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM for $SRR" >&2; exit 1; }
mv "$TBAM" "$OUTBAM"
cp "$SCRATCH/${SRR}.Log.final.out" "$OUTDIR/${SRR}.Log.final.out" 2>/dev/null || true
u=$(grep "Uniquely mapped reads %" "$OUTDIR/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
m=$(grep "% of reads mapped to multiple loci" "$OUTDIR/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
echo "DONE $SRR  unique=$u  multi=$m  $(date -Iseconds)"
