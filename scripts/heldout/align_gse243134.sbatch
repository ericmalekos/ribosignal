#!/usr/bin/env bash
#SBATCH --job-name=aln_243134
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-39%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/gse243134/aln_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/gse243134/aln_%A_%a.err
set -euo pipefail
# NB: the --output/--error directory MUST exist before submission -- SLURM opens those files
# before the script body runs, so an in-script `mkdir -p` is too late (all 40 tasks failed in
# 1s with exit 1:0 the first time). Create logs/gse243134/ at submit time.
# GSE243134 mouse liver (third independent liver arm): trim -> STAR -> per-run coverage (RNA) or
# transcriptome BAM (RPF). One run per array task; the _340_/_347_ flowcell pairs are pooled later at
# pack-build time (prepare_pack pools every hd5/psites file in a directory).
#
# BOTH assays are untrimmed and need cutadapt, but with DIFFERENT length bounds -- measured on the
# downloaded files, not a stream (see the dataset README):
#   RPF      93.1% adapter, median insert 29 nt, sharp 28-30 peak -> --minimum-length 20 --maximum-length 40
#   totalRNA 96.9% adapter, median insert 33 nt, flat 23-36       -> --minimum-length 20 (no cap)
# The footprint cap is the project Ribo-seq rule; applying it to RNA would silently truncate the library.
#
# STAR uses --outFilterMultimapNmax 1 for BOTH (project rule: Ribo-seq keeps unique mappers only, and the
# RNA coverage input is built the same way as every other dataset in this project).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
D=$NEW/data/external/gse243134_mouse_liver
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
mkdir -p "$NEW/logs/gse243134"

# run + assay from the manifest (skip header)
ROW=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+2{print $1"\t"$2}' "$D/manifest_wt.tsv")
SRR=$(echo "$ROW" | cut -f1); SAMP=$(echo "$ROW" | cut -f2)
[ -n "$SRR" ] || { echo "no run at index $SLURM_ARRAY_TASK_ID"; exit 0; }
case "$SAMP" in *RPF*) ASSAY=ribo ;; *) ASSAY=rna ;; esac

COVDIR=$NEW/data/heldout_rnaseq_coverage/mouse_gse243134_liver
BAMDIR=$NEW/data/heldout_bam/mouse_gse243134_liver
TRIMDIR=$D/trimmed
mkdir -p "$COVDIR" "$BAMDIR" "$TRIMDIR"
echo "run=$SRR sample=$SAMP assay=$ASSAY host=$(hostname) start=$(date -Iseconds)"

OUTBAM=$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam
OUTCOV=$COVDIR/${SRR}_coverage.hd5
if [ "$ASSAY" = rna ] && [ -s "$OUTCOV" ]; then echo "SKIP $SRR (coverage exists)"; exit 0; fi
if [ "$ASSAY" = ribo ] && [ -s "$OUTBAM" ]; then echo "SKIP $SRR (bam exists)"; exit 0; fi

FQ=$D/fastq/${SRR}.fastq.gz
[ -s "$FQ" ] || { echo "MISSING $FQ"; exit 1; }
SCRATCH=/data/tmp/emalekos/gse243134/$SRR; rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"
TRIM=$SCRATCH/${SRR}.trimmed.fq.gz

if [ "$ASSAY" = ribo ]; then
  "$ENV/bin/cutadapt" -j 8 -a AGATCGGAAGAGC --minimum-length 20 --maximum-length 40 \
    -o "$TRIM" "$FQ" > "$TRIMDIR/${SRR}.cutadapt.log" 2>&1
else
  "$ENV/bin/cutadapt" -j 8 -a AGATCGGAAGAGC -a "A{18}" --times 2 --minimum-length 20 \
    -o "$TRIM" "$FQ" > "$TRIMDIR/${SRR}.cutadapt.log" 2>&1
fi
grep -E "Reads written|Total reads processed" "$TRIMDIR/${SRR}.cutadapt.log" | sed 's/^/  /'

"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$TRIM" --readFilesCommand zcat \
  --runThreadN 8 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCRATCH/${SRR}."
grep -E "Uniquely mapped reads %|input reads" "$SCRATCH/${SRR}.Log.final.out" 2>/dev/null | sed 's/^/  /'

TBAM=$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM for $SRR"; exit 1; }
if [ "$ASSAY" = rna ]; then
  "$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$OUTCOV" "$SRR"
  echo "  coverage -> $OUTCOV"
fi
mv "$TBAM" "$OUTBAM"
cp "$SCRATCH/${SRR}.Log.final.out" "$BAMDIR/${SRR}.Log.final.out" 2>/dev/null || true
rm -rf "$SCRATCH"
echo "DONE $SRR ($ASSAY) $(date -Iseconds)"
