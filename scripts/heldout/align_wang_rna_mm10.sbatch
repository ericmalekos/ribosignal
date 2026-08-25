#!/usr/bin/env bash
#SBATCH --job-name=wang_rna_mm10
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout_mm10/wang_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout_mm10/wang_%A_%a.err
set -euo pipefail
# mm10 RNA coverage for the Wang 2021 mouse-liver HELD-OUT eval set.
#
# Derived from scripts/align_wang_rna.sbatch by changing ONLY --outFilterMultimapNmax 1 -> 10 and the
# output dirs. The cutadapt recipe is copied verbatim, including its two deliberate omissions
# documented upstream: NO --maximum-length (RNA inserts are not ~30 nt footprints) and NO
# --discard-untrimmed (an RNA read may simply not reach the adapter; applying it to the CAR-T RNA arm
# once kept only 16.4% of reads vs 99.0% without).
#
# FASTQs are ALREADY LOCAL and md5-verified against ENA at data/external/wang2021_mouse/fastq/
# mouse_liver/ -- no download. (They were briefly re-downloaded by mistake on 2026-08-24 because an
# availability survey used `find -maxdepth 1` and missed them at depth 3; the copy was deleted.)
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
D=$NEW/data/external/wang2021_mouse/fastq/mouse_liver
T=${SLURM_CPUS_PER_TASK:-12}
RUNS=(SRR5262874 SRR5262875)
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
FQ=$D/${SRR}.fastq.gz
COVDIR=$NEW/data/heldout_rnaseq_coverage_mm10/mouse_wang_liver
BAMDIR=$NEW/data/heldout_bam_mm10/mouse_wang_liver_rna
OUTCOV=$COVDIR/${SRR}_coverage.hd5
mkdir -p "$COVDIR" "$BAMDIR"
[ -s "$FQ" ] || { echo "missing $FQ" >&2; exit 1; }
[ -s "$OUTCOV" ] && { echo "SKIP $SRR (mm10 coverage exists)"; exit 0; }
SCR=/data/tmp/emalekos/wang_rna_mm10_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"
echo "run=$SRR host=$(hostname) start=$(date -Iseconds)"
"$ENV/bin/cutadapt" -j "$T" -a AGATCGGAAGAGC --minimum-length 20 \
  -o "$SCR/${SRR}.trimmed.fq.gz" "$FQ" > "$BAMDIR/${SRR}.cutadapt.log" 2>&1
grep -E "Total reads processed|Reads written" "$BAMDIR/${SRR}.cutadapt.log" | sed 's/^/  /'
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/${SRR}.trimmed.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/${SRR}."
cp "$SCR/${SRR}.Log.final.out" "$BAMDIR/${SRR}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|% of reads mapped to multiple loci|Number of input reads" \
  "$BAMDIR/${SRR}.STAR.Log.final.out" | sed 's/^/  /'
TBAM=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
"$ENV/bin/samtools" quickcheck "$TBAM" || { echo "TRUNCATED $TBAM" >&2; exit 2; }
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$OUTCOV" "$SRR"
mv "$TBAM" "$BAMDIR/" && rm -rf "$SCR"     # keep the BAM; do not delete intermediates
echo "$SRR mm10 done -> $OUTCOV"
