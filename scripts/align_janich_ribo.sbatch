#!/usr/bin/env bash
#SBATCH --job-name=janich_ribo
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich_ribo_%A_%a.out
set -euo pipefail
# Re-align the five Janich mouse-liver RIBO-SEQ runs that had been deleted.
#
# RECIPE IS COPIED FROM scripts/heldout/align_gse243134.sbatch, DELIBERATELY AND EXACTLY. The point of
# rebuilding Janich is to put it in the same 3x3 factorial as GSE243134 and Wang; if it were aligned
# with a different recipe the factorial would partly measure the pipeline difference. So: same STAR
# index, same cutadapt arguments, same --outFilterMultimapNmax 1, same transcriptome output, and NO
# bowtie2 contaminant stage -- the mouse path filters ncRNA at pack time via --ncrna-tx, unlike the
# human CAR-T path which depletes up front.
#
# Trimming verified against the fresh reads rather than assumed: 101 nt raw, 96.8% carry the TruSeq
# adapter, and cutting there gives mode 30 nt with 83.1% of inserts in 26-34 nt. A poly-A probe also
# fired on 93% of reads, but cutting there instead gives mode 94 nt and 0.0% in the footprint range --
# the A-run sits inside the adapter/downstream sequence, not as a tail. TruSeq is correct here.
#
# BAMs land in an ARM-SUFFIXED directory. The reason this dataset had to be re-downloaded at all is
# that `data/heldout_bam/mouse_janich_liver/` holds the RNA BAMs under a name that says only
# "janich_liver", so the Ribo half looked archived when it was gone.
#
#   sbatch --array=0-4 scripts/align_janich_ribo.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
D=$NEW/data/external/GSE72064_janich_ribo
BAMDIR=$NEW/data/heldout_bam/mouse_janich_liver_ribo
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$BAMDIR"

mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1}' "$D/manifest.tsv")
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
FQ=$D/${SRR}.fastq.gz
OUTBAM=$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$FQ" ]] || { echo "missing $FQ" >&2; exit 1; }
[[ -s "$OUTBAM" ]] && { echo "SKIP $SRR (bam exists)"; exit 0; }
SCR=/data/tmp/emalekos/janich_ribo_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "run=$SRR host=$(hostname) start=$(date -Iseconds)"

"$ENV/bin/cutadapt" -j "$T" -a AGATCGGAAGAGC --minimum-length 20 --maximum-length 40 \
  -o "$SCR/${SRR}.trimmed.fq.gz" "$FQ" > "$BAMDIR/${SRR}.cutadapt.log" 2>&1
grep -E "Total reads processed|Reads written" "$BAMDIR/${SRR}.cutadapt.log" | sed 's/^/  /'

"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/${SRR}.trimmed.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/${SRR}."
cp "$SCR/${SRR}.Log.final.out" "$BAMDIR/${SRR}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$BAMDIR/${SRR}.STAR.Log.final.out" | sed 's/^/  /'

TBAM=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$TBAM" ]] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$TBAM" "$OUTBAM"
echo "done=$(date -Iseconds) -> $OUTBAM"
