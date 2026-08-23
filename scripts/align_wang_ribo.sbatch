#!/usr/bin/env bash
#SBATCH --job-name=wang_ribo
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/wang_ribo_%A_%a.out
# SUPERSEDED by scripts/riboseq_align.sbatch (docs/PIPELINE_POLICY.md rule 1).
# Kept for provenance only -- do NOT use for new work. This script LACKED --alignEndsType EndToEnd
# and the ncRNA + cross-gene filter; its output was rebuilt canonically on 2026-08-13.
set -euo pipefail
# Align the two Wang 2021 mouse-liver RIBO-SEQ runs (SRR5262890, SRR5262891) that had been deleted.
#
# RECIPE IS COPIED FROM scripts/align_janich_ribo.sbatch, DELIBERATELY. Wang is the third arm of the
# mouse-liver 3x3 factorial; the step-0 gate showed that pooling one dataset through a different code
# path shifts P-site totals by ~1.7%, which is the same order as the effects the factorial is trying
# to measure. Same STAR index, same cutadapt arguments, same --outFilterMultimapNmax 1, same
# transcriptome output, no bowtie2 stage (mouse path filters ncRNA at pack time).
#
# ADAPTER IS NOT ASSUMED. The cutadapt argument below is a default carried over from Janich; the
# companion check (scripts/check_wang_adapter.sh) runs on the FULL downloaded files and must be read
# before this job is trusted. Adapter-checking a partial file reports nonsense
# (feedback_adapter_check_full_file): a streamed byte-range slice once reported 0% where the complete
# file was 93%. Sanity-check against insert biology too -- footprints are ~30 nt, so a "trimmed" read
# still 75 nt long is a contradiction, not a pass.
#
#   sbatch --array=0-1 scripts/align_wang_ribo.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
D=$NEW/data/external/wang2021_ribo_refetch
BAMDIR=$NEW/data/heldout_bam/mouse_wang_liver_ribo
ADAPTER=${ADAPTER:-AGATCGGAAGAGC}
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$BAMDIR"

mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1}' "$D/manifest.tsv")
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
FQ=$D/${SRR}.fastq.gz
OUTBAM=$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$FQ" ]] || { echo "missing $FQ" >&2; exit 1; }
[[ -s "$OUTBAM" ]] && { echo "SKIP $SRR (bam exists)"; exit 0; }
SCR=/data/tmp/emalekos/wang_ribo_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "run=$SRR adapter=$ADAPTER host=$(hostname) start=$(date -Iseconds)"

"$ENV/bin/cutadapt" -j "$T" -a "$ADAPTER" --minimum-length 20 --maximum-length 40 \
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

# `--quantMode TranscriptomeSAM` ALWAYS emits a read-ordered BAM with no `@HD SO:` tag, whatever
# --outSAMtype says, and prepare_from_bams calls `samtools index`. Sort here rather than in a
# follow-up job -- omitting this is what broke the Janich rebuild twice.
"$ENV/bin/samtools" sort -@ "$T" -m 3G -o "$SCR/${SRR}.sorted.bam" "$TBAM"
before=$("$ENV/bin/samtools" view -@ "$T" -c "$TBAM")
after=$("$ENV/bin/samtools" view -@ "$T" -c "$SCR/${SRR}.sorted.bam")
echo "  records before=$before after=$after"
[[ "$before" -eq "$after" ]] || { echo "ABORT: record count changed ($before -> $after)" >&2; exit 1; }
"$ENV/bin/samtools" view -H "$SCR/${SRR}.sorted.bam" | sed -n 1p | grep -q SO:coordinate \
  || { echo "ABORT: sorted BAM lacks SO:coordinate" >&2; exit 1; }

cp "$SCR/${SRR}.sorted.bam" "$OUTBAM"
"$ENV/bin/samtools" index -@ "$T" "$OUTBAM"
echo "done=$(date -Iseconds) -> $OUTBAM"
