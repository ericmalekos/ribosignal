#!/usr/bin/env bash
#SBATCH --job-name=wang_rna
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/wang_rna_%A_%a.out
set -euo pipefail
# Align the two Wang 2021 mouse-liver RNA-seq runs (SRR5262874, SRR5262875).
#
# WHY. The 3x3 factorial varies the RNA input across three datasets, so it needs Wang RNA coverage,
# not just Wang Ribo. `data/heldout_bam/` has NO Wang BAMs at all -- neither arm was ever aligned there
# (or both were deleted). The RNA FASTQs survive in data/external/wang2021_mouse/fastq/mouse_liver/
# and md5-verify against manifest_liver.tsv, so this needs no download and runs immediately.
#
# RNA ARM, NOT RPF. Differences from align_wang_ribo.sbatch, both deliberate:
#   * no --maximum-length: RNA inserts are not ~30 nt footprints and capping them would throw away
#     most of the library.
#   * no --discard-untrimmed: that flag is correct for RPF, where an untrimmed read means no insert
#     was seen, and wrong for RNA, where a read simply may not reach the adapter. Applying it to the
#     CAR-T RNA arm earlier kept only 16.4% of reads; without it, 99.0%.
# --outFilterMultimapNmax 1 is kept, matching every other arm in the factorial.
#
#   sbatch --array=0-1 scripts/align_wang_rna.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
D=$NEW/data/external/wang2021_mouse/fastq/mouse_liver
BAMDIR=$NEW/data/heldout_bam/mouse_wang_liver_rna
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$BAMDIR"

RUNS=(SRR5262874 SRR5262875)
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
FQ=$D/${SRR}.fastq.gz
OUTBAM=$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$FQ" ]] || { echo "missing $FQ" >&2; exit 1; }
[[ -s "$OUTBAM" ]] && { echo "SKIP $SRR (bam exists)"; exit 0; }
SCR=/data/tmp/emalekos/wang_rna_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "run=$SRR host=$(hostname) start=$(date -Iseconds)"

"$ENV/bin/cutadapt" -j "$T" -a AGATCGGAAGAGC --minimum-length 20 \
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
"$ENV/bin/samtools" sort -@ "$T" -m 3G -o "$SCR/${SRR}.sorted.bam" "$TBAM"
before=$("$ENV/bin/samtools" view -@ "$T" -c "$TBAM")
after=$("$ENV/bin/samtools" view -@ "$T" -c "$SCR/${SRR}.sorted.bam")
echo "  records before=$before after=$after"
[[ "$before" -eq "$after" ]] || { echo "ABORT: record count changed ($before -> $after)" >&2; exit 1; }
cp "$SCR/${SRR}.sorted.bam" "$OUTBAM"
"$ENV/bin/samtools" index -@ "$T" "$OUTBAM"
echo "done=$(date -Iseconds) -> $OUTBAM"
