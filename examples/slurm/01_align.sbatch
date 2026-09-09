#!/usr/bin/env bash
#SBATCH --job-name=ribosignal_align
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
# RNA-seq and Ribo-seq alignment. The two assays use opposite settings on purpose; see the
# tutorial page "Align and build coverage".
: "${RIBO_IMAGE:?set RIBO_IMAGE to the .sif}"
: "${RIBO_WORK:?set RIBO_WORK to the working directory}"
: "${RIBO_REF:=$RIBO_WORK/ref}"
T=${SLURM_CPUS_PER_TASK:-16}
S() { singularity exec --cleanenv -B "$RIBO_WORK" -B "$RIBO_REF" "$RIBO_IMAGE" "$@"; }
cd "$RIBO_WORK"

# RNA: keep multimappers, do NOT length-filter, do NOT discard untrimmed
S cutadapt --trim-n -m 20 -j "$T" -o rna/t1.fq.gz -p rna/t2.fq.gz \
  fastq/${RNA_SRR}_1.fastq.gz fastq/${RNA_SRR}_2.fastq.gz
S STAR --genomeDir "$RIBO_REF/star_index" --readFilesIn rna/t1.fq.gz rna/t2.fq.gz \
  --readFilesCommand zcat --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM --outFileNamePrefix "rna/$RNA_SRR."

# Ribo: unique only, length-selected, EndToEnd so the P-site is not shifted by soft-clipping.
# RIBO_ADAPTER is set from the measurement, and EMPTY is a legitimate answer: an already-trimmed
# library given -a --discard-untrimmed loses everything that lacks an adapter it never had.
# Export it as "" deliberately; leaving it unset is what the guard below catches.
[ "${RIBO_ADAPTER+set}" = set ] || { echo "ABORT: export RIBO_ADAPTER (\"\" if none); see the adapters page" >&2; exit 1; }
if [ -n "$RIBO_ADAPTER" ]; then
  S cutadapt -a "$RIBO_ADAPTER" -m 20 -M 40 --discard-untrimmed -j "$T" \
    -o ribo/t.fq.gz "fastq/${RIBO_SRR}.fastq.gz"
else
  S cutadapt --trim-n -m 20 -M 40 -j "$T" -o ribo/t.fq.gz "fastq/${RIBO_SRR}.fastq.gz"
fi
S STAR --genomeDir "$RIBO_REF/star_index" --readFilesIn ribo/t.fq.gz --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --alignEndsType EndToEnd \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "ribo/$RIBO_SRR."

# STAR can exit 0 having read nothing; assert otherwise
for f in "rna/$RNA_SRR.Log.final.out" "ribo/$RIBO_SRR.Log.final.out"; do
  n=$(awk -F'\t' '/Number of input reads/{gsub(/ /,"",$2);print $2}' "$f")
  [ "${n:-0}" -gt 0 ] || { echo "ABORT: STAR read 0 reads ($f)" >&2; exit 1; }
done
