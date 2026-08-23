#!/usr/bin/env bash
#SBATCH --job-name=hep_best
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --array=0-3%4
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/hepbest_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/hepbest_%A_%a.err
set -euo pipefail
# BEST-ALIGNMENT-ONLY posture: --outFilterMultimapNmax 25 --outSAMmultNmax 1
#
# Distinct from both arms already running:
#   --outFilterMultimapNmax 1   DISCARDS a read that maps to >1 locus  (~30% of reads lost)
#   --outSAMmultNmax 1          KEEPS the read, emits only its TOP-SCORING alignment
# STAR docs: "if this value is not equal to -1, the top scoring alignment will be output first".
#
# This is a genuine read-level posture C: every read committed to one place, nothing discarded, and
# no 25x count inflation -- plausibly the right setting for the model's RNA COVERAGE INPUT channel.
#
# TIE CAVEAT, which is why both orders are run. True paralogs/repeats produce EQUAL best scores;
# --outSAMprimaryFlag OneBestScore then picks one, ordered by --outMultimapperOrder. The default
# Old_2.4 is a quasi-random internal order; Random is at least unbiased and seedable. At exactly the
# loci this experiment is about, the choice among ties is arbitrary -- comparing the two orders
# measures how much that arbitrariness actually moves the coverage.
#
# BAMs KEPT (CLAUDE.md rule 3).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
FQD=$NEW/data/external/chothani_rna/fastq
OUT=$NEW/data/mm25_diagnostic/hep_rna_mm
declare -a SPEC
for ord in Old_2.4 Random; do for s in SRR15513268 SRR15513269; do SPEC+=("$s:$ord"); done; done
IFS=: read -r SRR ORD <<< "${SPEC[$SLURM_ARRAY_TASK_ID]}"
TAG=best_${ORD//./_}
BAM=$OUT/${SRR}.${TAG}.toTranscriptome.bam
[ -s "$BAM" ] && { echo "SKIP $SRR $TAG"; exit 0; }
R1=$FQD/${SRR}_1.fastq.gz; R2=$FQD/${SRR}_2.fastq.gz
for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done
SCR=/data/tmp/emalekos/hepbest_${SRR}_${TAG}.$$; mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --runThreadN "${SLURM_CPUS_PER_TASK:-12}" --outSAMtype None \
  --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 25 --outSAMmultNmax 1 \
  --outMultimapperOrder "$ORD" --outSAMprimaryFlag OneBestScore \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${SRR}."
T=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$T" ] || { echo "no transcriptome BAM for $SRR $TAG" >&2; exit 1; }
mv "$T" "$BAM"
cp "$SCR/${SRR}.Log.final.out" "$OUT/${SRR}.${TAG}.Log.final.out" 2>/dev/null || true
u=$(grep "Uniquely mapped reads %" "$OUT/${SRR}.${TAG}.Log.final.out" | grep -oE "[0-9.]+%")
m=$(grep "% of reads mapped to multiple loci" "$OUT/${SRR}.${TAG}.Log.final.out" | grep -oE "[0-9.]+%")
echo "DONE $SRR $TAG unique=$u multi=$m  records=$("$ENV/bin/samtools" view -c "$BAM")  BAM KEPT"
