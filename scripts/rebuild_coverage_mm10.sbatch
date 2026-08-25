#!/usr/bin/env bash
#SBATCH --job-name=cov_mm10
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --array=0-56
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm10/cov_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm10/cov_%A_%a.err
set -euo pipefail
# Rebuild the model's RNA COVERAGE INPUT at --outFilterMultimapNmax 10.
#
# WHY. GTF2I ENST00000901263.1 has a 1,617 nt CDS window where the mm1 coverage channel reads 9% of
# the CDS mean; at mm10 it reads 138%. The model predicted ~0 translation there because its INPUT
# said there was almost no RNA. That is a corrupted feature, not memorisation (results.md
# 2026-08-23). mm10 recovers everything mm25 does, so mm10 is the cheaper setting.
#
# THE TARGET IS NOT REBUILT. The Ribo-seq label stays at mm1 by standing rule: RiboCode does not
# drop multimappers and they destroy periodicity. Only the INPUT changes, via
# prepare_pack.py --reuse-target.
#
# BAMs ARE KEPT and archived afterwards (CLAUDE.md: archive every FASTQ and BAM; do not delete
# intermediates until the result is verified).
#
# DISK: ~57 BAMs. Check free space before launching -- see the guard below.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
FQD=$NEW/data/external/chothani_rna/fastq
BAMD=$NEW/data/rnaseq_bam_mm10;  COVD=$NEW/data/rnaseq_coverage/per_sample_mm10
mkdir -p "$BAMD" "$COVD" "$NEW/logs/mm10"
# every RNA run: the 25-run manifest plus Fibroblast's 32
mapfile -t RUNS < <( { awk -F'\t' 'NR>1{print $1}' "$NEW/data/external/chothani_rna/manifest_25.tsv"; \
                       awk -F'\t' 'NR>1{print $1}' "$NEW/data/external/chothani_rna/manifest_fib.tsv"; } | sort -u )
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
OUTCOV=$COVD/${SRR}_coverage.hd5
[ -s "$OUTCOV" ] && { echo "SKIP $SRR (coverage exists)"; exit 0; }
R1=$FQD/${SRR}_1.fastq.gz; R2=$FQD/${SRR}_2.fastq.gz
for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "MISSING $f" >&2; exit 1; }; done
SCR=/data/tmp/emalekos/covmm10_${SRR}.$$; mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "run=$SRR host=$(hostname) start=$(date -Iseconds)"
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --runThreadN "${SLURM_CPUS_PER_TASK:-12}" --outSAMtype None \
  --quantMode TranscriptomeSAM --outFilterMultimapNmax 10 \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${SRR}."
TBAM=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$TBAM" ] || { echo "no transcriptome BAM for $SRR" >&2; exit 1; }
"$ENV/bin/samtools" quickcheck "$TBAM" || { echo "$SRR BAM truncated" >&2; exit 1; }
# coverage first (the deliverable), THEN keep the BAM
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$OUTCOV" "$SRR"
[ -s "$OUTCOV" ] || { echo "no coverage written for $SRR" >&2; exit 1; }
mv "$TBAM" "$BAMD/${SRR}.Aligned.toTranscriptome.out.bam"
cp "$SCR/${SRR}.Log.final.out" "$BAMD/${SRR}.Log.final.out" 2>/dev/null || true
u=$(grep "Uniquely mapped reads %" "$BAMD/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
m=$(grep "% of reads mapped to multiple loci" "$BAMD/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
echo "DONE $SRR unique=$u multi=$m -> $(basename "$OUTCOV")"
