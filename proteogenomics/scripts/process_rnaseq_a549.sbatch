#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_a549
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/rnaseq_a549_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/rnaseq_a549_%j.err
set -euo pipefail
# A549 RNA-seq -> BOTH model inputs, from one FASTQ revert of the open hg19 BAM:
#   (1) salmon decoy-aware v49 isoform TPM  -> expressed-transcript filter + 3-frame null
#   (2) STAR --quantMode TranscriptomeSAM -> rnaseq_coverage.py -> per-nt coverage  -> THE MODEL INPUT
# Faithfully matches the training coverage recipe (align_rnaseq_sample.sh): same STAR params
# (--outFilterMultimapNmax 1, posture A, transcriptome axis), same rnaseq_coverage.py (auto-strand,
# no NH filter). The per-nt coverage is NOT optional -- it is what makes the model's translation
# prediction A549-specific rather than tissue-agnostic. Index paths corrected to STAR_indexes/ +
# Salmon_indexes/ (CLAUDE.md's genomes/ paths are stale/missing).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PILOT=$NEW/proteogenomics/data/A549_pilot
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
BAM=$PILOT/G20499.A549.2.bam
T=${SLURM_CPUS_PER_TASK:-8}
SCR=/data/tmp/emalekos/a549_rnaseq_${SLURM_JOB_ID:-manual}
mkdir -p "$SCR" "$PILOT/salmon" "$PILOT/coverage"
echo "host=$(hostname) start=$(date -Iseconds); BAM=$BAM"

# --- 1) revert coordinate-sorted hg19 BAM -> paired FASTQ (once, shared by salmon + STAR) ---
"$ENV/bin/samtools" collate -@ "$T" -O "$BAM" "$SCR/collate" \
  | "$ENV/bin/samtools" fastq -@ "$T" -1 "$SCR/R1.fq.gz" -2 "$SCR/R2.fq.gz" -0 /dev/null -s /dev/null -n
echo "fastq done $(date -Iseconds): $(du -h $SCR/R1.fq.gz $SCR/R2.fq.gz | cut -f1 | tr '\n' ' ')"

# --- 2) salmon decoy-aware v49 -> isoform TPM ---
"$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "$SCR/R1.fq.gz" -2 "$SCR/R2.fq.gz" \
  --validateMappings --gcBias --seqBias -p "$T" -o "$PILOT/salmon/A549"
echo "salmon done $(date -Iseconds); transcripts >1 TPM: $(awk 'NR>1 && $4>1' "$PILOT/salmon/A549/quant.sf" | wc -l)"

# --- 3) STAR transcriptome align -> per-nt coverage (THE model input), training recipe ---
"$ENV/bin/STAR" --genomeDir "$STARIDX" \
  --readFilesIn "$SCR/R1.fq.gz" "$SCR/R2.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/A549."
TBAM="$SCR/A549.Aligned.toTranscriptome.out.bam"
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCR/A549.Log.final.out" "$PILOT/coverage/A549.Log.final.out" 2>/dev/null || true

"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$PILOT/coverage/A549_coverage.hd5" A549
echo "coverage done $(date -Iseconds) -> $PILOT/coverage/A549_coverage.hd5"

rm -rf "$SCR"
echo "all done $(date -Iseconds)"
