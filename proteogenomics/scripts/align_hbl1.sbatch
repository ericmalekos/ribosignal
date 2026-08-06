#!/usr/bin/env bash
#SBATCH --job-name=align_hbl1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_hbl1_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_hbl1_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
HB=$NEW/proteogenomics/data/HBL1_pilot
R1=$HB/rnaseq/SRR12285172_1.fastq.gz; R2=$HB/rnaseq/SRR12285172_2.fastq.gz
T=${SLURM_CPUS_PER_TASK:-8}; SCR=/data/tmp/emalekos/hbl1_${SLURM_JOB_ID:-m}; mkdir -p "$SCR" "$HB/salmon" "$HB/coverage"
echo "host=$(hostname) start=$(date -Iseconds)"
"$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "$R1" -2 "$R2" --gcBias --seqBias -p "$T" -o "$HB/salmon/HBL1"
echo "salmon done: tx>1TPM $(awk 'NR>1 && $4>1' "$HB/salmon/HBL1/quant.sf" | wc -l)"
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/HBL1."
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$SCR/HBL1.Aligned.toTranscriptome.out.bam" "$HB/coverage/HBL1_coverage.hd5" HBL1
rm -rf "$SCR"; echo "done=$(date -Iseconds)"
