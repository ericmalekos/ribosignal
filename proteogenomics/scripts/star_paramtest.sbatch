#!/usr/bin/env bash
#SBATCH --job-name=star_ptest
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/star_ptest_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/star_ptest_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
RB=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/HBL1_pilot/riboseq
SCR=/data/tmp/emalekos/star_ptest_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"
# trim 2M reads two-step (|| true: head closes the pipe early -> zcat SIGPIPE under pipefail)
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -8000000 > raw.fq; set -o pipefail
cutadapt -j 12 -a AGATCGGAAGAGC --minimum-length 25 -o s1.fq raw.fq >/dev/null 2>&1
cutadapt -j 12 -u -5 --trim-n -m 20 -M 40 -o trim.fq s1.fq >/dev/null 2>&1
echo "trimmed reads: $(( $(wc -l < trim.fq) / 4 ))"
run() {  # name, extra STAR args
  local name=$1; shift
  STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn trim.fq \
    --outSAMtype None --outFilterMultimapNmax 20 --outSAMattributes NH HI \
    --outFileNamePrefix "$SCR/${name}." --outTmpDir "$SCR/_tmp_${name}" "$@" >/dev/null 2>&1
  echo "== $name =="
  grep -E "Uniquely mapped reads %|Uniquely mapped reads number|% of reads mapped to multiple loci|too many mismatches|too short" "$SCR/${name}.Log.final.out"
}
run EndToEnd_mm05 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7
run Local_mm10    --alignEndsType Local    --outFilterMismatchNoverLmax 0.10 --outFilterMatchNminOverLread 0.66
run Local_mm10_sensitive --alignEndsType Local --outFilterMismatchNoverLmax 0.10 --outFilterMatchNminOverLread 0.5 --seedSearchStartLmax 15 --outFilterScoreMinOverLread 0.5
rm -rf "$SCR"