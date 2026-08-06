#!/usr/bin/env bash
#SBATCH --job-name=fp_sweep
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/fp_sweep_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/fp_sweep_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
RB=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/HBL1_pilot/riboseq
SCR=/data/tmp/emalekos/fp_sweep_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -8000000 > raw.fq; set -o pipefail   # 2M reads
echo "=== 5' trim sweep; 3' = -a ATCGTAGATCGGAAGAGC; EndToEnd (mm0.07, matchNmin0.8) ==="
for U in 0 1 2 3 4 5; do
  cutadapt -j12 -u $U -a ATCGTAGATCGGAAGAGC --trim-n -m 20 -M 40 -o t$U.fq raw.fq >/dev/null 2>&1
  STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn t$U.fq \
    --outSAMtype None --outFilterMultimapNmax 20 --alignEndsType EndToEnd \
    --outFilterMismatchNoverLmax 0.07 --outFilterMatchNminOverLread 0.8 \
    --outFileNamePrefix "$SCR/u$U." --outTmpDir "$SCR/_t$U" >/dev/null 2>&1
  u=$(grep "Uniquely mapped reads %" $SCR/u$U.Log.final.out|grep -oE "[0-9.]+%")
  un=$(grep "Uniquely mapped reads number" $SCR/u$U.Log.final.out|grep -oE "[0-9]+")
  mm=$(grep "too many mismatches" $SCR/u$U.Log.final.out|grep -oE "[0-9.]+%"|head -1)
  sh=$(grep "unmapped: too short" $SCR/u$U.Log.final.out|grep -oE "[0-9.]+%"|head -1)
  printf "5'trim=%d  UNIQUE=%-7s (%s)  too_mm=%-6s  too_short=%s\n" "$U" "$u" "$un" "$mm" "$sh"
done
rm -rf "$SCR"