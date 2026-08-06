#!/usr/bin/env bash
#SBATCH --job-name=e2e_trim
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/e2e_trim_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/e2e_trim_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
RB=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/HBL1_pilot/riboseq
SCR=/data/tmp/emalekos/e2e_trim_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -8000000 > raw.fq; set -o pipefail   # 2M reads

# candidate trims -> fastq
cutadapt -j12 -a AGATCGGAAGAGC --minimum-length 25 -o A1.fq raw.fq >/dev/null 2>&1
cutadapt -j12 -u -5 --trim-n -m 20 -M 40 -o A.fq A1.fq >/dev/null 2>&1                       # A: current
cutadapt -j12 -a ATCGTAGATCGGAAGAGC --trim-n -m 20 -M 40 -o B.fq raw.fq >/dev/null 2>&1        # B: full tail
cutadapt -j12 -u 2 -a ATCGTAGATCGGAAGAGC --trim-n -m 20 -M 40 -o C.fq raw.fq >/dev/null 2>&1   # C: +5' trim 2
cutadapt -j12 -a ATCGTAGATCGGAAGAGC -e 0.2 --trim-n -m 20 -M 40 -o D.fq raw.fq >/dev/null 2>&1 # D: looser adapter err

e2e() {  # name fastq
  local name=$1 fq=$2
  STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn "$fq" \
    --outSAMtype None --outFilterMultimapNmax 20 --alignEndsType EndToEnd \
    --outFilterMismatchNoverLmax 0.07 --outFilterMatchNminOverLread 0.8 \
    --outFileNamePrefix "$SCR/${name}." --outTmpDir "$SCR/_t_${name}" >/dev/null 2>&1
  local u=$(grep "Uniquely mapped reads %" "$SCR/${name}.Log.final.out"|grep -oE "[0-9.]+%")
  local un=$(grep "Uniquely mapped reads number" "$SCR/${name}.Log.final.out"|grep -oE "[0-9]+")
  local mm=$(grep "too many mismatches" "$SCR/${name}.Log.final.out"|grep -oE "[0-9.]+%"|head -1)
  local sh=$(grep "unmapped: too short" "$SCR/${name}.Log.final.out"|grep -oE "[0-9.]+%"|head -1)
  local ml=$(grep "mapped to multiple loci" "$SCR/${name}.Log.final.out"|grep -oE "[0-9.]+%"|head -1)
  local nr=$(( $(wc -l < "$fq") / 4 ))
  printf "%-3s reads=%-9d  UNIQUE=%-7s (%s)  multi=%-6s  too_mm=%-6s  too_short=%s\n" "$name" "$nr" "$u" "$un" "$ml" "$mm" "$sh"
}
echo "=== EndToEnd (mm0.07, matchNmin0.8) unique% by trimming strategy ==="
echo "A = -a AGATCGGAAGAGC then -u -5 (current)"
echo "B = -a ATCGTAGATCGGAAGAGC (full tail as one adapter)"
echo "C = -u 2 (5' trim) then -a ATCGTAGATCGGAAGAGC"
echo "D = -a ATCGTAGATCGGAAGAGC -e 0.2 (looser adapter error)"
echo
e2e A A.fq
e2e B B.fq
e2e C C.fq
e2e D D.fq
rm -rf "$SCR"