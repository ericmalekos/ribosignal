#!/usr/bin/env bash
#SBATCH --job-name=sc_diag
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/sc_diag_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/sc_diag_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
RB=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/HBL1_pilot/riboseq
SCR=/data/tmp/emalekos/sc_diag_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -4000000 > raw.fq; set -o pipefail   # 1M reads
cutadapt -j 12 -a AGATCGGAAGAGC --minimum-length 25 -o s1.fq raw.fq >/dev/null 2>&1
cutadapt -j 12 -u -5 --trim-n -m 20 -M 40 -o trim.fq s1.fq >/dev/null 2>&1
STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn trim.fq \
  --outSAMtype SAM --outFilterMultimapNmax 20 --alignEndsType Local \
  --outFilterMismatchNoverLmax 0.10 --outFilterMatchNminOverLread 0.66 \
  --outFileNamePrefix "$SCR/d." --outTmpDir "$SCR/_tmp" >/dev/null 2>&1
echo "=== soft-clip diagnostic on + strand primary mapped reads ==="
# + strand primary mapped: flag & 16 == 0, & 256 == 0, & 4 == 0
awk 'BEGIN{FS="\t"}
  !/^@/ && and($2,16)==0 && and($2,256)==0 && and($2,4)==0 {
    n++; c=$6; seq=$10; lead=0; trail=0;
    if (match(c,/^[0-9]+S/)) { lead=substr(c,1,RLENGTH-1)+0 }
    if (match(c,/[0-9]+S$/)) { t=substr(c,RSTART,RLENGTH-1); trail=t+0 }
    LH[lead]++; TH[trail]++;
    if (lead>0) { LS[substr(seq,1,lead)]++; nlead++ }
    if (trail>0){ TS[substr(seq,length(seq)-trail+1,trail)]++; ntrail++ }
  }
  END{
    printf "primary +strand mapped: %d;  with 5prime clip: %d (%.1f%%);  with 3prime clip: %d (%.1f%%)\n", n, nlead, 100*nlead/n, ntrail, 100*ntrail/n
    print "\n-- 5prime clip length histogram --"; for(k in LH) printf "  %s\t%d\n",k,LH[k] | "sort -n | head -12"
    close("sort -n | head -12")
    print "\n-- 3prime clip length histogram --"; for(k in TH) printf "  %s\t%d\n",k,TH[k] | "sort -n | head -12"
    close("sort -n | head -12")
    print "\n-- top 15 5prime-clip sequences --"; for(s in LS) printf "  %d\t%s\n",LS[s],s | "sort -rn | head -15"
    close("sort -rn | head -15")
    print "\n-- top 15 3prime-clip sequences --"; for(s in TS) printf "  %d\t%s\n",TS[s],s | "sort -rn | head -15"
    close("sort -rn | head -15")
  }' "$SCR/d.Aligned.out.sam"
rm -rf "$SCR"