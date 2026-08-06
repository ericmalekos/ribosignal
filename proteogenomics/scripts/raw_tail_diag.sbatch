#!/usr/bin/env bash
#SBATCH --job-name=tail_diag
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/tail_diag_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/tail_diag_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
RB=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/HBL1_pilot/riboseq
SCR=/data/tmp/emalekos/tail_diag_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"
# RAW reads, NO trimming; just drop trailing N. Align Local so the footprint maps + the whole tail soft-clips.
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -2000000 > raw.fq; set -o pipefail   # 500k reads
STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn raw.fq \
  --outSAMtype SAM --outFilterMultimapNmax 20 --alignEndsType Local \
  --outFilterMismatchNoverLmax 0.10 --outFilterMatchNminOverLread 0.30 --outFilterScoreMinOverLread 0.30 \
  --outFileNamePrefix "$SCR/d." --outTmpDir "$SCR/_tmp" >/dev/null 2>&1
echo "=== RAW-read 3' tail consensus (+ strand primary mapped); reveals [linker][adapter] ==="
# For + strand primary mapped reads with a 5'-anchored match, dump the trailing soft-clip (= technical tail).
awk 'BEGIN{FS="\t"}
  !/^@/ && and($2,16)==0 && and($2,256)==0 && and($2,4)==0 {
    c=$6; seq=$10;
    if (match(c,/[0-9]+S$/)) {
      t=substr(c,RSTART,RLENGTH-1)+0;
      if (t>=10) {
        tail=substr(seq,length(seq)-t+1,t); n++;
        for(i=1;i<=18 && i<=length(tail);i++){ b=substr(tail,i,1); cnt[i","b]++ }
      }
    }
  }
  END{
    printf "reads with >=10nt 3prime tail: %d\n", n;
    print "pos  A%    C%    G%    T%   consensus";
    cons="";
    for(i=1;i<=18;i++){
      a=cnt[i",A"]; c=cnt[i",C"]; g=cnt[i",G"]; t=cnt[i",T"]; tot=a+c+g+t; if(tot==0)tot=1;
      split("A C G T",B," "); mx=0; mb="N";
      if(a>mx){mx=a;mb="A"} if(c>mx){mx=c;mb="C"} if(g>mx){mx=g;mb="G"} if(t>mx){mx=t;mb="T"}
      cons=cons mb;
      printf "%2d  %5.1f %5.1f %5.1f %5.1f    %s\n", i, 100*a/tot,100*c/tot,100*g/tot,100*t/tot, mb;
    }
    print "\nCONSENSUS TAIL (pos1-18): " cons;
  }' "$SCR/d.Aligned.out.sam"
rm -rf "$SCR"