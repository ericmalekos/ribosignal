#!/usr/bin/env bash
#SBATCH --job-name=hbl1_pilot
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_pilot_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_pilot_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
ANNOT=$ECH/data/ribocode_annot
RB=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/HBL1_pilot/riboseq
SCR=/data/tmp/emalekos/hbl1_pilot_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -80000000 > raw.fq; set -o pipefail   # 20M reads
echo "=== HBL-1 pilot: paper-matched STAR params (EndToEnd, seedSearchStartLmax 15, mm0.05, matchNmin 25) ==="
echo "=== 5' UMI trim sweep u=0..3; 3' = -a NNNNNAGATCGGAAGAGC (5nt UMI+adapter); footprint 26-34 ==="
echo "canonical P-site offset ~12; pick the u giving offset~12 + highest frame-0%"
echo
for U in 0 1 2 3; do
  cutadapt -j12 -u $U -a NNNNNAGATCGGAAGAGC --trim-n --minimum-length 26 --maximum-length 34 \
      -o t$U.fq raw.fq >/dev/null 2>&1
  STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn t$U.fq \
    --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
    --alignEndsType EndToEnd --seedSearchStartLmax 15 \
    --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 25 --outFilterMultimapNmax 20 \
    --outSAMattributes NH HI AS nM MD \
    --outFileNamePrefix "$SCR/u$U." --outTmpDir "$SCR/_t$U" >/dev/null 2>&1
  uniq=$(grep "Uniquely mapped reads %" u$U.Log.final.out|grep -oE "[0-9.]+%")
  uniqn=$(grep "Uniquely mapped reads number" u$U.Log.final.out|grep -oE "[0-9]+")
  # sort+index transcriptome BAM for metaplots
  samtools sort -@12 -o u$U.tx.bam u$U.Aligned.toTranscriptome.out.bam 2>/dev/null
  samtools index u$U.tx.bam
  echo "$SCR/u$U.tx.bam" > u$U.list
  "$RC/metaplots" -a "$ANNOT" -i u$U.list -o u$U >/dev/null 2>&1 || true
  echo "----- 5'trim=$U : unique=$uniq ($uniqn) -----"
  # metaplots periodic read lengths: read_len, offset(field3), f0%(field7)
  grep -E "^# [0-9]+\s" u$U"_pre_config.txt" 2>/dev/null | \
    awk -F'\t' '{printf "   len=%s  offset=%s  frame0=%s  (f0/f1/f2=%s/%s/%s)\n",$1,$3,$7,$4,$5,$6}' | sed 's/# //' | head
  [ -s u$U"_pre_config.txt" ] || echo "   (metaplots: no periodic length passed threshold)"
done
rm -rf "$SCR"