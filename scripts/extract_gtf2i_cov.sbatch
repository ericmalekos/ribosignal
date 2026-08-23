#!/usr/bin/env bash
#SBATCH --job-name=g2i_cov
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --array=0-9%10
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/g2icov_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/g2icov_%A_%a.err
set -euo pipefail
# Per-position RNA coverage over GTF2I ENST00000901263.1 for each multimap posture.
# Streams with samtools+awk rather than pysam: the BAMs are 11 GB each and a per-record Python
# loop over 110 GB is hours; grep-speed filtering is minutes, parallel across 10 jobs.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SAM=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
D=$NEW/data/mm25_diagnostic/hep_rna_mm
TX=ENST00000901263.1
declare -a SPEC
for tag in mm1 mm10 mm25 best_Old_2_4 best_Random; do
  for s in SRR15513268 SRR15513269; do SPEC+=("$s:$tag"); done; done
IFS=: read -r SRR TAG <<< "${SPEC[$SLURM_ARRAY_TASK_ID]}"
B=$D/${SRR}.${TAG}.toTranscriptome.bam
OUT=$D/gtf2i/${SRR}.${TAG}.cov.tsv
[ -s "$OUT" ] && { echo "SKIP"; exit 0; }
[ -s "$B" ] || { echo "missing $B" >&2; exit 1; }
"$SAM" view -@ 4 "$B" \
  | awk -v tx="$TX" -F'\t' '$3==tx{
      pos=$4; cig=$6; n="";
      for(i=1;i<=length(cig);i++){ c=substr(cig,i,1);
        if(c ~ /[0-9]/){ n=n c } else {
          if(c=="M"||c=="D"||c=="="||c=="X"){ for(k=0;k<n+0;k++) cov[pos+k]++ ; pos+=n }
          else if(c=="N"){ pos+=n }
          n="" } } }
     END{ for(p=1;p<=5178;p++) print p"\t"(cov[p]+0) }' > "$OUT"
echo "DONE $SRR $TAG  total=$(awk -F'\t' '{s+=$2}END{print s}' "$OUT")"
