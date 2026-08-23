#!/usr/bin/env bash
#SBATCH --job-name=cho_nh
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --array=0-73%10
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/cho_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/cho_%A_%a.err
set -euo pipefail
# EXPERIMENT (a): how much of the TRAINING SET's Ribo-seq signal did --outFilterMultimapNmax 1
# throw away, per transcript?
#
# Re-runs the Chothani training alignment at mm25 -- identical cutadapt, identical TWO-STAGE bowtie2
# contaminant + snoRNA depletion, identical STAR flags including --alignEndsType EndToEnd -- and
# differs ONLY in the multimap cap. Then counts, per transcript, alignments with NH==1 vs NH>1.
#
# NO BAM IS KEPT. At 92.8% quota, 79 mm25 BAMs would be several hundred GB. The BAM lives in
# node-local scratch, is counted, and dies with the job; only a small per-run TSV of nonzero
# transcripts persists. Re-deriving a BAM costs one STAR run, so nothing irreproducible is lost.
#
# DIAGNOSTIC. mm25 output must never feed a pack or a training run -- see
# data/mm25_diagnostic/README.md.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
CONTAM=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/contam_grch38_v49
SNORNA=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/snorna_grch38_v49
D=$NEW/data/external/chothani_ribo_refetch
OUT=$NEW/data/mm25_diagnostic/chothani_nh_counts
T=${SLURM_CPUS_PER_TASK:-12}
mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1}' "$D/manifest.tsv")
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
TSV=$OUT/${SRR}.nh.tsv.gz
[[ -s "$TSV" ]] && { echo "SKIP $SRR (counts exist)"; exit 0; }
FQ=$D/${SRR}.fastq.gz
[[ -s "$FQ" ]] || { echo "missing $FQ" >&2; exit 1; }
SCR=/data/tmp/emalekos/chonh_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "run=$SRR host=$(hostname) start=$(date -Iseconds)"
"$ENV/bin/cutadapt" --trim-n -m 20 -M 40 -j "$T" -o "$SCR/trim.fq.gz" "$FQ" > "$SCR/cutadapt.log" 2>&1
"$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$CONTAM" -U "$SCR/trim.fq.gz" \
  --un-gz "$SCR/c1.fq.gz" -S /dev/null 2> "$SCR/bt_contam.log"
"$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$SNORNA" -U "$SCR/c1.fq.gz" \
  --un-gz "$SCR/clean.fq.gz" -S /dev/null 2> "$SCR/bt_snorna.log"
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/clean.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --alignEndsType EndToEnd --outFilterMultimapNmax 25 \
  --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${SRR}."
TBAM=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$TBAM" ]] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
u=$(grep "Uniquely mapped reads %" "$SCR/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
m=$(grep "% of reads mapped to multiple loci" "$SCR/${SRR}.Log.final.out" | grep -oE "[0-9.]+%")
# per-transcript: alignments from uniquely-mapped reads vs from multimappers
"$ENV/bin/samtools" view -@ "$T" "$TBAM" \
  | awk -F'\t' '{nh=0; for(i=12;i<=NF;i++) if($i ~ /^NH:i:/){split($i,a,":"); nh=a[3]; break}
                 if(nh==1) u[$3]++; else if(nh>1) m[$3]++}
       END{for(t in u) print t"\t"u[t]"\t"(t in m ? m[t] : 0);
           for(t in m) if(!(t in u)) print t"\t0\t"m[t]}' \
  | gzip -c > "$TSV"
echo "DONE $SRR unique=$u multi=$m  tx_with_signal=$(zcat "$TSV" | wc -l)  bam DISCARDED"
