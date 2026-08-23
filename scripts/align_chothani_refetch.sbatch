#!/usr/bin/env bash
#SBATCH --job-name=cho_align
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/cho_align_%A_%a.out
# SUPERSEDED by scripts/riboseq_align.sbatch (docs/PIPELINE_POLICY.md rule 1).
# Kept for provenance only -- do NOT use for new work. Its Chothani re-alignment is now a
# samplesheet: data/chothani_ribo_samplesheet.tsv.
set -euo pipefail
# Re-align the 74 Chothani Ribo-seq runs the RELEASED MODELS WERE TRAINED ON.
#
# WHY. audit_pack_inputs.py found every `*_psites.hd5` behind the training panel had been deleted:
# 14 packs referencing 120 missing source files, 100 unique accessions gone. The packs still work and
# the models are already trained, so nothing is broken today -- what is missing is the ABILITY to
# rebuild the training data. The FASTQs were re-fetched from ENA (74/74, md5-verified); this restores
# the alignment layer.
#
# READ FORMAT, verified empirically rather than assumed. Reads are fixed 35 nt (98.7%) with NO ligated
# adapter (TruSeq 0.1%, Ingolia 0%, smallRNA 0%). The 3' terminal k-mers are dominated by N runs: the
# reads are N-PADDED to a fixed length, the same layout as GSE208041's RPF arm. `cutadapt --trim-n`
# recovers mode 28 nt with 98.8% of reads in 20-40 nt, a clean footprint distribution. Trimming for a
# non-existent adapter, or skipping the trim, would both be wrong.
#
# EXPECT DIVERGENCE FROM THE EXISTING PACKS, AND MEASURE IT. The original recipe is unrecoverable --
# data/packed has no provenance.json at all. This uses the project's current human Ribo standard.
# Independently, RiboCode metaplots re-derives periodic read lengths per sample, which alone cost
# 1.69% on the Janich rebuild even with identical inputs. So the deliverable is a documented,
# reproducible path to a training pack plus a measurement of how far it lands from the pack the
# released models actually saw -- NOT a claim of byte-identity.
#
#   sbatch --array=0-73%8 scripts/align_chothani_refetch.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
CONTAM=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/contam_grch38_v49
SNORNA=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/snorna_grch38_v49
D=$NEW/data/external/chothani_ribo_refetch
BAMDIR=$NEW/data/chothani_bam
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$BAMDIR"

mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1}' "$D/manifest.tsv")
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
FQ=$D/${SRR}.fastq.gz
OUTBAM=$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$FQ" ]] || { echo "missing $FQ" >&2; exit 1; }
[[ -s "$OUTBAM" ]] && { echo "SKIP $SRR (bam exists)"; exit 0; }
SCR=/data/tmp/emalekos/cho_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "run=$SRR host=$(hostname) start=$(date -Iseconds)"

"$ENV/bin/cutadapt" --trim-n -m 20 -M 40 -j "$T" -o "$SCR/trim.fq.gz" "$FQ" \
  > "$BAMDIR/${SRR}.cutadapt.log" 2>&1
grep -E "Total reads processed|Reads written" "$BAMDIR/${SRR}.cutadapt.log" | sed 's/^/  /'

"$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$CONTAM" -U "$SCR/trim.fq.gz" \
  --un-gz "$SCR/c1.fq.gz" -S /dev/null 2> "$BAMDIR/${SRR}.bowtie2_contam.log"
grep -E "overall alignment rate" "$BAMDIR/${SRR}.bowtie2_contam.log" | sed 's/^/  contam: /'
"$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$SNORNA" -U "$SCR/c1.fq.gz" \
  --un-gz "$SCR/clean.fq.gz" -S /dev/null 2> "$BAMDIR/${SRR}.bowtie2_snorna.log"
grep -E "overall alignment rate" "$BAMDIR/${SRR}.bowtie2_snorna.log" | sed 's/^/  snoRNA: /'

"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/clean.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --alignEndsType EndToEnd --outFilterMultimapNmax 1 \
  --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${SRR}."
cp "$SCR/${SRR}.Log.final.out" "$BAMDIR/${SRR}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$BAMDIR/${SRR}.STAR.Log.final.out" | sed 's/^/  /'

TBAM=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
[[ -s "$TBAM" ]] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
# TranscriptomeSAM is read-ordered; prepare_from_bams calls samtools index. Sort HERE, in the same
# job, rather than in a follow-up -- omitting this cost two runs on the Janich rebuild.
"$ENV/bin/samtools" sort -@ "$T" -m 3G -o "$SCR/sorted.bam" "$TBAM"
before=$("$ENV/bin/samtools" view -@ "$T" -c "$TBAM")
after=$("$ENV/bin/samtools" view -@ "$T" -c "$SCR/sorted.bam")
[[ "$before" -eq "$after" ]] || { echo "ABORT: record count changed ($before -> $after)" >&2; exit 1; }
cp "$SCR/sorted.bam" "$OUTBAM"
"$ENV/bin/samtools" index -@ "$T" "$OUTBAM"
echo "done=$(date -Iseconds) -> $OUTBAM"
