#!/usr/bin/env bash
#SBATCH --job-name=ribo_bw
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/ribo_bw_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/ribo_bw_%A_%a.err
#
#   sbatch --array=0-66 --export=ALL,SHEET=<ribo_ALL_genome.tsv> scripts/xspecies/ribo_genome_bigwig.sbatch
#
# GENOME-COORDINATE Ribo-seq BAMs + bigwigs, for browser viewing.
#
# WHY THIS EXISTS AT ALL. scripts/riboseq_align.sbatch runs STAR with
# `--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM` and copies out ONLY the
# transcriptome BAM; the genomic one is written into node-local $SCR and removed by that
# job's own EXIT trap. So no genome-coordinate Ribo BAM exists for any of the 67 runs and
# bamCoverage cannot be pointed at the ones on disk: their references are TRANSCRIPTS
# (201,583 @SQ lines for human), which no genome browser can show against an assembly.
# This job re-aligns from the same FASTQs under the SAME posture and keeps the genomic BAM.
#
# POSTURE MUST MATCH scripts/riboseq_align.sbatch OR THE TRACKS DO NOT DESCRIBE THE CALLS:
#   --outFilterMultimapNmax 1   unique alignments only (project-wide Ribo-seq rule)
#   --alignEndsType EndToEnd    no soft-clipping, footprint ends are the signal
#   identical cutadapt per run, taken from the same samplesheet column
# The one deliberate difference is the ncRNA drop. The transcriptome pipeline drops ncRNA
# by TRANSCRIPT ID (filter_tx_heldout.py); in genome space the equivalent is dropping
# reads that fall in rRNA/tRNA/miRNA/Mt loci, which is what data/annot/<sp>/ncrna_postfilter.bed
# is for. Same intent, different coordinate system.
#
# NO NORMALIZATION and --binSize 1, as requested: the bigwig carries raw per-base read
# depth, directly comparable to the P-site profiles the model is scored against. Note this
# means tracks are NOT comparable ACROSS runs without accounting for library size.
set -euo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
SHEET=${SHEET:?set SHEET=<ribo samplesheet.tsv>}
BAMROOT=${BAMROOT:-$NEW/data/xspecies_ribo_bam_mm1_genome}
BWROOT=${BWROOT:-$NEW/data/xspecies_ribo_bigwig_genome}
T=${SLURM_CPUS_PER_TASK:-12}
export PYTHONNOUSERSITE=1

mapfile -t ROWS < <(awk -F'\t' 'NR>1' "$SHEET")
IFS=$'\t' read -r DS RUN FQ ADAPTER SPECIES UMIPAT CUTEXTRA <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
read -r -a CUTX <<< "${CUTEXTRA:-}"

# Species resolution mirrors riboseq_align.sbatch: xspecies rows come from species.tsv so
# an eighth species needs a table row and no edit here.
SPTSV=$NEW/scripts/xspecies/species.tsv
ANNVER=$(awk -F'\t' -v s="$SPECIES" '!/^#/ && $1==s {print $2; exit}' "$SPTSV" 2>/dev/null)
[[ -n "$ANNVER" ]] || { echo "unknown species '$SPECIES' (no row in $SPTSV)" >&2; exit 1; }
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_${SPECIES}_${ANNVER}
NCBED=$NEW/data/annot/$SPECIES/ncrna_postfilter.bed
# An empty interval file would drop nothing and report success, which is the exact silent
# failure the normalizer's assertions exist to prevent. Refuse to run on one.
[[ -s "$NCBED" ]] && [[ $(wc -l < "$NCBED") -ge 10 ]] || {
  echo "ncRNA interval BED $NCBED missing or implausibly small" >&2; exit 1; }
for f in "$FQ" "$IDX"; do [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }; done

BAMDIR=$BAMROOT/$DS; BWDIR=$BWROOT/$DS; mkdir -p "$BAMDIR" "$BWDIR"
OUTBAM=$BAMDIR/${RUN}.genome.bam
OUTBW=$BWDIR/${RUN}.bw
[[ -s "$OUTBW" && -s "$OUTBAM.bai" ]] && { echo "SKIP $RUN (done)"; exit 0; }

SCR=/data/tmp/emalekos/ribobw_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "dataset=$DS run=$RUN species=$SPECIES adapter=$ADAPTER host=$(hostname) start=$(date -Iseconds)"

# ---- trim, identically to the transcriptome run -----------------------------
if [[ "$ADAPTER" == "--trim-n" ]]; then
  "$ENV/bin/cutadapt" --trim-n "${CUTX[@]}" -m 20 -M 40 -j "$T" -o "$SCR/t.fq.gz" "$FQ" > "$BAMDIR/${RUN}.cutadapt.log" 2>&1
else
  "$ENV/bin/cutadapt" -a "$ADAPTER" "${CUTX[@]}" -m 20 -M 40 -j "$T" -o "$SCR/t.fq.gz" "$FQ" > "$BAMDIR/${RUN}.cutadapt.log" 2>&1
fi
grep -E "Total reads processed|Reads written" "$BAMDIR/${RUN}.cutadapt.log" | sed 's/^/  /'

# ---- align, genome coordinates this time ------------------------------------
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/t.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --alignEndsType EndToEnd --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM \
  --outBAMsortingThreadN "$T" --outFileNamePrefix "$SCR/${RUN}."
cp "$SCR/${RUN}.Log.final.out" "$BAMDIR/${RUN}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$BAMDIR/${RUN}.STAR.Log.final.out" | sed 's/^/  /'

GB=$SCR/${RUN}.Aligned.sortedByCoord.out.bam
[[ -s "$GB" ]] || { echo "STAR produced no genomic BAM" >&2; exit 1; }
"$ENV/bin/samtools" quickcheck "$GB" || { echo "ABORT: STAR genomic BAM failed quickcheck" >&2; exit 1; }

# ---- drop rRNA / tRNA / miRNA / Mt loci -------------------------------------
# -L selects reads overlapping the BED; -U writes the ones that do NOT, which is what we keep.
before=$("$ENV/bin/samtools" view -@ "$T" -c "$GB")
"$ENV/bin/samtools" index -@ "$T" "$GB"
"$ENV/bin/samtools" view -@ "$T" -b -L "$NCBED" -U "$SCR/keep.bam" -o "$SCR/ncrna.bam" "$GB"
after=$("$ENV/bin/samtools" view -@ "$T" -c "$SCR/keep.bam")
echo "  ncRNA-locus drop: $before -> $after records ($(awk -v a=$before -v b=$after 'BEGIN{printf "%.2f", (a?100*(a-b)/a:0)}')% dropped)"
[[ "$after" -gt 0 ]] || { echo "ABORT: ncRNA filter removed every read" >&2; exit 1; }

cp "$SCR/keep.bam" "$OUTBAM"
"$ENV/bin/samtools" quickcheck "$OUTBAM" || { echo "ABORT: $OUTBAM truncated" >&2; exit 1; }
"$ENV/bin/samtools" view -H "$OUTBAM" | sed -n 1p | grep -q SO:coordinate \
  || { echo "ABORT: output BAM is not coordinate-sorted" >&2; exit 1; }
"$ENV/bin/samtools" index -@ "$T" "$OUTBAM"

# ---- bigwig -----------------------------------------------------------------
"$ENV/bin/bamCoverage" --bam "$OUTBAM" --outFileName "$OUTBW" \
  --binSize 1 --numberOfProcessors "$T" --outFileFormat bigwig
[[ -s "$OUTBW" ]] || { echo "ABORT: no bigwig written" >&2; exit 1; }
echo "  bigwig $(du -h "$OUTBW" | cut -f1)   bam $(du -h "$OUTBAM" | cut -f1)"
echo "=== $DS/$RUN DONE=$(date -Iseconds)"
