#!/usr/bin/env bash
#SBATCH --job-name=cart_proc
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/cart_proc_%A_%a.out
# SUPERSEDED by scripts/riboseq_align.sbatch (docs/PIPELINE_POLICY.md rule 1).
# PARTIALLY superseded: the Ribo ALIGNMENT half is replaced (see data/human_ribo_canon_samplesheet.tsv).
# The RNA arm and salmon steps have no canonical equivalent yet. TWO defects, both fixed in the
# canonical path: no filter_tx_heldout, and umi_tools dedup ran on the GENOME bam while the pack
# was built from the raw STAR transcriptome bam (~29% PCR duplicates reached the pack).
set -euo pipefail
# GSE304796 (19-28z CAR-T) read processing, to the recipe the user specified:
#
#   umi_tools extract --bc-pattern=NNNNNNNNNNNN   12 nt 5' UMI -> read name
#   cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --discard-untrimmed -m 20 -M 40
#   bowtie2 --very-fast-local                     rRNA/tRNA depletion
#   STAR --alignEndsType EndToEnd --outFilterMultimapNmax 1
#   umi_tools dedup                               post-alignment
#   canonical chroms + snoRNA removal
#
# The 12 nt UMI was verified against the data before anything was run: stripping 12 nt from the 5' end
# and trimming this adapter gives mode 27 nt, 74.3% of kept reads in 26-34 nt and 89.7% inside the
# -m20/-M40 window -- a clean footprint distribution. GEO's own note ("UMIs were pruned using
# UMI-tools", Eclipsebio pipeline v1) confirms UMIs exist but never states the length; SRA carries the
# RAW reads, so they are still attached.
#
# ARM LABELS ARE THE USER'S, NOT GEO'S. GEO titles all six samples "Ribo-seq"; three are RNA. Verified
# empirically after trimming: the RPF trio shows 72-75% adapter read-through and 51.7% of inserts in
# 20-40 nt, the RNA trio 7-8% read-through and 92.6% full-length 72 nt. manifest.tsv carries the
# corrected `arm` column and everything here keys off that.
#
# The RNA arm runs the SAME UMI extraction and adapter trim -- both library types went through one
# adapter ligation -- but WITHOUT the -M 40 cap, which is footprint-specific and would discard nearly
# every RNA read.
#
#   sbatch --array=0-5 scripts/process_cart_gse304796.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
# umi_tools needs its OWN env: it does not support the riboseq env's python 3.13, and that env also
# contains a stale python3.5 umi_tools 0.2.3 egg whose broken shim shadows PATH. Absolute path only.
UMI=/private/groups/carpenterlab/emalekos/conda_envs/umitools/bin/umi_tools
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
CONTAM=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/contam_grch38_v49
SNORNA=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/snorna_grch38_v49
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
SALMONIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
EXT=$NEW/data/external/GSE304796_cart
# SNODEP_SUFFIX writes to a parallel output dir so a re-run with the snoRNA stage does not overwrite,
# or race, the first-pass outputs. Set it to "_snodep" for the depleted re-run.
OUT=$NEW/data/processed/GSE304796_cart${SNODEP_SUFFIX:-}
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
T=${SLURM_CPUS_PER_TASK:-16}
SCR=/data/tmp/emalekos/cart_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR" "$OUT"; trap 'rm -rf "$SCR"' EXIT

mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1"\t"$4}' "$EXT/manifest.tsv")
IFS=$'\t' read -r RUN ARM <<< "${RUNS[$SLURM_ARRAY_TASK_ID]}"
FQ=$EXT/${RUN}.fastq.gz
[[ -s "$FQ" ]] || { echo "missing $FQ" >&2; exit 1; }
for f in "$STARIDX" "$CONTAM.1.bt2" "$UMI"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
echo "run=$RUN arm=$ARM host=$(hostname) start=$(date -Iseconds)"

# ---- 1. UMI extract: 12 nt from the 5' end into the read name ----
"$UMI" extract --bc-pattern=NNNNNNNNNNNN --stdin="$FQ" --stdout="$SCR/umi.fq.gz" \
  --log="$OUT/${RUN}.umi_extract.log"

# ---- 2. adapter trim ----
# --discard-untrimmed applies to the RPF ARM ONLY. For a footprint, a read with no adapter never
# reached the end of its insert, so its length is a floor rather than a measurement and it must go
# (~30% of RPF reads, expected).
#
# For RNA it is actively WRONG and was corrected here after it fired: RNA inserts are longer than the
# 72 nt read, so only ~17% ever reach the adapter, and --discard-untrimmed threw away the other 83% as
# "failures" when they are ordinary full-length reads. Measured on the first pass: SRR34901749 kept
# 16.8% of 68.1 M reads, SRR34901751 kept 16.4% of 62.9 M. The -M 40 cap is likewise footprint-only.
LEN=(-m 20 -M 40); DISCARD=(--discard-untrimmed)
if [[ "$ARM" == "rna" ]]; then LEN=(-m 20); DISCARD=(); fi
"$ENV/bin/cutadapt" -a "$ADAPTER" "${DISCARD[@]}" "${LEN[@]}" -j "$T" \
  -o "$SCR/trim.fq.gz" "$SCR/umi.fq.gz" > "$OUT/${RUN}.cutadapt.log" 2>&1
grep -E "Reads written|Total reads processed|Reads that were too" "$OUT/${RUN}.cutadapt.log" | sed 's/^/  /'

# ---- 3. rRNA / tRNA depletion ----
"$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$CONTAM" -U "$SCR/trim.fq.gz" \
  --un-gz "$SCR/contam_clean.fq.gz" -S /dev/null 2> "$OUT/${RUN}.bowtie2_contam.log"
grep -E "overall alignment rate|reads; of these" "$OUT/${RUN}.bowtie2_contam.log" | sed 's/^/  contam: /'

# ---- 3b. snoRNA / scaRNA depletion ----
# A SECOND, separate stage rather than snoRNA folded into the contaminant index: it keeps the two
# depletion rates independently visible in the logs, and it avoided rebuilding an index that running
# jobs held open.
#
# The post-alignment BED filter at step 6 already removes snoRNA from the GENOME bam, so why deplete
# at read level? Two reasons, both load-bearing:
#   * the transcriptome BAM copied at the end never passes through that filter, so the reads feeding
#     P-sites and pack coverage still contained snoRNA;
#   * the salmon quant that DEFINES the transcript universe runs on FASTQ, upstream of every BAM
#     filter. Measured on this dataset: 43.1% of TPM sat in snoRNA host genes, which GENCODE types as
#     lncRNA so no biotype filter can reach them (15.3% in the poly(A) control GSE208041).
if [[ -s "$SNORNA.1.bt2" ]]; then
  "$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$SNORNA" -U "$SCR/contam_clean.fq.gz" \
    --un-gz "$SCR/clean.fq.gz" -S /dev/null 2> "$OUT/${RUN}.bowtie2_snorna.log"
  grep -E "overall alignment rate|reads; of these" "$OUT/${RUN}.bowtie2_snorna.log" | sed 's/^/  snoRNA: /'
else
  echo "ABORT: $SNORNA.1.bt2 missing -- build it with scripts/build_snorna_index.sbatch" >&2
  exit 1
fi

# ---- 4. STAR, end-to-end, unique alignments only ----
# EndToEnd forbids soft-clipping: a soft-clipped footprint end shifts the inferred P-site, which is
# the one thing this whole pipeline exists to measure.
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$SCR/clean.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --alignEndsType EndToEnd --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${RUN}."
cp "$SCR/${RUN}.Log.final.out" "$OUT/${RUN}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$OUT/${RUN}.STAR.Log.final.out" | sed 's/^/  /'

# ---- 5. dedup on the genome BAM, using the UMIs carried in the read names ----
"$ENV/bin/samtools" index -@ "$T" "$SCR/${RUN}.Aligned.sortedByCoord.out.bam"
"$UMI" dedup -I "$SCR/${RUN}.Aligned.sortedByCoord.out.bam" -S "$SCR/${RUN}.dedup.bam" \
  --log="$OUT/${RUN}.umi_dedup.log"
"$ENV/bin/samtools" index -@ "$T" "$SCR/${RUN}.dedup.bam"

# ---- 6. canonical chromosomes, minus snoRNA and the other post-filter classes ----
# chrM is dropped here too: the 13 chrM protein-coding mRNAs are excluded project-wide
# (feedback_riboseq_exclude_mito_genes), so letting them through only to drop them later invites a
# mismatch between what was aligned and what is scored.
CANON=$(printf "chr%s " {1..22} X Y)
"$ENV/bin/samtools" view -@ "$T" -b "$SCR/${RUN}.dedup.bam" $CANON > "$SCR/${RUN}.canon.bam"
"$ENV/bin/samtools" index -@ "$T" "$SCR/${RUN}.canon.bam"
POST=$NEW/data/ncrna_intervals/human_v49_postfilter.bed
if [[ -s "$POST" ]]; then
  "$ENV/bin/samtools" view -@ "$T" -b -L "$POST" -U "$OUT/${RUN}.final.bam" \
      "$SCR/${RUN}.canon.bam" > /dev/null
else
  echo "  WARNING: $POST missing -- no snoRNA/snRNA removal applied" >&2
  cp "$SCR/${RUN}.canon.bam" "$OUT/${RUN}.final.bam"
fi
"$ENV/bin/samtools" index -@ "$T" "$OUT/${RUN}.final.bam"
cp "$SCR/${RUN}.Aligned.toTranscriptome.out.bam" "$OUT/${RUN}.toTranscriptome.bam"

{ echo "run=$RUN arm=$ARM"
  for s in "$SCR/${RUN}.Aligned.sortedByCoord.out.bam" "$SCR/${RUN}.dedup.bam" \
           "$SCR/${RUN}.canon.bam" "$OUT/${RUN}.final.bam"; do
    printf "%-14s %s\n" "$(basename "$s" .bam)" "$("$ENV/bin/samtools" view -@ "$T" -c "$s")"
  done; } | tee "$OUT/${RUN}.read_counts.txt" | sed 's/^/  /'
# ---- 7. RNA arm only: salmon on the DEPLETED reads, for universe definition ----
# Run here, in the same job, rather than as a later pass: clean.fq.gz lives in node-local scratch and
# the EXIT trap deletes it, so any later job on another node would silently see nothing
# (/data/tmp is per-node). Quantifying the depleted reads is the whole point -- a salmon run on raw
# FASTQ is what let the snoRNA host genes into the universe in the first place.
if [[ "$ARM" == "rna" && -d "$SALMONIDX" ]]; then
  SQ=$OUT/salmon/${RUN}
  mkdir -p "$SQ"
  "$SALMON" quant -i "$SALMONIDX" -l A -r "$SCR/clean.fq.gz" -p "$T" --validateMappings -o "$SQ" \
    > "$OUT/${RUN}.salmon.log" 2>&1
  [[ -s "$SQ/quant.sf" ]] || { echo "ABORT: salmon produced no quant.sf" >&2; exit 1; }
  echo "  salmon: $(grep -o '\"percent_mapped\": *[0-9.]*' "$SQ/aux_info/meta_info.json" | head -1)"
fi
echo "done=$(date -Iseconds) -> $OUT/${RUN}.final.bam"
