#!/usr/bin/env bash
#SBATCH --job-name=thp1_proc
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/thp1_proc_%A_%a.out
set -euo pipefail
# Both THP-1 datasets. Same skeleton as the CAR-T pipeline (bowtie2 contaminant depletion -> STAR
# EndToEnd mm1 -> canonical chroms minus snoRNA), but the READ-PREP step differs per dataset because
# the two were built with completely different chemistries. Each recipe was derived from the data, not
# assumed -- see DATA_PROVENANCE.md.
#
#   GSE208041 RPF   N-PADDED. 99.5% of reads are exactly 35 nt with 0% adapter, which reads as
#                   "already trimmed". The submitter trimmed AND padded back to fixed length with N
#                   (34% of reads end in 8 Ns). Aligning as-is charges every pad N as a mismatch
#                   under EndToEnd. Strip the trailing N run -> 84.4% of inserts land in 26-34 nt,
#                   mode 28. No UMI: the insert is already footprint-length, so a UMI would have to
#                   be inside the footprint, which it is not.
#   GSE208041 RNA   Paired 2x101, poly(A)-selected (TruSeq Stranded mRNA), no adapter read-through.
#                   Nothing to trim.
#   GSE39561        POLY(A)-TAILED, 2012 Ingolia-era protocol: the footprint is poly(A)-tailed before
#                   reverse transcription, so there is no ligated adapter to find -- 0% TruSeq, and
#                   the dominant 12-mer at read position 29 is AAAAAAAAAAAA. Cut at the first run of
#                   >=6 A -> 69.6% in 26-34 nt. Standard TruSeq trimming leaves the poly(A) attached.
#
# NO umi_tools ANYWHERE HERE. GSE208041's submitter already removed UMIs (the stripped insert is
# footprint-length, leaving no room for one) and GSE39561 predates UMIs entirely. Running dedup on
# these would collapse genuinely independent footprints that share a start position.
#
#   sbatch --array=0-11%4 scripts/process_thp1.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
CONTAM=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/contam_grch38_v49
SNORNA=/private/groups/carpenterlab/emalekos/Bowtie2_indexes/snorna_grch38_v49
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
SALMONIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
POST=$NEW/data/ncrna_intervals/human_v49_postfilter.bed
T=${SLURM_CPUS_PER_TASK:-16}

# Flat task table across both datasets: dataset<TAB>run<TAB>arm
mapfile -t TASKS < <(
  for ds in GSE208041_thp1 GSE39561_thp1; do
    awk -F'\t' -v d="$ds" 'NR>1{print d"\t"$1"\t"$4}' "$NEW/data/external/$ds/manifest.tsv"
  done)
IFS=$'\t' read -r DS RUN ARM <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"
EXT=$NEW/data/external/$DS
# SNODEP_SUFFIX writes to a parallel output dir so the snoRNA-depleted re-run neither overwrites nor
# races the first-pass outputs. Set it to "_snodep" for the depleted re-run.
OUT=$NEW/data/processed/${DS}${SNODEP_SUFFIX:-}
SCR=/data/tmp/emalekos/thp1_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR" "$OUT"; trap 'rm -rf "$SCR"' EXIT
echo "dataset=$DS run=$RUN arm=$ARM host=$(hostname) start=$(date -Iseconds)"

R1=$EXT/${RUN}.fastq.gz; [[ -s "$R1" ]] || R1=$EXT/${RUN}_1.fastq.gz
R2=$EXT/${RUN}_2.fastq.gz
[[ -s "$R1" ]] || { echo "missing reads for $RUN" >&2; exit 1; }
PAIRED=0; [[ "$ARM" == "rna" && -s "$R2" ]] && PAIRED=1

# ---- 1. dataset-specific read prep ----
if [[ "$DS" == "GSE208041_thp1" && "$ARM" == "ribo" ]]; then
  # Strip the trailing N pad, then keep footprint-length inserts. cutadapt's --trim-n removes N from
  # BOTH ends; only the 3' pad is artificial here, but a leading N is an uncalled base and is no loss.
  "$ENV/bin/cutadapt" --trim-n -m 20 -M 40 -j "$T" \
    -o "$SCR/prep.fq.gz" "$R1" > "$OUT/${RUN}.cutadapt.log" 2>&1
elif [[ "$DS" == "GSE39561_thp1" ]]; then
  # Poly(A) tail, not a ligated adapter. -a "A{6}" cuts at the first run of 6+ A.
  "$ENV/bin/cutadapt" -a "A{6}" -m 20 -M 40 -j "$T" \
    -o "$SCR/prep.fq.gz" "$R1" > "$OUT/${RUN}.cutadapt.log" 2>&1
elif [[ "$PAIRED" -eq 1 ]]; then
  "$ENV/bin/cutadapt" -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 25 -j "$T" \
    -o "$SCR/prep_1.fq.gz" -p "$SCR/prep_2.fq.gz" "$R1" "$R2" > "$OUT/${RUN}.cutadapt.log" 2>&1
else
  "$ENV/bin/cutadapt" -a AGATCGGAAGAGC -m 25 -j "$T" \
    -o "$SCR/prep.fq.gz" "$R1" > "$OUT/${RUN}.cutadapt.log" 2>&1
fi
grep -E "Total reads processed|Total read pairs|Reads written|pairs written" "$OUT/${RUN}.cutadapt.log" | sed 's/^/  /'

# ---- 2. contaminant depletion (Ribo arm only; RNA keeps its ncRNA for honest quantification) ----
if [[ "$ARM" == "ribo" ]]; then
  "$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$CONTAM" -U "$SCR/prep.fq.gz" \
    --un-gz "$SCR/contam_clean.fq.gz" -S /dev/null 2> "$OUT/${RUN}.bowtie2_contam.log"
  grep -E "overall alignment rate" "$OUT/${RUN}.bowtie2_contam.log" | sed 's/^/  contam: /'
  READS=("$SCR/contam_clean.fq.gz")
elif [[ "$PAIRED" -eq 1 ]]; then
  READS=("$SCR/prep_1.fq.gz" "$SCR/prep_2.fq.gz")
else
  READS=("$SCR/prep.fq.gz")
fi

# ---- 2b. snoRNA / scaRNA depletion, BOTH arms ----
# Applied to the RNA arm as well as Ribo, which is the change: snoRNA host genes are typed lncRNA in
# GENCODE, so no biotype filter downstream can reach them, and the salmon quant that defines the
# transcript universe reads FASTQ upstream of every BAM filter. Paired input needs --un-conc-gz, not
# --un-gz; bowtie2 writes .1/.2 from the %-substituted pattern.
if [[ -s "$SNORNA.1.bt2" ]]; then
  if [[ "${#READS[@]}" -eq 2 ]]; then
    "$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$SNORNA" \
      -1 "${READS[0]}" -2 "${READS[1]}" \
      --un-conc-gz "$SCR/sno_clean_%.fq.gz" -S /dev/null 2> "$OUT/${RUN}.bowtie2_snorna.log"
    READS=("$SCR/sno_clean_1.fq.gz" "$SCR/sno_clean_2.fq.gz")
  else
    "$ENV/bin/bowtie2" --very-fast-local -p "$T" -x "$SNORNA" -U "${READS[0]}" \
      --un-gz "$SCR/sno_clean.fq.gz" -S /dev/null 2> "$OUT/${RUN}.bowtie2_snorna.log"
    READS=("$SCR/sno_clean.fq.gz")
  fi
  grep -E "overall alignment rate" "$OUT/${RUN}.bowtie2_snorna.log" | sed 's/^/  snoRNA: /'
  for f in "${READS[@]}"; do
    [[ -s "$f" ]] || { echo "ABORT: bowtie2 wrote no $f" >&2; exit 1; }
  done
else
  echo "ABORT: $SNORNA.1.bt2 missing -- build it with scripts/build_snorna_index.sbatch" >&2
  exit 1
fi

# ---- 3. STAR ----
# EndToEnd for the Ribo arm: soft-clipping a footprint end shifts the inferred P-site. RNA keeps the
# default Local, where clipping a stray adapter base is harmless and dropping the read is not.
ENDS=(--alignEndsType EndToEnd); [[ "$ARM" == "rna" ]] && ENDS=()
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "${READS[@]}" --readFilesCommand zcat \
  --runThreadN "$T" "${ENDS[@]}" --outFilterMultimapNmax 1 \
  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${RUN}."
cp "$SCR/${RUN}.Log.final.out" "$OUT/${RUN}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$OUT/${RUN}.STAR.Log.final.out" | sed 's/^/  /'

# ---- 4. canonical chroms, minus the post-filter ncRNA classes (Ribo arm) ----
BAM=$SCR/${RUN}.Aligned.sortedByCoord.out.bam
"$ENV/bin/samtools" index -@ "$T" "$BAM"
CANON=$(printf "chr%s " {1..22} X Y)
"$ENV/bin/samtools" view -@ "$T" -b "$BAM" $CANON > "$SCR/${RUN}.canon.bam"
"$ENV/bin/samtools" index -@ "$T" "$SCR/${RUN}.canon.bam"
if [[ "$ARM" == "ribo" && -s "$POST" ]]; then
  "$ENV/bin/samtools" view -@ "$T" -b -L "$POST" -U "$OUT/${RUN}.final.bam" \
      "$SCR/${RUN}.canon.bam" > /dev/null
else
  cp "$SCR/${RUN}.canon.bam" "$OUT/${RUN}.final.bam"
fi
"$ENV/bin/samtools" index -@ "$T" "$OUT/${RUN}.final.bam"
cp "$SCR/${RUN}.Aligned.toTranscriptome.out.bam" "$OUT/${RUN}.toTranscriptome.bam"

{ echo "dataset=$DS run=$RUN arm=$ARM"
  for s in "$BAM" "$SCR/${RUN}.canon.bam" "$OUT/${RUN}.final.bam"; do
    printf "%-34s %s\n" "$(basename "$s" .bam)" "$("$ENV/bin/samtools" view -@ "$T" -c "$s")"
  done; } | tee "$OUT/${RUN}.read_counts.txt" | sed 's/^/  /'
# ---- 5. RNA arm only: salmon on the DEPLETED reads, for universe definition ----
# In-job, not a later pass: these reads live in node-local /data/tmp and the EXIT trap removes them,
# so a job on another node would silently find nothing. Quantifying the depleted reads is the point --
# salmon on raw FASTQ is what admitted the snoRNA host genes to the universe.
if [[ "$ARM" == "rna" && -d "$SALMONIDX" ]]; then
  SQ=$OUT/salmon/${RUN}
  mkdir -p "$SQ"
  if [[ "${#READS[@]}" -eq 2 ]]; then
    "$SALMON" quant -i "$SALMONIDX" -l A -1 "${READS[0]}" -2 "${READS[1]}" \
      -p "$T" --validateMappings -o "$SQ" > "$OUT/${RUN}.salmon.log" 2>&1
  else
    "$SALMON" quant -i "$SALMONIDX" -l A -r "${READS[0]}" \
      -p "$T" --validateMappings -o "$SQ" > "$OUT/${RUN}.salmon.log" 2>&1
  fi
  [[ -s "$SQ/quant.sf" ]] || { echo "ABORT: salmon produced no quant.sf" >&2; exit 1; }
  echo "  salmon: $(grep -o '\"percent_mapped\": *[0-9.]*' "$SQ/aux_info/meta_info.json" | head -1)"
fi
echo "done=$(date -Iseconds) -> $OUT/${RUN}.final.bam"
