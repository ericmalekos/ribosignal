#!/usr/bin/env bash
#SBATCH --job-name=xsp_rna
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/rna_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/rna_%A_%a.err
#
# RNA-seq alignment for the cross-species arm, at ONE multimap posture per invocation.
#
#   sbatch --array=0-N --export=ALL,SHEET=<rna_<sp>.tsv>,MM=10 scripts/xspecies/align_rna_xspecies.sbatch
#
# WHY THIS SCRIPT EXISTS RATHER THAN REUSING THE MOUSE TEMPLATES: both
# scripts/heldout/align_{janich,wang}_rna_mm10.sbatch are SINGLE-END ONLY -- one
# `--readFilesIn`, and a cutadapt call with no -A/-p. All primate and zebrafish RNA in
# this expansion is paired, so they cannot be used as-is. The PE pattern here follows
# scripts/heldout/screen_liver_polya_salmon.sbatch and scripts/rebuild_coverage_mm10.sbatch.
#
# RNA-SPECIFIC FLAGS, and why they differ from the Ribo recipe:
#   * NO -M/--maximum-length. RNA inserts are not ~30 nt footprints; capping length
#     here would discard nearly everything.
#   * NO --discard-untrimmed. Correct for RPF, wrong for RNA: applying it to a CAR-T
#     RNA arm once kept 16.4% of reads against 99.0% without.
#   * NO --alignEndsType EndToEnd. That exists for the Ribo arm because soft-clipping
#     shifts the inferred P-site; RNA coverage has no P-site to shift.
#
# The adapter comes from the samplesheet, which comes from the MEASURED probe. It is
# never defaulted.
set -euo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
SHEET=${SHEET:?set SHEET=<samplesheet.tsv>}
MM=${MM:?set MM=<multimap posture, e.g. 1 or 10>}
T=${SLURM_CPUS_PER_TASK:-12}

mapfile -t ROWS < <(awk -F'\t' 'NR>1' "$SHEET")
IFS=$'\t' read -r DS RUN R1 R2 ADAPTER SPECIES LAYOUT <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
# The samplesheet writes "-" for an absent R2 rather than an empty field: bash treats
# TAB as IFS whitespace, so consecutive tabs collapse into one delimiter and shift every
# later column left. A single-end row then parses SPECIES as "SINGLE".
[[ "$R2" == "-" ]] && R2=""
# Guard the shift explicitly rather than trusting the placeholder.
[[ "$LAYOUT" == "PAIRED" || "$LAYOUT" == "SINGLE" ]] || {
  echo "samplesheet row parsed wrong: layout='$LAYOUT' species='$SPECIES'" >&2; exit 1; }

SPTSV=$NEW/scripts/xspecies/species.tsv
ANNVER=$(awk -F'\t' -v s="$SPECIES" '!/^#/ && $1==s {print $2; exit}' "$SPTSV")
[[ -n "$ANNVER" ]] || { echo "species '$SPECIES' has no row in $SPTSV" >&2; exit 1; }
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_${SPECIES}_${ANNVER}
[[ -s "$IDX/SA" ]] || { echo "missing STAR index $IDX" >&2; exit 1; }

OUTBAMDIR=$NEW/data/xspecies_rna_bam_mm${MM}/$DS
OUTCOVDIR=$NEW/data/xspecies_rna_coverage_mm${MM}/$DS
mkdir -p "$OUTBAMDIR" "$OUTCOVDIR"
OUTBAM=$OUTBAMDIR/${RUN}.Aligned.toTranscriptome.out.bam
OUTCOV=$OUTCOVDIR/${RUN}_coverage.hd5

[[ -s "$OUTCOV" ]] && { echo "SKIP $RUN (coverage exists)"; exit 0; }

SCR=/data/tmp/emalekos/xsprna_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT

echo "dataset=$DS run=$RUN species=$SPECIES mm=$MM layout=$LAYOUT adapter=$ADAPTER"
echo "  host=$(hostname) cpus=$T start=$(date -Iseconds)"
for f in "$R1" "$IDX"; do [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }; done

# ---- trim -------------------------------------------------------------------
if [[ "$LAYOUT" == "PAIRED" ]]; then
  [[ -e "$R2" ]] || { echo "PAIRED but R2 missing: $R2" >&2; exit 1; }
  if [[ "$ADAPTER" == "--trim-n" ]]; then
    "$ENV/bin/cutadapt" --trim-n -m 20 -j "$T" \
      -o "$SCR/t1.fq.gz" -p "$SCR/t2.fq.gz" "$R1" "$R2" > "$OUTBAMDIR/${RUN}.cutadapt.log" 2>&1
  else
    "$ENV/bin/cutadapt" -a "$ADAPTER" -A "$ADAPTER" -m 20 -j "$T" \
      -o "$SCR/t1.fq.gz" -p "$SCR/t2.fq.gz" "$R1" "$R2" > "$OUTBAMDIR/${RUN}.cutadapt.log" 2>&1
  fi
  READFILES=("$SCR/t1.fq.gz" "$SCR/t2.fq.gz")
else
  if [[ "$ADAPTER" == "--trim-n" ]]; then
    "$ENV/bin/cutadapt" --trim-n -m 20 -j "$T" \
      -o "$SCR/t1.fq.gz" "$R1" > "$OUTBAMDIR/${RUN}.cutadapt.log" 2>&1
  else
    "$ENV/bin/cutadapt" -a "$ADAPTER" -m 20 -j "$T" \
      -o "$SCR/t1.fq.gz" "$R1" > "$OUTBAMDIR/${RUN}.cutadapt.log" 2>&1
  fi
  READFILES=("$SCR/t1.fq.gz")
fi
grep -E "Total reads processed|Total read pairs processed|Reads written|Pairs written" \
  "$OUTBAMDIR/${RUN}.cutadapt.log" | sed 's/^/  /'

# ---- align ------------------------------------------------------------------
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "${READFILES[@]}" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax "$MM" --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/${RUN}."
cp "$SCR/${RUN}.Log.final.out" "$OUTBAMDIR/${RUN}.STAR.Log.final.out"
grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci" \
  "$OUTBAMDIR/${RUN}.STAR.Log.final.out" | sed 's/^/  /'

TXBAM="$SCR/${RUN}.Aligned.toTranscriptome.out.bam"
# quickcheck, not a header grep: a BAM truncated mid-write keeps a VALID header, so
# `view -H | grep SO:coordinate` accepts it. This cost 92 minutes of RiboCode once.
"$ENV/bin/samtools" quickcheck "$TXBAM" || { echo "ABORT: $TXBAM failed quickcheck" >&2; exit 2; }

# ---- per-nt coverage --------------------------------------------------------
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TXBAM" "$OUTCOV" "$RUN"
[[ -s "$OUTCOV" ]] || { echo "ABORT: no coverage written" >&2; exit 3; }
mv "$TXBAM" "$OUTBAM"

echo "  coverage $(du -h "$OUTCOV" | cut -f1)   bam $(du -h "$OUTBAM" | cut -f1)"
echo "=== $RUN DONE=$(date -Iseconds)"
