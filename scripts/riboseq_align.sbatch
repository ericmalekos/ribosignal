#!/usr/bin/env bash
#SBATCH --job-name=ribo_align
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/ribo_align_%A_%a.out
set -euo pipefail
# THE canonical Ribo-seq aligner. One script for every dataset and both species.
# Replaces align_janich_ribo / align_wang_ribo / align_chothani_refetch / the STAR half of
# process_thp1 and process_cart_gse304796. Per docs/PIPELINE_POLICY.md, no new align_*.sbatch.
#
# FIXED FOR EVERY RUN, not tunable per dataset:
#   --alignEndsType EndToEnd        soft-clipping shifts the inferred P-site, which is the exact
#                                   quantity the model predicts. Measured on HUVEC against its
#                                   reference pack, alignment recipe the only variable:
#                                   EndToEnd r=0.935 vs soft-clipped r=0.583.
#   --outFilterMultimapNmax 1       RiboCode does NOT drop multi-mappers itself; each record
#                                   contributes a P-site to every transcript it points at.
#   --quantMode TranscriptomeSAM    downstream P-site calling is transcriptome-coordinate.
#   then: coordinate-sort + index + filter_tx_heldout (ncRNA drop + cross-gene read drop).
#
# LEFT AT STAR DEFAULTS, deliberately and with evidence. A 13-arm sweep on HUVEC
# (results/star_param_sweep/) varied outFilterMismatchNoverLmax, outFilterMatchNminOverLread,
# outFilterScoreMinOverLread and seedSearchStartLmax. In-frame fraction was 55.95-56.05% in EVERY
# arm -- a 0.10 point spread across a 3x change in mismatch rate. Tightening costs up to 5.3% of
# reads and buys no frame purity, and would widen divergence from the training packs, which were
# built under these defaults. seedSearchStartLmax is inert at these read lengths (20 and 30 give
# byte-identical counts to the default 50, because the default already exceeds the read length).
#
# PER-DATASET VARIATION comes only from the samplesheet: adapter, UMI, and species (which selects
# index, ncRNA list and tx->gene map). Nothing else.
#
# Samplesheet TSV: dataset  run  fastq  adapter  species  [umi]  [cutadapt_extra]
#   sbatch --array=0-N --export=ALL,SHEET=<tsv>,OUTROOT=<dir> scripts/riboseq_align.sbatch
#
# Columns 6-7 are OPTIONAL and default to empty; 5-column samplesheets keep working unchanged.
#   umi             a umi_tools bc-pattern (e.g. NNNNNNNNNNNN for GSE304796's 12 nt 5' UMI), or `none`
#   cutadapt_extra  extra cutadapt flags, space-separated (e.g. --discard-untrimmed)
#
# `cutadapt_extra` exists so the four datasets this project already has can ALL be expressed as
# samplesheet rows rather than as branches in a script. Before it, the CAR-T RPF arm could not be:
# it needs --discard-untrimmed, because with a 12 nt UMI on a 72 nt read only ~17% of reads reach the
# adapter, and a footprint read with no adapter was never a complete footprint. That is a property of
# the library, which is exactly what belongs in a samplesheet. Keep it to genuine library properties;
# it is not a hatch for re-tuning the fixed STAR filters.
#
# UMI ORDER MATTERS, and getting it wrong is silent. dedup runs AFTER filter_tx_heldout, not before:
#   filter_tx_heldout  needs QUERY-GROUPED input (cross-gene detection compares a read's alignments)
#   umi_tools dedup    needs COORDINATE-SORTED + INDEXED input
# The filter consumes STAR's read-ordered output and leaves a sorted, indexed BAM behind, so the two
# requirements chain in exactly one workable order. Deduping first would hand the filter a sorted BAM,
# which it now refuses outright -- before that guard existed it silently reported ~0.30% dropped
# instead of the true 9-22%.
#
# DEDUP THE BAM THE PACK ACTUALLY READS. process_cart_gse304796.sbatch deduplicated the GENOME bam and
# then built its pack from the raw STAR TRANSCRIPTOME bam, so ~29% PCR duplicates went into the pack
# while the dedup step looked like it had run. Here there is only one BAM and the filter and dedup
# both act on it.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
RPY=$ENV/bin/python3
SHEET=${SHEET:?set SHEET=<samplesheet.tsv>}
OUTROOT=${OUTROOT:?set OUTROOT=<output dir>}
# Multimap posture. DEFAULTS TO 1 and must stay that way: every pack, P-site set and model
# target in this project is built on uniquely-mapped reads, per the project-wide Ribo-seq
# rule (RiboCode and Ribo-TISH do not drop multi-mappers themselves, so keeping them
# silently inflates P-site counts and breaks the periodicity test).
#
# MM=25 exists ONLY as the diagnostic arm that measures what mm1 discards. Its BAMs are
# NOT a valid model target. When comparing the two postures, compare TOTAL RECORDS PER
# TRANSCRIPT via bamlib.multimap_cost(); do NOT split on NH, because --quantMode
# TranscriptomeSAM projects one genomic alignment onto every compatible isoform, so a
# genomically unique read routinely carries NH:i:17. An NH split once reported "99.8% of
# reads discarded" against STAR's own 30%.
MM=${MM:-1}
[[ "$MM" =~ ^[0-9]+$ ]] || { echo "MM must be an integer, got '$MM'" >&2; exit 1; }
T=${SLURM_CPUS_PER_TASK:-12}

mapfile -t ROWS < <(awk -F'\t' 'NR>1' "$SHEET")
IFS=$'\t' read -r DS RUN FQ ADAPTER SPECIES UMIPAT CUTEXTRA <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
UMIPAT=${UMIPAT:-none}
# Word-split deliberately: this column carries flags, not a filename.
read -r -a CUTX <<< "${CUTEXTRA:-}"
# umi_tools needs its OWN env: it does not support the riboseq env's python 3.13, and that env also
# carries a stale python3.5 umi_tools 0.2.3 egg whose broken shim shadows PATH. Absolute path only.
UMITOOLS=/private/groups/carpenterlab/emalekos/conda_envs/umitools/bin/umi_tools
if [[ "$UMIPAT" != "none" ]]; then
  [[ -x "$UMITOOLS" ]] || { echo "umi=$UMIPAT requested but $UMITOOLS is missing" >&2; exit 1; }
  [[ "$UMIPAT" =~ ^[NCGATX]+$ ]] || { echo "umi pattern '$UMIPAT' is not a bc-pattern" >&2; exit 1; }
fi
# human and mouse keep their historical hardcoded paths: those packs were built against
# these exact files and re-pointing them is not a refactor, it is a re-validation.
# Every other species resolves through scripts/xspecies/species.tsv, so adding an eighth
# needs a row in that table and no edit here.
case "$SPECIES" in
  mouse) IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
         NCRNA=$NEW/data/ncrna_tx_mouse_vM38.txt; T2G=$NEW/data/tx_to_gene_mouse_vM38.tsv ;;
  human) IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
         NCRNA=$NEW/data/ncrna_tx_human_v49.txt;  T2G=$NEW/data/tx_to_gene_human_v49.tsv ;;
  *)
    SPTSV=$NEW/scripts/xspecies/species.tsv
    ANNVER=$(awk -F'\t' -v s="$SPECIES" '!/^#/ && $1==s {print $2; exit}' "$SPTSV" 2>/dev/null)
    [[ -n "$ANNVER" ]] || { echo "unknown species '$SPECIES' (not hardcoded, and no row in $SPTSV)" >&2; exit 1; }
    IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_${SPECIES}_${ANNVER}
    NCRNA=$NEW/data/annot/$SPECIES/ncrna_tx.txt
    T2G=$NEW/data/annot/$SPECIES/tx_to_gene.tsv
    # An empty ncRNA list is the exact silent failure normalize_annotation.py exists to
    # prevent: the drop step would run, remove nothing, and report success.
    [[ -s "$NCRNA" ]] && [[ $(wc -l < "$NCRNA") -ge 10 ]] || {
      echo "ncRNA drop list $NCRNA is missing or implausibly small; run scripts/xspecies/build_species_refs.sbatch" >&2
      exit 1; }
    ;;
esac
BAMDIR=$OUTROOT/$DS
mkdir -p "$BAMDIR"
OUTBAM=$BAMDIR/${RUN}.Aligned.toTranscriptome.out.bam
for f in "$FQ" "$IDX" "$NCRNA" "$T2G"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
[[ -s "$OUTBAM" && -s "$OUTBAM.bai" ]] && { echo "SKIP $RUN (done)"; exit 0; }
SCR=/data/tmp/emalekos/ribal_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "dataset=$DS run=$RUN species=$SPECIES adapter=$ADAPTER umi=$UMIPAT mm=$MM host=$(hostname) start=$(date -Iseconds)"

# UMI extract BEFORE trimming: the UMI sits at the 5' end of the raw read, so it must come off before
# cutadapt applies the -m/-M length window, or the window is measured on insert+UMI.
CUTIN=$FQ
if [[ "$UMIPAT" != "none" ]]; then
  "$UMITOOLS" extract --bc-pattern="$UMIPAT" --stdin "$FQ" --stdout "$SCR/umi.fq.gz" \
    --log "$BAMDIR/${RUN}.umi_extract.log"
  CUTIN=$SCR/umi.fq.gz
fi

# Adapter: a literal sequence, or the token --trim-n for N-padded fixed-length reads.
if [[ "$ADAPTER" == "--trim-n" ]]; then
  "$ENV/bin/cutadapt" --trim-n "${CUTX[@]}" -m 20 -M 40 -j "$T" -o "$SCR/t.fq.gz" "$CUTIN" > "$BAMDIR/${RUN}.cutadapt.log" 2>&1
else
  "$ENV/bin/cutadapt" -a "$ADAPTER" "${CUTX[@]}" -m 20 -M 40 -j "$T" -o "$SCR/t.fq.gz" "$CUTIN" > "$BAMDIR/${RUN}.cutadapt.log" 2>&1
fi
grep -E "Total reads processed|Reads written" "$BAMDIR/${RUN}.cutadapt.log" | sed 's/^/  /'

"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/t.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --alignEndsType EndToEnd --outFilterMultimapNmax "$MM" \
  --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${RUN}."
cp "$SCR/${RUN}.Log.final.out" "$BAMDIR/${RUN}.STAR.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$BAMDIR/${RUN}.STAR.Log.final.out" | sed 's/^/  /'

TB=$SCR/${RUN}.Aligned.toTranscriptome.out.bam
[[ -s "$TB" ]] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
# TranscriptomeSAM is read-ordered; sort here, in the same job. filter_tx_heldout needs query
# grouping (which STAR provides) and re-sorts + indexes in place afterwards.
before=$("$ENV/bin/samtools" view -@ "$T" -c "$TB")
cp "$TB" "$OUTBAM"
"$RPY" "$NEW/scripts/heldout/filter_tx_heldout.py" "$OUTBAM" "$NCRNA" "$T2G"
after=$("$ENV/bin/samtools" view -@ "$T" -c "$OUTBAM")
echo "  filter_tx_heldout: $before -> $after records ($(awk -v a=$before -v b=$after 'BEGIN{printf "%.2f",100*(a-b)/a}')% dropped)"
# quickcheck as well as the header -- a truncated BAM keeps a valid header (see 2026-08-21).
"$ENV/bin/samtools" quickcheck "$OUTBAM" 2>/dev/null \
  || { echo "$OUTBAM is TRUNCATED (quickcheck failed)" >&2; exit 1; }
"$ENV/bin/samtools" view -H "$OUTBAM" | sed -n 1p | grep -q SO:coordinate \
  || { echo "ABORT: output BAM is not coordinate-sorted" >&2; exit 1; }
[[ -s "$OUTBAM.bai" ]] || "$ENV/bin/samtools" index -@ "$T" "$OUTBAM"

# UMI dedup LAST, on the coordinate-sorted+indexed BAM the filter just produced -- and on the same
# file the pack will read. Each transcript is its own reference here, so with posture A (within-gene
# isoform multimapping retained) a read present on three isoforms is deduplicated independently on
# each, which is what per-transcript P-site counting wants.
if [[ "$UMIPAT" != "none" ]]; then
  predup=$("$ENV/bin/samtools" view -@ "$T" -c "$OUTBAM")
  "$UMITOOLS" dedup --stdin "$OUTBAM" --stdout "$SCR/dedup.bam" \
    --log "$BAMDIR/${RUN}.umi_dedup.log"
  mv "$SCR/dedup.bam" "$OUTBAM"
  "$ENV/bin/samtools" index -@ "$T" "$OUTBAM"
  postdup=$("$ENV/bin/samtools" view -@ "$T" -c "$OUTBAM")
  echo "  umi_tools dedup: $predup -> $postdup records ($(awk -v a=$predup -v b=$postdup 'BEGIN{printf "%.2f",100*(a-b)/a}')% duplicates)"
fi
echo "done=$(date -Iseconds) -> $OUTBAM"
