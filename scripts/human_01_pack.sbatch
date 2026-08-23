#!/usr/bin/env bash
#SBATCH --job-name=hs_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hs_pack_%A_%a.out
set -euo pipefail
# Pool each human dataset from BAMs onto its own universe, and build its ORF track.
#
#   task 0  GSE208041 (THP-1)   5 RPF + 4 RNA, own universe
#   task 1  GSE39561  (THP-1)   3 RPF + GSE208041's 4 RNA  <- BORROWED RNA, see below
#   task 2  GSE304796 (CAR-T)   3 RPF + 3 RNA, own universe (ncRNA_host genes excluded)
#
# BORROWED RNA FOR GSE39561. That dataset has three RPF runs and no RNA-seq at all, so it cannot
# supply the coverage input the model requires, nor an expression-restricted universe. It uses
# GSE208041's RNA for both. The two are the same cell line (THP-1), different study. The
# per-dataset-universe rule was downgraded from a hard bar to a warning on 2026-08-10 specifically to
# allow this, on condition it is labelled everywhere -- the universe TSV carries a `# BORROWED
# UNIVERSE` header and this pack's provenance records it.
#
# The mouse 3x3 factorial is the direct evidence that this is tolerable: swapping in an unrelated
# experiment's RNA-seq moved ORF-call F1 by +0.0030 on average (range -0.0057 to +0.0079) across 12
# model x reference x arm combinations. Borrowing RNA from the SAME cell line should cost less than
# that. It is still a caveat to report, not a free pass.
#
#   sbatch --array=0-2 scripts/human_01_pack.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
# rnaseq_coverage.py needs pysam, which cas12a lacks; riboseq has 0.23.3. Getting this wrong kills the
# coverage stage AFTER the expensive Ribo stage has already run.
PYSAM=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
ANNOT=$NEW/data/human_ribocode_annot_primary
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
REF=$NEW/data/heldout_refs
P=$NEW/data/processed

DS=(gse208041 gse39561 cart)
D=${DS[$SLURM_ARRAY_TASK_ID]}

A=$P/GSE208041_thp1_snodep
case "$D" in
  gse208041)
    RIBO=($A/SRR20106164.toTranscriptome.bam $A/SRR20106165.toTranscriptome.bam
          $A/SRR20106166.toTranscriptome.bam $A/SRR20106167.toTranscriptome.bam
          $A/SRR20106168.toTranscriptome.bam)
    RNA=($A/SRR21228005.toTranscriptome.bam $A/SRR21228006.toTranscriptome.bam
         $A/SRR21228015.toTranscriptome.bam $A/SRR21228016.toTranscriptome.bam)
    UNI=$REF/human_thp1_gse208041_tx.txt; FA=$REF/human_thp1_gse208041_universe.fa ;;
  gse39561)
    B=$P/GSE39561_thp1_snodep
    RIBO=($B/SRR525273.toTranscriptome.bam $B/SRR525274.toTranscriptome.bam
          $B/SRR525275.toTranscriptome.bam)
    RNA=($A/SRR21228005.toTranscriptome.bam $A/SRR21228006.toTranscriptome.bam
         $A/SRR21228015.toTranscriptome.bam $A/SRR21228016.toTranscriptome.bam)
    UNI=$REF/human_thp1_gse39561_tx.txt; FA=$REF/human_thp1_gse39561_universe.fa ;;
  cart)
    C=$P/GSE304796_cart_snodep
    RIBO=($C/SRR34901743.toTranscriptome.bam $C/SRR34901744.toTranscriptome.bam
          $C/SRR34901745.toTranscriptome.bam)
    RNA=($C/SRR34901749.toTranscriptome.bam $C/SRR34901750.toTranscriptome.bam
         $C/SRR34901751.toTranscriptome.bam)
    UNI=$REF/human_cart_gse304796_tx.txt; FA=$REF/human_cart_gse304796_universe.fa ;;
esac
OUT=$NEW/data/packed_heldout_human_$D
mkdir -p "$OUT"

for f in "${RIBO[@]}" "${RNA[@]}" "$UNI" "$FA" "$ANNOT/transcripts.pickle"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# prepare_from_bams calls samtools index; TranscriptomeSAM output is read-ordered until sorted.
S=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
for f in "${RIBO[@]}" "${RNA[@]}"; do
  # quickcheck as well as the header -- a truncated BAM keeps a valid header (see 2026-08-21).
  "$S" quickcheck "$f" 2>/dev/null || { echo "$f is TRUNCATED" >&2; exit 1; }
  [[ "$("$S" view -H "$f" 2>/dev/null | sed -n 1p)" == *SO:coordinate* ]] \
    || { echo "ABORT: $f is not coordinate-sorted -- run sort_human_tx_bams.sbatch first" >&2; exit 1; }
done
echo "dataset=$D ribo=${#RIBO[@]} rna=${#RNA[@]} universe=$(wc -l < "$UNI") host=$(hostname) start=$(date -Iseconds)"

# REUSE COVERAGE ACROSS DATASETS THAT SHARE RNA. GSE39561 borrows GSE208041's four RNA BAMs and its
# universe is byte-identical (md5 38ed6010a757deac80060eece323c68d), so the per-sample coverage hd5
# are the SAME FILES computed twice. On the first attempt both tasks derived them in parallel: ~2.5 h
# per sample on a paired-end RNA transcriptome BAM (1.48 BILLION records for an 82 M-read sample --
# TranscriptomeSAM emits one record per compatible isoform, the ~18x multimap inflation), so the
# duplication was about 5 h of wasted compute.
#
# HS_REUSE_COVERAGE_FROM=<pack dir> takes that dataset's coverage hd5 instead of re-deriving them.
# Only valid when the universes are identical -- asserted below, not assumed.
REUSE=${HS_REUSE_COVERAGE_FROM:-}
if [[ -n "$REUSE" ]]; then
  mapfile -t RC < <(ls "$REUSE"/_work/coverage/*_coverage.hd5 2>/dev/null | sort)
  [[ ${#RC[@]} -eq ${#RNA[@]} ]] || { echo "ABORT: $REUSE has ${#RC[@]} coverage hd5, need ${#RNA[@]}" >&2; exit 1; }
  # The donor may be an INCOMPLETE pack -- coverage hd5 exist but tx_order.txt does not, because the
  # run that produced them failed later. Under `set -e` an unguarded `x=$(md5sum missing)` aborts the
  # script with NO message, which is exactly how this failed the first time (dead in 5 s, empty log).
  md5a=$(md5sum "$UNI" | cut -d" " -f1)
  md5b=$(md5sum "$REUSE/tx_order.txt" 2>/dev/null | cut -d" " -f1 || true)
  if [[ -n "$md5b" ]]; then
    [[ "$md5a" == "$md5b" ]] || { echo "ABORT: universe mismatch -- $UNI ($md5a) vs \
$REUSE/tx_order.txt ($md5b); coverage is not reusable" >&2; exit 1; }
    echo "  reusing ${#RC[@]} coverage hd5 from $REUSE (tx_order md5 matches)"
  elif [[ -n "${HS_REUSE_UNIVERSE:-}" ]]; then
    md5c=$(md5sum "$HS_REUSE_UNIVERSE" | cut -d" " -f1)
    [[ "$md5a" == "$md5c" ]] || { echo "ABORT: universe mismatch -- $UNI ($md5a) vs \
$HS_REUSE_UNIVERSE ($md5c)" >&2; exit 1; }
    echo "  reusing ${#RC[@]} coverage hd5 from $REUSE (donor pack incomplete; universe md5 matches \
via HS_REUSE_UNIVERSE)"
  else
    echo "ABORT: $REUSE has no tx_order.txt and HS_REUSE_UNIVERSE is unset, so the coverage axis \
cannot be verified against $UNI. Pass HS_REUSE_UNIVERSE=<donor tx list> to assert they match." >&2
    exit 1
  fi
  RNA_ARGS=(--rna-coverage "${RC[@]}")
else
  RNA_ARGS=(--rna-bam "${RNA[@]}")
fi

# Resumable: the Ribo stage (metaplots + RiboCode) writes its hd5 before coverage runs, so a
# coverage-side failure must not cost it twice.
# COVERAGE-ONLY. GSE39561 (THP-1, 2012) has NO detectable 3-nt periodicity: RiboCode metaplots
# reports "No obviously periodicity are detected" for all three runs and selects zero read lengths,
# so P-sites cannot be placed and there is no observed-call arm for that dataset. That is a property
# of the library, not of the processing -- its aligned inserts peak at 30 nt with 27-32 dominant,
# a normal footprint distribution, and the poly(A) trim is correct (no ligated adapter is present in
# the reads at all). Such a pack carries RNA coverage only and supports `pred_preddepth`, the
# standalone arm, which needs no Ribo-seq.
if [[ "${HS_COVERAGE_ONLY:-0}" == "1" ]]; then
  echo "  COVERAGE-ONLY pack (no Ribo target): $D has no usable P-site signal"
  "$PY" "$NEW/scripts/prepare/prepare_pack.py" \
    --coverage-only "${RNA_ARGS[@]}" \
    --universe-tx "$UNI" --fasta "$FA" --build-orf-track --kozak none \
    --group "human_$D" --species human --out "$OUT" \
    --tx2biotype "$NEW/data/tx2biotype.tsv"
  "$PY" -c "
import numpy as np
c=np.load('$OUT/coverage.npy'); o=np.load('$OUT/offsets.npy')
n=sum(1 for _ in open('$OUT/tx_order.txt'))
assert len(c)==int(o[-1]), 'axis mismatch'
print(f'  coverage-only pack: {n:,} tx, {c.sum():,} coverage over {len(c):,} nt')
"
  echo "done=$(date -Iseconds) -> $OUT (coverage-only)"
  exit 0
fi

CACHED=("$OUT"/_work/psites/*_psites.hd5)
if [[ -e "${CACHED[0]}" && "${#CACHED[@]}" -eq "${#RIBO[@]}" ]]; then
  echo "  reusing ${#CACHED[@]} cached psites hd5"
  RIBO_ARGS=(--ribo-psites "${CACHED[@]}")
else
  RIBO_ARGS=(--ribo-bam "${RIBO[@]}")
fi

"$PY" "$NEW/scripts/prepare/prepare_from_bams.py" \
  "${RIBO_ARGS[@]}" "${RNA_ARGS[@]}" \
  --pysam-python "$PYSAM" \
  --universe-tx "$UNI" --fasta "$FA" --build-orf-track --kozak none \
  --annot "$ANNOT" --tx2biotype "$NEW/data/tx2biotype.tsv" \
  --ribocode-bin "$RBIN" \
  --group "human_$D" --species human --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-12}"

"$PY" -c "
import numpy as np
t=np.load('$OUT/target_counts.npy'); c=np.load('$OUT/coverage.npy'); o=np.load('$OUT/offsets.npy')
n=sum(1 for _ in open('$OUT/tx_order.txt'))
assert len(t)==len(c)==int(o[-1]), 'axis mismatch'
assert n==$(wc -l < "$UNI"), 'pack tx != universe tx'
print(f'  pack: {n:,} tx, {t.sum():,} P-sites, {c.sum():,} coverage over {len(t):,} nt')
"
TRACK=$OUT/orf_track_v2_nokozak.npy
[[ -s "$TRACK" ]] && echo "  ORF track: $(du -h "$TRACK" | cut -f1)" || echo "  WARNING: no ORF track written"
echo "done=$(date -Iseconds) -> $OUT"
