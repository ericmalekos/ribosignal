#!/usr/bin/env bash
#SBATCH --job-name=rebuild_canon
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rebuild_canon_%A_%a.out
set -euo pipefail
# Rebuild a held-out pack on CANONICAL Ribo alignments, keeping its universe exactly.
#
# These three leukocyte packs were built before docs/PIPELINE_POLICY.md: their Ribo BAMs lacked
# --alignEndsType EndToEnd and the ncRNA + cross-gene filter. The Ribo reads have since been
# re-aligned canonically (scripts/riboseq_align.sbatch -> data/ribo_bam_e2e/<arm>/), so only the
# pooling needs redoing.
#
# THE UNIVERSE IS NOT REDEFINED. Each pack's existing tx_order.txt IS its universe (pack contract),
# and that universe came from salmon TPM on the RNA-seq, which is untouched -- RNA-seq stays Local,
# only Ribo needed EndToEnd. Re-deriving it would import today's flags into a definition that was
# already correct, and would break axis alignment with the in-pack ORF track. Each arm keeps ITS OWN
# universe; these are separate tests and must not be forced onto a shared one
# (feedback_never_reuse_universe_across_datasets).
#
# THE ORF TRACK IS REUSED, not rebuilt. It is a function of sequence and annotation, not of reads,
# so an unchanged universe means an unchanged track. The script ASSERTS tx_order and lengths are
# byte-identical to the old pack before copying it across; if they are not, the track would be
# silently misaligned and the run aborts instead.
#
#   sbatch --array=0-2 scripts/rebuild_heldout_canon.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
# rnaseq_coverage.py imports pysam, which cas12a lacks; riboseq carries pysam 0.23.3.
PYSAM=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANN=/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations
REF=$NEW/data/heldout_refs

# SPECIES is a knob, not a fork (docs/PIPELINE_POLICY.md rule 1). It selects the arm list, the
# RiboCode annotation, the tx2biotype and the GENCODE FASTAs -- never which steps run.
#   SPECIES=mouse sbatch --array=0-2 ...   (gse120762 lps/nt, gse155087 tcell)
#   SPECIES=human sbatch --array=0-2 ...   (CAR-T, GSE208041 THP-1, GSE39561 THP-1)
#
# GSE39561 IS EXPECTED TO PRODUCE A ZERO-P-SITE PACK. Its library has no detectable 3-nt periodicity,
# so metaplots selects zero read lengths. That is the point of rebuilding it: if the CANONICAL recipe
# (EndToEnd + mm1 + ncRNA filter) recovers periodicity the old off-recipe alignment destroyed, that is
# a finding about the recipe; if it stays at zero, GSE39561 is a clean negative control. Do not treat
# a zero-P-site result here as a failure -- check the metaplots output.
SPECIES=${SPECIES:-mouse}
case "$SPECIES" in
  mouse) ARMS=(mouse_gse120762_lps mouse_gse120762_nt mouse_gse155087_tcell)
         RIBODIRS=(gse120762_lps gse120762_nt gse155087_tcell)
         ANNOT=$NEW/data/mouse_ribocode_annot
         T2B=$NEW/data/tx2biotype_mouse.tsv
         SRCFA=("$ANN/gencode.vM38.pc_transcripts.fa" "$ANN/gencode.vM38.lncRNA_transcripts.fa")
         RIBOROOT=$NEW/data/ribo_bam_e2e
         RNADIR=""   # per-arm, resolved below
         RNASUF=".toTranscriptome.out.bam" ;;
  human) ARMS=(human_cart human_gse208041 human_gse39561)
         RIBODIRS=(cart gse208041 gse39561)
         ANNOT=$NEW/data/human_ribocode_annot_primary
         T2B=$NEW/data/tx2biotype.tsv
         SRCFA=("$ANN/gencode.v49.pc_transcripts.fa" "$ANN/gencode.v49.lncRNA_transcripts.fa")
         RIBOROOT=$NEW/data/ribo_bam_e2e_human
         # RNA dir is PER-ARM for human: each dataset's RNA lives under its own processed/ dir,
         # and GSE39561 has none of its own so it borrows GSE208041's (its rna_runs.txt lists the
         # same four SRRs -- verified byte-identical coverage and identical tx_order in the existing
         # packs). Resolved after A is known, below.
         RNADIRS=("$NEW/data/processed/GSE304796_cart"
                  "$NEW/data/processed/GSE208041_thp1"
                  "$NEW/data/processed/GSE208041_thp1")
         RNASUF=".toTranscriptome.bam" ;;
  *) echo "unknown SPECIES=$SPECIES (want mouse|human)" >&2; exit 1 ;;
esac
A=${ARMS[$SLURM_ARRAY_TASK_ID]}
RD=${RIBODIRS[$SLURM_ARRAY_TASK_ID]}
# Per-arm RNA dir (human); mouse keeps the one-dir-per-arm default resolved further down.
# `${#arr[@]:-0}` is not valid bash (bad substitution); and under `set -u` a bare ${#arr[@]} on an
# unset array errors too, which is the mouse branch's case since it sets no RNADIRS.
if [[ -n "${RNADIRS+x}" ]] && (( ${#RNADIRS[@]} > 0 )); then
  RNADIR=${RNADIRS[$SLURM_ARRAY_TASK_ID]}
fi

OLD=$NEW/data/packed_heldout_$A
OUT=$NEW/data/packed_heldout_${A}_canon
UNI=$REF/${A}_universe_tx.txt
FASTA=$REF/${A}_universe.fa
mkdir -p "$OUT"

[[ -s "$OLD/tx_order.txt" ]] || { echo "missing old pack $OLD" >&2; exit 1; }
# The old pack's tx_order IS the universe. Materialise it under heldout_refs so the rebuild has a
# named universe file like every other dataset, rather than reaching into a pack directory.
[[ -s "$UNI" ]] || cp "$OLD/tx_order.txt" "$UNI"

if [[ ! -s "$FASTA" ]]; then
  "$PY" "$NEW/scripts/prepare/extract_universe_fasta.py" --universe-tx "$UNI" \
    --source-fasta "${SRCFA[@]}" \
    --out "$FASTA"
fi

mapfile -t RIBO < <(ls "$RIBOROOT/$RD"/*.Aligned.toTranscriptome.out.bam | sort)
# RNA comes from an EXPLICIT run list, never a glob of the directory. data/heldout_bam/<arm>/ is a
# COMBINED dir: it holds this arm's RNA alignments AND the old off-recipe Ribo alignments, plus (for
# the gse120762 arms) genome-coordinate BAMs alongside transcriptome ones. A bare `*.bam` glob
# therefore feeds Ribo runs in as RNA coverage and double-counts samples that have both BAM types --
# an error that yields a plausible pack rather than a crash. Caught before it ran; the run lists
# were derived as (transcriptome BAMs in the dir) minus (canonical Ribo runs) and written out so the
# selection is auditable rather than re-derived each time.
#
# Those derived lists were then checked against scripts/heldout/heldout_config.sh, which is the
# authoritative record of how these datasets were originally split, and match its RNA_SRRS exactly
# (lps 41/42, nt 38/39/40, tcell 20/21/22). So the ORIGINAL packs were built on the correct RNA
# runs: their only deviation is the alignment recipe, which is what this rebuild corrects.
# RNA location and suffix are per-species knobs: the mouse arms keep one dir per arm under
# data/heldout_bam/, while CAR-T's RNA sits in data/processed/GSE304796_cart/ with a different
# suffix. Assuming one layout is how the earlier glob bug happened; state both.
[[ -n "$RNADIR" ]] || RNADIR=$NEW/data/heldout_bam/$A
mapfile -t RNA < <(while read -r s; do
  echo "$RNADIR/${s}"*"$RNASUF"
done < "$REF/${A}_rna_runs.txt")
[[ ${#RIBO[@]} -ge 1 && ${#RNA[@]} -ge 1 ]] || { echo "no BAMs for $A" >&2; exit 1; }
# Belt and braces: no RNA input may be one of this arm's Ribo runs.
for r in "${RIBO[@]}"; do
  rid=$(basename "$r" | sed 's/\..*//')
  for q in "${RNA[@]}"; do
    [[ "$(basename "$q" | sed 's/\..*//')" == "$rid" ]] && {
      echo "ABORT: Ribo run $rid also listed as RNA input" >&2; exit 1; }
  done
done
for f in "${RIBO[@]}" "${RNA[@]}"; do [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }; done
echo "arm=$A ribo=${#RIBO[@]} rna=${#RNA[@]} universe=$(wc -l < "$UNI") host=$(hostname) start=$(date -Iseconds)"
printf '  ribo: %s\n' "${RIBO[@]##*/}"
printf '  rna:  %s\n' "${RNA[@]##*/}"

# --no-ncrna-filter is a DECLARED deviation, not a skip: the canonical aligner already applied
# filter_tx_heldout to the raw query-grouped STAR output, which is the only point at which it can be
# applied. These BAMs are coordinate-sorted now, so re-running it would double-filter and be refused
# by the sorted-input guard.
"$PY" "$NEW/scripts/prepare/prepare_from_bams.py" \
  --ribo-bam "${RIBO[@]}" \
  --rna-bam "${RNA[@]}" \
  --pysam-python "$PYSAM" \
  --universe-tx "$UNI" \
  --annot "$ANNOT" --tx2biotype "$T2B" \
  --ribocode-bin "$RBIN" \
  --no-ncrna-filter \
  --group "${A}_canon" --species "$SPECIES" --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-12}"

# Gate axis alignment BEFORE copying the ORF track across. A track whose rows no longer match the
# pack is the kind of error that produces plausible numbers rather than a crash.
"$PY" - <<PYEOF
import numpy as np, sys
old, new = "$OLD", "$OUT"
o = [l.strip() for l in open(old + "/tx_order.txt")]
n = [l.strip() for l in open(new + "/tx_order.txt")]
if o != n:
    sys.exit(f"ABORT: tx_order changed ({len(o):,} -> {len(n):,}); ORF track would be misaligned")
lo, ln = np.load(old + "/lengths.npy"), np.load(new + "/lengths.npy")
if not np.array_equal(lo, ln):
    sys.exit("ABORT: lengths changed; ORF track would be misaligned")
t, c = np.load(new + "/target_counts.npy"), np.load(new + "/coverage.npy")
to = np.load(old + "/target_counts.npy")
print(f"  {len(n):,} tx  P-sites {to.sum():,} -> {t.sum():,} "
      f"({100*(t.sum()-to.sum())/to.sum():+.2f}%)  coverage {c.sum():,}")
print("  OK: tx_order and lengths identical -- ORF track stays valid")
PYEOF

# Explicit if, not `[[ ]] && cp && echo`: a short-circuiting && list under `set -e` has bitten this
# project more than once, and not every arm has both tracks (gse155087_tcell ships nokozak only).
for T in orf_track_v2_nokozak.npy orf_track_v2.npy; do
  if [[ -s "$OLD/$T" && ! -s "$OUT/$T" ]]; then
    cp "$OLD/$T" "$OUT/$T"
    echo "  copied $T"
  fi
done
echo "done=$(date -Iseconds) -> $OUT"
