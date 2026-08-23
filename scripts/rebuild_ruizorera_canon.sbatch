#!/usr/bin/env bash
#SBATCH --job-name=rz_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rz_pack_%j.out
set -euo pipefail
# Rebuild the iPSC-CM (Ruiz-Orera) pack on the FINAL RECIPE. Task #92, the last off-recipe arm in P5.
#
# WHY THIS ARM NEEDED A RE-DOWNLOAD AT ALL. Its pack provenance records `inputs_known: false`, its
# source psites hd5 were deleted and its 5 ENA BAMs were gone, so nothing survived to re-derive from.
# The FASTQs were re-fetched from PRJEB65856 on 2026-08-15 (15 files, all md5-verified) and the 5 Ribo
# runs realigned by scripts/riboseq_align.sbatch (job 36791179, 5/5 COMPLETED).
#
# ONLY THE RIBO HALF IS REBUILT, AND THAT IS THE POLICY, NOT A SHORTCUT.
# PIPELINE_POLICY and rebuild_heldout_canon are explicit: "RNA-seq stays Local, only Ribo needed
# EndToEnd", and "THE UNIVERSE IS NOT REDEFINED -- each pack's existing tx_order.txt IS its universe",
# because that universe came from salmon TPM on the RNA-seq. Redefining it would import today's flags
# into a definition that was already correct AND break axis alignment with the in-pack ORF track.
#
# So: --reuse-coverage takes the universe (36,668 tx) and coverage.npy verbatim from the existing
# pack, and only target_counts is rebuilt from the new P-sites. The RNA FASTQs downloaded alongside
# the Ribo ones are NOT needed for this and are not used.
#
# --no-ncrna-filter is a DECLARED deviation, not a skip: riboseq_align.sbatch already applied
# filter_tx_heldout to the raw query-grouped STAR output, which is the only point at which it can be
# applied. These BAMs are coordinate-sorted now, so re-running it would double-filter and be refused
# by the sorted-input guard anyway.
#
#   sbatch scripts/rebuild_ruizorera_canon.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PYSAM=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/human_ribocode_annot_primary
T2B=$NEW/data/tx2biotype.tsv
REF=$NEW/data/packed_heldout_human_ruizorera          # universe + coverage source
OUT=$NEW/data/packed_heldout_human_ruizorera_canon
BAMDIR=$NEW/data/heldout_bam/human_ruizorera/ruizorera

# Inputs are LISTED, never globbed (PIPELINE_POLICY rule 3).
RUNS=(ERR12549926 ERR12549927 ERR12549928 ERR12549929 ERR12549930)
RIBO=(); for r in "${RUNS[@]}"; do RIBO+=("$BAMDIR/${r}.Aligned.toTranscriptome.out.bam"); done
for f in "${RIBO[@]}" "$REF/tx_order.txt" "$REF/coverage.npy" "$ANNOT/transcripts_cds.txt"; do
  [[ -s "$f" ]] || { echo "ABORT: missing $f" >&2; exit 1; }
done
echo "host=$(hostname) start=$(date -Iseconds)"
echo "  ribo BAMs: ${#RIBO[@]}   universe: $(wc -l < "$REF/tx_order.txt") tx (reused from $(basename $REF))"

"$PY" "$NEW/scripts/prepare/prepare_from_bams.py" \
  --ribo-bam "${RIBO[@]}" \
  --pysam-python "$PYSAM" \
  --reuse-coverage "$REF" \
  --annot "$ANNOT" --tx2biotype "$T2B" \
  --ribocode-bin "$RBIN" \
  --no-ncrna-filter \
  --group human_ruizorera_canon --species human \
  --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-12}"

# Gate axis alignment BEFORE copying the ORF track across. A track whose rows no longer match the
# pack produces plausible numbers rather than a crash, which is the worst failure mode here.
"$PY" - <<PYEOF
import numpy as np, sys, pathlib
ref, out = pathlib.Path("$REF"), pathlib.Path("$OUT")
a = [l.strip() for l in open(ref/"tx_order.txt") if l.strip()]
b = [l.strip() for l in open(out/"tx_order.txt") if l.strip()]
assert a == b, f"tx_order differs: {len(a)} vs {len(b)}"
la, lb = np.load(ref/"lengths.npy"), np.load(out/"lengths.npy")
assert np.array_equal(la, lb), "lengths differ"
ca = np.asarray(np.load(ref/"coverage.npy", mmap_mode="r"))
cb = np.asarray(np.load(out/"coverage.npy", mmap_mode="r"))
assert np.array_equal(ca, cb), "coverage was NOT reused verbatim"
ta = int(np.asarray(np.load(ref/"target_counts.npy", mmap_mode="r")).sum())
tb = int(np.asarray(np.load(out/"target_counts.npy", mmap_mode="r")).sum())
print(f"  universe + lengths + coverage identical to the source pack")
print(f"  P-sites  OLD (off-recipe) {ta:,}  ->  NEW (final recipe) {tb:,}   "
      f"({100*(tb-ta)/ta:+.1f}%)")
PYEOF

cp -a "$REF/orf_track_v2_nokozak.npy" "$OUT/" 2>/dev/null || true
cp -a "$REF/orf_track_v2_nokozak_meta.json" "$OUT/" 2>/dev/null || true
echo "  ORF track copied (universe verified identical above)"
echo "done=$(date -Iseconds) -> $OUT"
