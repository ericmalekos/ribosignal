#!/usr/bin/env bash
#SBATCH --job-name=merged_liver_rc
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/merged_liver_rc_%j.out
set -euo pipefail
# One RiboCode run over ALL THREE mouse-liver Ribo-seq datasets pooled: Janich (5) + GSE243134 (21)
# + Wang (2) = 28 libraries, aligned with the final recipe.
#
# WHY. Every observed reference in this project so far is one dataset deep, and the deepest single
# arm (GSE243134, 143M P-sites) still calls only ~12,584 ORFs on the shared universe. The UpSet
# (figures/P9_upset_mouse3x3) shows 1,942 ORFs called by the model and by NO single experiment. Some
# of those are false positives and some are real ORFs that no one experiment was deep enough to see,
# and a single-dataset reference cannot separate the two. A merged 28-library reference is the
# deeper ground truth that can: an ORF the model called that the MERGED reference also misses is a
# much stronger false-positive candidate than one merely absent from Wang's 2 libraries.
#
# This is a joint call, not a union of three call sets. `metaplots` derives a P-site offset per
# LIBRARY (28 config lines), then `RiboCode -g` tests periodicity on the POOLED signal, so an ORF
# with weak but consistent signal in all three datasets can clear significance here while failing in
# each dataset alone. Unioning three separately-called sets could not do that.
#
# BAM SELECTION MATCHES THE POOLS EXACTLY, including GSE243134's run list rather than a directory
# glob -- that dataset ships 21 WT RPF libraries among a larger submission, and the list is the
# authoritative record of which (liver3x3_01_pool.sbatch uses the same file). If this merge silently
# used a different set the merged reference would not be comparable to the per-dataset ones it exists
# to be compared against.
#
# NOT a replacement for the per-dataset references. The 3x3 factorial needs three INDEPENDENT
# observed arms; this is a fourth, deeper arm for over-call validation
# (feedback_orf_call_eval_standing).
#
#   sbatch scripts/merged_liver_ribocode.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/mouse_ribocode_annot
E2E=$NEW/data/ribo_bam_e2e
L=$NEW/data/liver3x3
OUT=$NEW/results/merged_liver_ribocode
mkdir -p "$OUT"

mapfile -t RIBO < <({
  ls "$E2E/janich"/*.Aligned.toTranscriptome.out.bam | sort
  while read -r s; do echo "$E2E/gse243134/${s}.Aligned.toTranscriptome.out.bam"; done < "$L/gse243134_ribo_runs.txt"
  ls "$E2E/wang"/*.Aligned.toTranscriptome.out.bam | sort
})
for f in "${RIBO[@]}"; do [[ -s "$f" ]] || { echo "missing BAM: $f" >&2; exit 1; }; done
[[ ${#RIBO[@]} -eq 28 ]] || { echo "ABORT: ${#RIBO[@]} BAMs, expected 28 (5 janich + 21 gse243134 + 2 wang)" >&2; exit 1; }
[[ -d "$ANNOT" ]] || { echo "missing annotation $ANNOT" >&2; exit 1; }

printf '%s\n' "${RIBO[@]}" > "$OUT/bam_list.txt"
echo "merged liver RiboCode: ${#RIBO[@]} libraries host=$(hostname) start=$(date -Iseconds)"
printf '  %s\n' "${RIBO[@]##*/}"

# Resumable: metaplots + RiboCode is the expensive half and writes its hd5 before anything reads
# them, so a later failure should not cost the whole run again.
shopt -s nullglob
_h=("$OUT"/*_psites.hd5)
shopt -u nullglob
if [[ ${#_h[@]} -eq 0 ]]; then
  ( cd "$OUT" && "$RBIN/metaplots" -a "$ANNOT" -i "$OUT/bam_list.txt" -o merged_liver )
  CFG=$(ls "$OUT"/merged_liver_pre_config.txt | head -1)
  # Count REAL sample rows: the config interleaves per-read-length comment tables and blank lines,
  # so `grep -vc '^#'` counts the blanks too and reported 84 for 28 libraries.
  echo "  metaplots selected offsets for $(awk -F'\t' 'NF>=5 && $1 !~ /^#/ && $1 != "" && $1 != "SampleName"' "$CFG" | wc -l) libraries"
  ( cd "$OUT" && "$RBIN/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o merged_liver )
else
  echo "  reusing ${#_h[@]} cached psites hd5"
fi

# Report what the merge bought, against the deepest single arm. Raw counts only -- the
# universe-restricted, filtered comparison belongs in the scorer, not here.
"$PY" - <<PYEOF
import glob, pathlib
out = pathlib.Path("$OUT")
coll = sorted(out.glob("merged_liver*_collapsed.txt"))
if not coll:
    raise SystemExit("ABORT: RiboCode produced no *_collapsed.txt")
n = sum(1 for _ in open(coll[0])) - 1
print(f"  merged reference: {n:,} raw ORF calls -> {coll[0].name}")
for lab, p in [("janich (5 libs)", "$NEW/results/mouse_liver_3x3_canon/mamba4_ribo-janich_rna-janich/real_collapsed.txt"),
               ("gse243134 (21)", "$NEW/results/mouse_liver_3x3_canon/mamba4_ribo-gse243134_rna-gse243134/real_collapsed.txt"),
               ("wang (2 libs)",  "$NEW/results/mouse_liver_3x3_canon/mamba4_ribo-wang_rna-wang/real_collapsed.txt")]:
    q = pathlib.Path(p)
    if q.exists():
        m = sum(1 for _ in open(q)) - 1
        print(f"    vs {lab:<16} {m:,} raw calls   ({100*(n-m)/m:+.1f}%)")
print(f"  psites hd5: {len(list(out.glob('*_psites.hd5')))}")
PYEOF
echo "done=$(date -Iseconds) -> $OUT"
