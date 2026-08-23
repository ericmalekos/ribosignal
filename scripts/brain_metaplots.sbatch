#!/usr/bin/env bash
#SBATCH --job-name=brain_meta
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/brain_meta_%j.out
set -euo pipefail
# QC metaplots for Brain -- the 9th Chothani tissue, and the only one with no QC artifact.
#
# WHY IT WAS MISSING. Brain was excluded from the source study's usable panel for low periodicity
# (period_obs 0.044) and was never fetched, aligned, or plotted here. So the exclusion rested on a
# single recorded number with nothing to point at. Its 5 runs (SRR15513148-152) were downloaded and
# aligned with the final recipe 2026-08-15; this produces the missing plots.
#
# OUTPUT PATH IS DELIBERATE. results/chothani_regeneration/_canonical_Brain/ is exactly the pattern
# `figures/S_riboseq_qc` globs (`_canonical_<T>/canon_<T>_pre_config.txt`), so the QC panel picks
# Brain up with no code change and becomes complete at 9 tissues.
#
# A ZERO-READ-LENGTH RESULT IS AN ANSWER, NOT A FAILURE. Brain's unique mapping is 29.55%, the same
# range as GSE39561 (29.5-32.7%), whose final-recipe rebuild selected ZERO periodic read lengths. If
# Brain does the same it is a second verified negative control; if it selects read lengths but with
# weak periodicity, that is a different and milder failure mode. Do not treat an empty config as a
# broken run -- read the PDFs.
#
#   sbatch scripts/brain_metaplots.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/human_ribocode_annot_primary
BAMDIR=$NEW/data/chothani_bam/Brain
OUT=$NEW/results/chothani_regeneration/_canonical_Brain
mkdir -p "$OUT"

mapfile -t BAMS < <(ls "$BAMDIR"/*.Aligned.toTranscriptome.out.bam | sort)
[[ ${#BAMS[@]} -eq 5 ]] || { echo "ABORT: ${#BAMS[@]} BAMs, expected 5" >&2; exit 1; }
printf '%s\n' "${BAMS[@]}" > "$OUT/bam_list.txt"
echo "Brain metaplots: ${#BAMS[@]} libraries host=$(hostname) start=$(date -Iseconds)"
printf '  %s\n' "${BAMS[@]##*/}"

( cd "$OUT" && "$RBIN/metaplots" -a "$ANNOT" -i "$OUT/bam_list.txt" -o canon_Brain )

C=$OUT/canon_Brain_pre_config.txt
if [[ -s "$C" ]]; then
  # Count REAL sample rows: the config interleaves per-read-length comment tables and blank lines,
  # so `grep -vc '^#'` counts the blanks too.
  n=$(awk -F'\t' 'NF>=5 && $1 !~ /^#/ && $1 != "" && $1 != "SampleName"' "$C" | wc -l)
  echo "  metaplots selected offsets for $n of ${#BAMS[@]} libraries"
  if [[ "$n" -eq 0 ]]; then
    echo "  RESULT: ZERO periodic read lengths -- Brain is a NEGATIVE CONTROL like GSE39561,"
    echo "          and the source study's exclusion is confirmed on our own pipeline."
  else
    echo "  RESULT: $n libraries have a usable P-site offset; Brain is weak, not unusable."
  fi
else
  echo "  no pre_config written"
fi
echo "  PDFs: $(ls "$OUT"/*.pdf 2>/dev/null | wc -l)"
echo "done=$(date -Iseconds) -> $OUT"
