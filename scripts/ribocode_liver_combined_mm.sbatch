#!/usr/bin/env bash
#SBATCH --job-name=liv_comb
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/comb_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/comb_%j.err
set -euo pipefail
# DIAGNOSTIC (see data/mm25_diagnostic/README.md). COMBINED mouse-liver RiboCode over all three
# datasets pooled -- gse243134 (21) + janich (5) + wang (2) = 28 RPF libraries.
#
#   MM=25  the maximum-sensitivity reference: every liver library, multimappers kept.
#   MM=1   the matched control at the SAME recipe, so the mm1-vs-mm25 delta is not confounded by
#          the differing historical recipes of the three datasets.
#
# This is the reference the model-specific set was actually defined against (the UNION of all three
# real liver datasets), so it is the right denominator -- gse243134 alone tests only one third of
# that definition.
#
#   MM=25 sbatch scripts/ribocode_liver_combined_mm.sbatch
#   MM=1  sbatch scripts/ribocode_liver_combined_mm.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
SAMTOOLS=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
ANNOT=$NEW/data/mouse_ribocode_annot
MM=${MM:-25}
D=$NEW/data/mm25_diagnostic
OUT=$D/psites_combined_mm${MM}; mkdir -p "$OUT"
T=${SLURM_CPUS_PER_TASK:-8}
export PYTHONNOUSERSITE=1
BAMLIST=$OUT/rpf_bams.txt; : > "$BAMLIST"; n=0
add() {  # $1 = dir of raw STAR transcriptome BAMs
  [ -d "$1" ] || { echo "missing dir $1" >&2; exit 1; }
  for raw in "$1"/*.Aligned.toTranscriptome.out.bam; do
    [ -s "$raw" ] || continue
    srt=${raw%.Aligned.toTranscriptome.out.bam}.sorted.bam
    # STAR TranscriptomeSAM is READ-ordered with no @HD SO: line; metaplots' pysam fetch needs
    # coordinate-sorted+indexed. sed -n 1p, not `head -1|grep -q`, which SIGPIPEs under pipefail.
    # quickcheck, NOT just the header: a BAM truncated mid-write still has a valid header and
    # would be silently accepted. That is exactly how job 36834132 died -- a concurrent sorter
    # truncated a file this guard had already waved through. quickcheck verifies the EOF block.
    if [ ! -s "$srt" ] || ! "$SAMTOOLS" quickcheck "$srt" 2>/dev/null \
       || [ "$("$SAMTOOLS" view -H "$srt" | sed -n 1p)" != *SO:coordinate* ]; then
      "$SAMTOOLS" sort -@ "$T" -m 2G -o "$srt" "$raw"
    fi
    [ -s "$srt.bai" ] || "$SAMTOOLS" index -@ "$T" "$srt"
    echo "$srt" >> "$BAMLIST"; n=$((n+1))
  done
}
if [ "$MM" = 25 ]; then add "$D/gse243134_ribo_mm25"; fi
add "$D/janich_ribo_mm${MM}"
add "$D/wang_ribo_mm${MM}"
if [ "$MM" = 1 ]; then
  # gse243134 mm1 already exists in the archive at the recipe its own script defined
  while read -r s; do
    b=$NEW/data/heldout_bam/mouse_gse243134_liver/${s}.Aligned.toTranscriptome.out.bam
    [ -s "$b" ] && { echo "$b" >> "$BAMLIST"; n=$((n+1)); }
  done < $NEW/data/liver3x3/gse243134_ribo_runs.txt
fi
echo "[$(date +%T)] MM=$MM  pooled $n RPF BAMs across 3 liver datasets"
[ "$n" -eq 28 ] || { echo "expected 28, got $n" >&2; exit 1; }
cd "$OUT"
"$RC/metaplots" -a "$ANNOT" -i "$BAMLIST" -o liver_combined_mm${MM} > "$OUT/metaplots.log" 2>&1
CFG=$OUT/liver_combined_mm${MM}_pre_config.txt
[ -f "$CFG" ] || { echo "no pre_config" >&2; tail -20 "$OUT/metaplots.log"; exit 2; }
echo "  P-site offsets: $(grep -c '^' "$CFG")"
"$RC/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o liver_combined_mm${MM} > "$OUT/ribocode.log" 2>&1
echo "[$(date +%T)] DONE MM=$MM: $(($(wc -l < "$OUT"/liver_combined_mm${MM}_collapsed.txt) - 1)) collapsed calls"
