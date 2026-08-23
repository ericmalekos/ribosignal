#!/usr/bin/env bash
#SBATCH --job-name=mm25_rc
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/rc_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/rc_%j.err
set -euo pipefail
# DIAGNOSTIC (see align_gse243134_mm25.sbatch): RiboCode on the mm25 BAMs, invoked IDENTICALLY to
# the mm1 call in scripts/ribocode_gse243134.sbatch -- same annotation, same metaplots -> RiboCode
# flow, same flags. The ONLY upstream difference is --outFilterMultimapNmax 25 vs 1.
#
# EXPECT INFLATED COUNTS. RiboCode does not drop multimappers, so P-sites are over-counted here by
# construction. That is the point: the question is not "are these counts correct" but "do the ORFs
# the model calls, and mm1 Ribo-seq misses, APPEAR when multimappers are kept".
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
SAMTOOLS=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
ANNOT=$NEW/data/mouse_ribocode_annot
BAMDIR=$NEW/data/mm25_diagnostic/gse243134_ribo_mm25
OUT=$NEW/data/mm25_diagnostic/psites_mm25; mkdir -p "$OUT"
T=${SLURM_CPUS_PER_TASK:-8}
export PYTHONNOUSERSITE=1
BAMLIST=$OUT/rpf_bams.txt; : > "$BAMLIST"; n=0
while read -r srr; do
  raw=$BAMDIR/${srr}.Aligned.toTranscriptome.out.bam
  [[ -s "$raw" ]] || { echo "missing $raw" >&2; exit 1; }
  srt=$BAMDIR/${srr}.sorted.bam
  # STAR emits TranscriptomeSAM in READ order (no SO: line); metaplots' pysam fetch needs
  # coordinate-sorted + indexed. sed -n 1p, not `head -1 | grep -q`, which SIGPIPEs under pipefail.
  if [[ ! -s "$srt" ]] || ! "$SAMTOOLS" quickcheck "$srt" 2>/dev/null \
     || [[ "$("$SAMTOOLS" view -H "$srt" | sed -n 1p)" != *SO:coordinate* ]]; then
    "$SAMTOOLS" sort -@ "$T" -m 2G -o "$srt" "$raw"
  fi
  [[ -s "$srt.bai" ]] || "$SAMTOOLS" index -@ "$T" "$srt"
  echo "$srt" >> "$BAMLIST"; n=$((n+1))
done < "$NEW/data/liver3x3/gse243134_ribo_runs.txt"
echo "[$(date +%T)] $n mm25 RPF BAMs pooled"
[[ $n -eq 21 ]] || { echo "expected 21, got $n" >&2; exit 1; }
cd "$OUT"
echo "[$(date +%T)] metaplots"
"$RC/metaplots" -a "$ANNOT" -i "$BAMLIST" -o gse243134_mm25 > "$OUT/metaplots.log" 2>&1
CONFIG=$OUT/gse243134_mm25_pre_config.txt
[[ -f "$CONFIG" ]] || { echo "no pre_config" >&2; tail -20 "$OUT/metaplots.log"; exit 2; }
echo "  libraries with a P-site offset: $(grep -c '^' "$CONFIG")"
echo "[$(date +%T)] RiboCode joint ORF call"
"$RC/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o gse243134_mm25 > "$OUT/ribocode.log" 2>&1
echo "[$(date +%T)] DONE: $(($(wc -l < "$OUT"/gse243134_mm25_collapsed.txt) - 1)) collapsed calls (mm1 gave 12,584 filtered)"
