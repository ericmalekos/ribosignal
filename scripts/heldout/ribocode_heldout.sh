#!/usr/bin/env bash
# Pooled RiboCode P-site build for one held-out dataset: metaplots (aggregate per-read-length
# P-site offsets over the dataset's cleaned transcriptome BAMs) -> RiboCode (joint ORF call +
# per-nt *_psites.hd5). Mirrors expression_context_human phase30 run_ribocode_per_tissue.sh.
# Usage: ribocode_heldout.sh <DATASET>
set -euo pipefail
DATASET="${1:?DATASET required}"
HELDOUT_DIR=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/heldout
export DATASET
source "$HELDOUT_DIR/heldout_config.sh"
mkdir -p "$PSITE_DIR" "$LOG_DIR"

BAM_LIST="$PSITE_DIR/${DATASET}_bam_list.txt"; : > "$BAM_LIST"
for srr in $RIBO_SRRS; do
  bam="$BAM_DIR/${srr}.Aligned.toTranscriptome.out.bam"
  [[ -f "$bam" ]] || { echo "MISSING bam: $bam" >&2; exit 2; }
  [[ -f "$bam.bai" ]] || "$RIBOSEQ_ENV/bin/samtools" index "$bam"
  echo "$bam" >> "$BAM_LIST"
done
echo "[$(date +%T)] $DATASET: $(wc -l < "$BAM_LIST") pooled BAMs"

cd "$PSITE_DIR"
echo "[$(date +%T)] $DATASET: metaplots (P-site offsets)"
"$RIBOCODE_BIN/metaplots" -a "$ANNOT" -i "$BAM_LIST" -o "$DATASET" \
  > "$LOG_DIR/${DATASET}.metaplots.log" 2>&1
CONFIG="$PSITE_DIR/${DATASET}_pre_config.txt"
[[ -f "$CONFIG" ]] || { echo "metaplots produced no config: $CONFIG" >&2; exit 4; }
echo "[$(date +%T)] $DATASET: pre_config has $(wc -l < "$CONFIG") rows"

echo "[$(date +%T)] $DATASET: RiboCode joint ORF call + per-nt psites"
"$RIBOCODE_BIN/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o "$DATASET" \
  > "$LOG_DIR/${DATASET}.ribocode.log" 2>&1

echo "[$(date +%T)] $DATASET: DONE"
ls -lh "$PSITE_DIR"/*_psites.hd5 2>/dev/null | sed 's/^/  /'
n=$(($(wc -l < "$PSITE_DIR/${DATASET}_collapsed.txt" 2>/dev/null || echo 1) - 1))
echo "  RiboCode collapsed ORF calls: $n"
