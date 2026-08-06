#!/usr/bin/env bash
#SBATCH --job-name=tcell_rna_map
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/tcell_rna_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/tcell_rna_%A_%a.err
set -euo pipefail
# Map GSE155087 WT CD4+ T-cell mRNA-seq (3 reps, SINGLE-END, mouse GRCm39/vM38). Per sample: salmon decoy-aware
# -> isoform TPM (for the T-cell expressed universe); STAR unique-mapper -> transcriptome BAM -> rnaseq_coverage.py
# per-nt coverage (the model RNA input, mm1). RNA is adapter-clean (0.2%) so it aligns raw (macrophage recipe).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
T=${SLURM_CPUS_PER_TASK:-8}
FQDIR=$NEW/data/heldout_raw/gse155087_tcell/fastq
COVDIR=$NEW/data/heldout_rnaseq_coverage/mouse_gse155087_tcell
SALDIR=$COVDIR/salmon
BAMDIR=$NEW/data/heldout_bam/mouse_gse155087_tcell
mkdir -p "$COVDIR" "$SALDIR" "$BAMDIR" "$NEW/logs/heldout"

SRRS=(SRR12318320 SRR12318321 SRR12318322)      # WT mRNA Ctrl 1-3
SRR=${SRRS[$SLURM_ARRAY_TASK_ID]:?array index out of range}
FQ="$FQDIR/$SRR.fastq.gz"
[ -s "$FQ" ] || { echo "missing $FQ" >&2; exit 1; }
if [ -s "$COVDIR/${SRR}_coverage.hd5" ] && [ -s "$SALDIR/$SRR/quant.sf" ]; then
  echo "SKIP $SRR: already mapped"; exit 0
fi
SCR=/data/tmp/emalekos/tcell_rna_${SRR}; rm -rf "$SCR"; mkdir -p "$SCR"
echo "MAP $SRR host=$(hostname) start=$(date -Iseconds)"

# --- salmon: decoy-aware isoform TPM (single-end: -r, -l A, --seqBias; no --gcBias which is PE-only) ---
"$ENV/bin/salmon" quant -i "$SALIDX" -l A -r "$FQ" --seqBias -p "$T" -o "$SALDIR/$SRR"
echo "  salmon done: tx>1TPM $(awk 'NR>1 && $4>1' "$SALDIR/$SRR/quant.sf" | wc -l)"

# --- STAR: unique-mapper transcriptome BAM (model coverage input); no genome BAM ---
"$ENV/bin/STAR" --genomeDir "$STARIDX" \
  --readFilesIn "$FQ" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/$SRR."
TBAM="$SCR/$SRR.Aligned.toTranscriptome.out.bam"
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM for $SRR" >&2; exit 1; }
cp "$SCR/$SRR.Log.final.out" "$COVDIR/${SRR}.Log.final.out" 2>/dev/null || true

# --- per-nt coverage hd5 (the model input) ---
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$COVDIR/${SRR}_coverage.hd5" "$SRR"
mv "$TBAM" "$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam"   # keep BAM (project rule)
rm -rf "$SCR"
echo "DONE $SRR $(date -Iseconds)"
