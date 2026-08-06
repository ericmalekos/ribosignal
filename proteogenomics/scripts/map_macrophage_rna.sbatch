#!/usr/bin/env bash
#SBATCH --job-name=macro_rna_map
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-35%6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_rna_map_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_rna_map_%A_%a.err
set -euo pipefail
# Map PRJNA482293 macrophage RNA-seq (SINGLE-END 1x150, mouse GRCm39/vM38) for the proteogenomics search DB.
# Per sample: salmon decoy-aware -> isoform TPM (expressed universe); STAR unique-mapper -> transcriptome BAM
# -> rnaseq_coverage.py per-nt coverage (THE model input). Matches the A549/GSE proteogenomics recipe, adapted
# to single-end (salmon -r, no --gcBias; STAR one readFilesIn). Self-skips samples not fully downloaded or
# already mapped, so it can be re-run as the remaining fastqs land.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
STUDY=$NEW/proteogenomics/data/macrophage_tissue
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
T=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$STUDY/salmon" "$STUDY/coverage" "$STUDY/bam" "$NEW/proteogenomics/logs"

# pick this task's SRR + expected size from the manifest (row order = array index)
mapfile -t ROWS < <(tail -n +2 "$STUDY/rnaseq_manifest.tsv")
[ "$SLURM_ARRAY_TASK_ID" -lt "${#ROWS[@]}" ] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 0; }
IFS=$'\t' read -r SRR POP LAYOUT EXPB MD5 FTP <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
FQ="$STUDY/rnaseq/$SRR.fastq.gz"

# gate: only map a fully-downloaded fastq -- ENA size must equal the published byte count, OR the file came
# from NCBI prefetch (marker .prefetched; prefetch validates the .sra, so the fastq is trusted but its bytes
# differ from ENA's).
if [ ! -s "$FQ" ] || { [ "$(stat -c%s "$FQ")" != "$EXPB" ] && [ ! -f "$STUDY/rnaseq/$SRR.prefetched" ]; }; then
  echo "SKIP $SRR ($POP): fastq not fully downloaded ($(stat -c%s "$FQ" 2>/dev/null || echo 0)/$EXPB)"; exit 0
fi
if [ -s "$STUDY/coverage/${SRR}_coverage.hd5" ] && [ -s "$STUDY/salmon/$SRR/quant.sf" ]; then
  echo "SKIP $SRR ($POP): already mapped"; exit 0
fi
SCR=/data/tmp/emalekos/macro_map_${SRR}
mkdir -p "$SCR"
echo "MAP $SRR ($POP) host=$(hostname) start=$(date -Iseconds)"

# --- salmon: decoy-aware isoform TPM (single-end: -r, -l A auto lib type, --seqBias; --gcBias is PE-only) ---
"$ENV/bin/salmon" quant -i "$SALIDX" -l A -r "$FQ" --seqBias -p "$T" \
  -o "$STUDY/salmon/$SRR"
echo "  salmon done: tx>1TPM $(awk 'NR>1 && $4>1' "$STUDY/salmon/$SRR/quant.sf" | wc -l)"

# --- STAR: unique-mapper transcriptome BAM (model coverage input); no genome BAM ---
"$ENV/bin/STAR" --genomeDir "$STARIDX" \
  --readFilesIn "$FQ" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/$SRR."
TBAM="$SCR/$SRR.Aligned.toTranscriptome.out.bam"
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM for $SRR" >&2; exit 1; }
cp "$SCR/$SRR.Log.final.out" "$STUDY/coverage/${SRR}.Log.final.out" 2>/dev/null || true

# --- per-nt coverage hd5 (the model input); guard fails loudly on near-empty coverage ---
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" \
  "$STUDY/coverage/${SRR}_coverage.hd5" "$SRR"
mv "$TBAM" "$STUDY/bam/${SRR}.Aligned.toTranscriptome.out.bam"   # keep BAM (project rule)
rm -rf "$SCR"
echo "DONE $SRR ($POP) $(date -Iseconds)"
