#!/usr/bin/env bash
#SBATCH --job-name=align_a549enc
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_a549enc_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_a549enc_%j.err
set -euo pipefail
# UC-A: A549 model inputs from the ENCODE raw FASTQ (ENCSR000CON, PE 2x76, stranded-reverse, polyA+).
# Both bio-reps combined for depth. GENUINE raw FASTQ (no BAM revert). Faithful training recipe:
#   STAR --quantMode TranscriptomeSAM -> rnaseq_coverage.py -> per-nt coverage (THE model input, A549-specific)
#   salmon decoy-aware v49 -> isoform TPM (expressed filter + 3-frame null)
# rnaseq_coverage.py auto-detects strand (handles stranded-reverse); salmon -l A auto-detects (-> ISR).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RNA=$NEW/proteogenomics/data/A549_pilot/encode_rna
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
OUT=$NEW/proteogenomics/data/A549_pilot
T=${SLURM_CPUS_PER_TASK:-8}
SCR=/data/tmp/emalekos/a549enc_${SLURM_JOB_ID:-manual}
mkdir -p "$SCR" "$OUT/salmon" "$OUT/coverage"
R1=("$RNA/ENCFF000EJJ.fastq.gz" "$RNA/ENCFF000EJW.fastq.gz")   # rep1 R1, rep2 R1
R2=("$RNA/ENCFF000EJV.fastq.gz" "$RNA/ENCFF000EKB.fastq.gz")   # rep1 R2, rep2 R2
for f in "${R1[@]}" "${R2[@]}"; do [ -s "$f" ] || { echo "MISSING $f" >&2; exit 1; }; done
R1CSV=$(IFS=,; echo "${R1[*]}"); R2CSV=$(IFS=,; echo "${R2[*]}")
echo "host=$(hostname) start=$(date -Iseconds)"

# salmon isoform TPM (skip if already done -- it completed in the timed-out run)
if [ -s "$OUT/salmon/A549_encode/quant.sf" ]; then
  echo "salmon quant.sf exists -- skipping salmon (tx >1 TPM: $(awk 'NR>1 && $4>1' "$OUT/salmon/A549_encode/quant.sf" | wc -l))"
else
  "$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "${R1[@]}" -2 "${R2[@]}" \
    --gcBias --seqBias -p "$T" -o "$OUT/salmon/A549_encode"
  echo "salmon done $(date -Iseconds); tx >1 TPM: $(awk 'NR>1 && $4>1' "$OUT/salmon/A549_encode/quant.sf" | wc -l)"
fi

# STAR -> transcriptome BAM -> per-nt coverage (comma-separated multi-file)
"$ENV/bin/STAR" --genomeDir "$STARIDX" \
  --readFilesIn "$R1CSV" "$R2CSV" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/A549enc."
TBAM="$SCR/A549enc.Aligned.toTranscriptome.out.bam"
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCR/A549enc.Log.final.out" "$OUT/coverage/A549_encode.Log.final.out" 2>/dev/null || true
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$OUT/coverage/A549_encode_coverage.hd5" A549_encode
echo "coverage done -> $OUT/coverage/A549_encode_coverage.hd5"
rm -rf "$SCR"
echo "all done $(date -Iseconds)"
