#!/usr/bin/env bash
#SBATCH --job-name=align_line
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_line_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_line_%j.err
set -euo pipefail
# Generalized RNA-seq -> salmon TPM (expressed universe) + STAR coverage (model input), mirroring align_hbl1.
# Pass via --export: LINE (dir name), R1, R2 (paired fastqs), SPECIES (human|mouse, default human).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
SPECIES=${SPECIES:-human}
if [ "$SPECIES" = "mouse" ]; then
  STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
  SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
else
  STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
  SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
fi
P=$NEW/proteogenomics/data/${LINE}_pilot
T=${SLURM_CPUS_PER_TASK:-8}; SCR=/data/tmp/emalekos/${LINE}_${SLURM_JOB_ID:-m}; mkdir -p "$SCR" "$P/salmon" "$P/coverage"
echo "host=$(hostname) LINE=$LINE SPECIES=$SPECIES start=$(date -Iseconds)"
[[ -f "$R1" && -f "$R2" ]] || { echo "missing fastqs $R1 $R2" >&2; exit 1; }
# Subsample RNA depth: the model only needs per-nt coverage + salmon TPM, both saturate well below full depth.
# Deep libraries (89-130M frags) make STAR + rnaseq_coverage.py take 4h; ~30M is plenty and ~4x faster.
SUBSAMPLE=${SUBSAMPLE:-30000000}   # read pairs to keep; 0 = full depth
if [ "$SUBSAMPLE" -gt 0 ]; then
  echo "[$(date +%T)] subsampling to first ${SUBSAMPLE} read pairs (head, R1/R2 in sync)"
  set +o pipefail
  zcat "$R1" | head -n $((SUBSAMPLE * 4)) | gzip > "$SCR/sub_1.fq.gz"
  zcat "$R2" | head -n $((SUBSAMPLE * 4)) | gzip > "$SCR/sub_2.fq.gz"
  set -o pipefail
  R1="$SCR/sub_1.fq.gz"; R2="$SCR/sub_2.fq.gz"
  echo "[$(date +%T)] subsample done: $(( $(zcat "$R1" | wc -l) / 4 )) pairs"
fi
# SKIP_SALMON=1 reuses the existing expressed universe (salmon is unaffected by the STAR multimap posture,
# so an mm1 re-run of coverage does not need to re-quant) -- saves ~10 min per line.
if [ "${SKIP_SALMON:-0}" != "1" ]; then
  "$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "$R1" -2 "$R2" --gcBias --seqBias -p "$T" -o "$P/salmon/${LINE}"
  echo "salmon done: tx>1TPM $(awk 'NR>1 && $4>1' "$P/salmon/${LINE}/quant.sf" | wc -l)"
else
  echo "SKIP_SALMON=1: reusing existing $P/salmon/${LINE} + universe (unchanged by mm1)"
fi
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${LINE}."
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$SCR/${LINE}.Aligned.toTranscriptome.out.bam" "$P/coverage/${LINE}_coverage.hd5" "${LINE}"
rm -rf "$SCR"; echo "done=$(date -Iseconds)"