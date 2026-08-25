#!/usr/bin/env bash
#SBATCH --job-name=align_line_mm10
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_line_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_line_%j.err
set -euo pipefail
# mm10 VARIANT. Derived from align_line.sbatch by changing ONLY --outFilterMultimapNmax 1 -> 10 and
# the coverage output dir (coverage/ -> coverage_mm10/). Subsampling, salmon flags and the PE
# invocation are copied verbatim -- SUBSAMPLE must match the original or the two coverage arms are
# not comparable.
# Run with SKIP_SALMON=1: salmon is unaffected by the STAR multimap posture, so the expressed
# universe is reused unchanged. That keeps the universe IDENTICAL between the mm1 and mm10 arms,
# which is required for the ORF-call / DB comparison to isolate the coverage change.
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
T=${SLURM_CPUS_PER_TASK:-8}; SCR=/data/tmp/emalekos/${LINE}_${SLURM_JOB_ID:-m}; mkdir -p "$SCR" "$P/salmon" "$P/coverage_mm10"
echo "host=$(hostname) LINE=$LINE SPECIES=$SPECIES start=$(date -Iseconds)"
[[ -f "$R1" && -f "$R2" ]] || { echo "missing fastqs $R1 $R2" >&2; exit 1; }
# SUBSAMPLING IS OFF AND MUST STAY OFF (user directive 2026-08-24).
# It defaulted to the first 30M read pairs to save ~3h of STAR time, justified by a claim that
# coverage and TPM "saturate well below full depth" -- a claim asserted in this comment and never
# measured anywhere in results.md or methods.md. Three problems:
#   1. `zcat | head` takes the FIRST N pairs in file order, not a random sample. SRA FASTQs are
#      typically ordered by flowcell tile, so this samples a spatial subset of the flowcell.
#   2. It was applied inconsistently -- THP1 ran at full depth (82M reads) while HBL1/SUDHL4/DoHH2
#      ran at 30M, so the Fig 2a panel already mixes depths across datasets.
#   3. It suppresses exactly what a multimap-posture comparison measures: discarding reads while
#      asking how many ALIGNMENTS PER READ are retained.
# SUBSAMPLE is kept as a knob only so old runs remain reproducible; do not set it.
SUBSAMPLE=${SUBSAMPLE:-0}   # 0 = FULL DEPTH (default). Do not set. See note above.
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
# so an mm1/mm10 re-run of coverage does not need to re-quant) -- saves ~10 min per line.
if [ "${SKIP_SALMON:-0}" != "1" ]; then
  "$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "$R1" -2 "$R2" --gcBias --seqBias -p "$T" -o "$P/salmon/${LINE}"
  echo "salmon done: tx>1TPM $(awk 'NR>1 && $4>1' "$P/salmon/${LINE}/quant.sf" | wc -l)"
else
  echo "SKIP_SALMON=1: reusing existing $P/salmon/${LINE} + universe (unchanged by the STAR multimap posture)"
fi
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${LINE}."
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$SCR/${LINE}.Aligned.toTranscriptome.out.bam" "$P/coverage_mm10/${LINE}_coverage.hd5" "${LINE}"
rm -rf "$SCR"; echo "done=$(date -Iseconds)"