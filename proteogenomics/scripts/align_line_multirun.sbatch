#!/usr/bin/env bash
#SBATCH --job-name=align_line_mr
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_mr_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/align_mr_%j.err
set -euo pipefail
# RNA-seq -> salmon TPM (expressed universe) + STAR per-nt coverage (the model input), for a line with
# MULTIPLE runs, at FULL DEPTH, at a chosen multimap posture.
#
# Supersedes align_line.sbatch for the immunopeptidome lines. Two things it fixes:
#   1. NO SUBSAMPLING. align_line.sbatch defaulted to the first 30M read pairs (see its header for why
#      that is wrong). Full depth only.
#   2. ALL RUNS, not one. HBL1/DoHH2/SU-DHL-4 each have FOUR RNA runs in PRJNA647736; the original
#      panel used one per line. Runs are comma-joined into a single STAR invocation and passed
#      together to salmon as one sample, exactly as process_rnaseq_thp1.sbatch already did for THP1.
#
# Env: LINE, RUNS (space-separated accessions), MM (1|10), SPECIES (default human),
#      FQDIR (default <pilot>/rnaseq), SKIP_SALMON (default 0 -- salmon MUST re-run when the run set
#      changes, because the expressed universe is derived from it).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
LINE=${LINE:?set LINE}; MM=${MM:?set MM=1 or 10}
read -r -a RUNARR <<< "${RUNS:?set RUNS='SRR1 SRR2 ...'}"
SPECIES=${SPECIES:-human}
if [ "$SPECIES" = "mouse" ]; then
  STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
  SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
else
  STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
  SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
fi
P=$NEW/proteogenomics/data/${LINE}_pilot
FQDIR=${FQDIR:-$P/rnaseq}
T=${SLURM_CPUS_PER_TASK:-12}
COVDIR=$P/coverage_mm${MM}_full
SALDIR=$P/salmon_full/${LINE}
# SCRATCH ON THE GROUP FS, NOT /data/tmp. These are full-depth multi-run transcriptome BAMs
# (tens of GB each, larger at mm10), and SLURM packed 5 of 8 jobs onto one node -- node-local
# /data/tmp filled and samtools died with "bgzf_flush: File write failed (wrong size)" AFTER STAR
# had finished mapping. The group fs has TB free. Node-local scratch is the right default for small
# intermediates; it is the wrong default when one node may host several multi-GB writers.
SCR=$NEW/data/_scratch_align/${LINE}_mm${MM}_${SLURM_JOB_ID:-m}
mkdir -p "$SCR" "$COVDIR" "$(dirname "$SALDIR")"
echo "host=$(hostname) LINE=$LINE MM=$MM runs=${RUNARR[*]} start=$(date -Iseconds)"
R1=(); R2=()
for s in "${RUNARR[@]}"; do
  a=$FQDIR/${s}_1.fastq.gz; b=$FQDIR/${s}_2.fastq.gz
  [[ -s "$a" && -s "$b" ]] || { echo "missing fastq for $s under $FQDIR" >&2; exit 1; }
  R1+=("$a"); R2+=("$b")
done
echo "  ${#R1[@]} run pairs, full depth, no subsampling"
# salmon: all runs as ONE sample. Only skipped if an identical-run-set quant already exists.
if [ "${SKIP_SALMON:-0}" != "1" ] && [ ! -s "$SALDIR/quant.sf" ]; then
  "$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "${R1[@]}" -2 "${R2[@]}" \
    --validateMappings --gcBias --seqBias -p "$T" -o "$SALDIR"
  echo "  salmon done; tx>1TPM $(awk 'NR>1 && $4>1' "$SALDIR/quant.sf" | wc -l)"
else
  echo "  salmon skipped (quant.sf present or SKIP_SALMON=1)"
fi
IFS=,; RR1="${R1[*]}"; RR2="${R2[*]}"; unset IFS
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$RR1" "$RR2" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax "$MM" --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/${LINE}."
cp "$SCR/${LINE}.Log.final.out" "$COVDIR/${LINE}.Log.final.out"
grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci" \
  "$COVDIR/${LINE}.Log.final.out" | sed 's/^/  /'
TBAM=$SCR/${LINE}.Aligned.toTranscriptome.out.bam
"$ENV/bin/samtools" quickcheck "$TBAM" || { echo "TRUNCATED $TBAM" >&2; exit 2; }
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$COVDIR/${LINE}_coverage.hd5" "$LINE"
BAMDIR=$P/bam_mm${MM}_full; mkdir -p "$BAMDIR"
mv "$TBAM" "$BAMDIR/" && rm -rf "$SCR"
echo "done=$(date -Iseconds) -> $COVDIR/${LINE}_coverage.hd5"
