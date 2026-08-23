#!/usr/bin/env bash
#SBATCH --job-name=hep_rna
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --array=0-5%6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/heprna_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/heprna_%A_%a.err
set -euo pipefail
# Does the RNA-coverage INPUT CHANNEL carry the same multimap hole as the Ribo-seq label?
#
# GTF2I ENST00000901263.1 has a 1,617 nt stretch of its CDS where observed Ribo-seq is ~0 AND the
# pack's RNA coverage is 11x below the CDS mean (12.6 vs 141.5 per nt). The model predicts near-zero
# there. That is only "memorisation" if the model had no input signal telling it to; if the RNA
# channel is itself multimap-depleted, the model is reading a corrupted feature correctly.
#
# Realigns the 2 Hepatocytes RNA libraries at mm1 (control), mm10 and mm25, recipe otherwise
# identical to scripts/align_rnaseq_sample.sh.
#
# BAMs ARE KEPT (CLAUDE.md rule 3: do not delete intermediates until the result is verified). The
# earlier in-job discard cost 74 STAR runs when the analysis turned out to measure the wrong thing.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
FQD=$NEW/data/external/chothani_rna/fastq
OUT=$NEW/data/mm25_diagnostic/hep_rna_mm
declare -a SPEC
for mm in 1 10 25; do for s in SRR15513268 SRR15513269; do SPEC+=("$s:$mm"); done; done
IFS=: read -r SRR MM <<< "${SPEC[$SLURM_ARRAY_TASK_ID]}"
BAM=$OUT/${SRR}.mm${MM}.toTranscriptome.bam
[ -s "$BAM" ] && { echo "SKIP $SRR mm$MM"; exit 0; }
R1=$FQD/${SRR}_1.fastq.gz; R2=$FQD/${SRR}_2.fastq.gz
for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done
SCR=/data/tmp/emalekos/heprna_${SRR}_mm${MM}.$$; mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --runThreadN "${SLURM_CPUS_PER_TASK:-12}" --outSAMtype None \
  --quantMode TranscriptomeSAM --outFilterMultimapNmax "$MM" \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${SRR}."
T=$SCR/${SRR}.Aligned.toTranscriptome.out.bam
[ -s "$T" ] || { echo "no transcriptome BAM for $SRR mm$MM" >&2; exit 1; }
mv "$T" "$BAM"
cp "$SCR/${SRR}.Log.final.out" "$OUT/${SRR}.mm${MM}.Log.final.out" 2>/dev/null || true
u=$(grep "Uniquely mapped reads %" "$OUT/${SRR}.mm${MM}.Log.final.out" | grep -oE "[0-9.]+%")
m=$(grep "% of reads mapped to multiple loci" "$OUT/${SRR}.mm${MM}.Log.final.out" | grep -oE "[0-9.]+%")
echo "DONE $SRR mm$MM unique=$u multi=$m  BAM KEPT -> $(basename "$BAM")"
