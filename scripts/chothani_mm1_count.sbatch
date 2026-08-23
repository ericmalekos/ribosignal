#!/usr/bin/env bash
#SBATCH --job-name=cho_mm1c
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --array=0-73%12
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/mm1c_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm25/mm1c_%A_%a.err
set -euo pipefail
# Per-transcript ALIGNMENT-RECORD counts from the EXISTING mm1 BAMs, to pair with the mm25 totals.
#
# WHY NOT NH. In a --quantMode TranscriptomeSAM BAM, NH counts compatible TRANSCRIPTS, not genomic
# loci: one genomic alignment is expanded across every isoform of the gene, so a genomically-unique
# read routinely carries NH>1. Splitting on NH measures isoform multiplicity (~99% under posture A)
# and says nothing about what --outFilterMultimapNmax 1 removed. The sound comparison is the TOTAL
# record count per transcript, mm1 vs mm25: isoform expansion inflates both sides identically.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SAM=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
D=$NEW/data/external/chothani_ribo_refetch
OUT=$NEW/data/mm25_diagnostic/chothani_mm1_counts
mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1}' "$D/manifest.tsv")
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
TSV=$OUT/${SRR}.mm1.tsv.gz
[[ -s "$TSV" ]] && { echo "SKIP $SRR"; exit 0; }
B=""
for c in "$NEW/data/chothani_bam/${SRR}.Aligned.toTranscriptome.out.bam" \
         "$NEW"/data/ribo_bam_e2e_human/*/"${SRR}.Aligned.toTranscriptome.out.bam"; do
  [[ -s "$c" ]] && { B="$c"; break; }
done
[[ -n "$B" ]] || { echo "no mm1 BAM for $SRR" >&2; exit 1; }
"$SAM" quickcheck "$B" || { echo "$B truncated" >&2; exit 1; }
"$SAM" view -@ "${SLURM_CPUS_PER_TASK:-4}" "$B" \
  | awk -F'\t' '{c[$3]++} END{for(t in c) print t"\t"c[t]}' | gzip -c > "$TSV"
echo "DONE $SRR  src=$(basename $(dirname $B))  tx=$(zcat "$TSV" | wc -l)"
