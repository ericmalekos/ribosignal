#!/usr/bin/env bash
#SBATCH --job-name=salmon_cho
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-55%8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/salmon_chothani/%x_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/salmon_chothani/%x_%A_%a.err
set -euo pipefail
# Recover per-tissue Chothani TPM WITHOUT re-downloading. The original decoy-aware quant died with
# the biotype_probe tree; ENA and AWS are both serving at 0.3-0.6 MB/s today, making a 67 GB refetch
# a 30-50 hour proposition. But 55 of the 57 RNA runs still have STAR transcriptome BAMs on disk
# WITH multimappers retained (NH up to 17), which is what salmon's EM needs.
#
# METHOD CAVEAT, STATED LOUDLY: this is salmon ALIGNMENT-BASED mode (-a) against the transcriptome
# only. The original was DECOY-AWARE selective alignment from FASTQ. Decoys exist to stop intronic
# and repeat-derived reads being misassigned to transcripts, so these numbers are NOT drop-in
# replacements for the originals and must never be quoted as such. Fibroblast is quantified here
# too precisely so the two methods can be compared against the surviving decoy-aware
# fibroblast_salmon_mean_tpm.tsv -- that comparison is what tells us how much the difference costs.
#
#   sbatch scripts/salmon_chothani_alnmode.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
TXFA=$NEW/data/human_ribocode_annot_primary/transcripts_sequence.fa   # matches the BAM @SQ, 509,650
MAP=$NEW/data/tpm/chothani_alnmode/srr_tissue_have.tsv
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAP")
SRR=$(echo "$row" | cut -f1); TIS=$(echo "$row" | cut -f2)
BAM=$NEW/data/rnaseq_bam_mm1/${SRR}.Aligned.toTranscriptome.out.bam
OUT=$NEW/data/tpm/chothani_alnmode/${TIS}/${SRR}
[ -f "$BAM" ] || { echo "FATAL no BAM $BAM" >&2; exit 2; }
mkdir -p "$OUT"
echo "[$(date +%T)] $SRR ($TIS) -> $OUT"
"$SALMON" quant -t "$TXFA" -l A -a "$BAM" -o "$OUT" -p "${SLURM_CPUS_PER_TASK:-8}" --noErrorModel
[ -s "$OUT/quant.sf" ] || { echo "FATAL no quant.sf" >&2; exit 3; }
echo "[$(date +%T)] $SRR done -- $(( $(wc -l < "$OUT/quant.sf") - 1 )) transcripts"
