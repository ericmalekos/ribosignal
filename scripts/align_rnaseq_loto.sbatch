#!/usr/bin/env bash
#SBATCH --job-name=loto_rnaseq_cov
#SBATCH --partition=medium
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --array=1-25%10
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaseq_star/loto_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaseq_star/loto_%A_%a.err
set -euo pipefail
# STAR-align the 25 non-Fibroblast Chothani RNA-seq PE samples to the transcriptome and write
# per-nt coverage (posture A isoform projection, --outFilterMultimapNmax 1 = unique genomic mappers).
# mm1 is now the standard for RNA-seq coverage too: the DoHH2 mm1-vs-mm20 proxy (results.md Task 25)
# showed coverage is identical on 85.7% of tx (median per-tx Pearson 1.0000; only paralog/repeat loci
# differ), so unique mappers give smaller intermediates + faster STAR and one consistent multimap
# posture with the Ribo-seq alignment. align_rnaseq_sample.sh is tissue-agnostic + resumable
# (skips a sample whose coverage hd5 exists). Per-tissue pooling is a separate downstream step.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
LIST=$NEW/data/loto_rnaseq_srr.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -n "$SRR" ] || { echo "no SRR at row ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
mkdir -p "$NEW/logs/rnaseq_star"
bash "$NEW/scripts/align_rnaseq_sample.sh" "$SRR"
