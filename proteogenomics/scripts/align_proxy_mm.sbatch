#!/usr/bin/env bash
#SBATCH --job-name=proxy_mm
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/proxy_mm_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/proxy_mm_%j.err
set -euo pipefail
# mm1-vs-mm20 proxy: align the SAME 30M DoHH2 read pairs with --outFilterMultimapNmax=$MULTIMAP, build coverage.
# Only variable is the multimap posture (same reads, same everything else). Pass via --export: MULTIMAP (1 or 20).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
R1=$NEW/proteogenomics/data/DoHH2_pilot/rnaseq/SRR12285182_1.fastq.gz
R2=$NEW/proteogenomics/data/DoHH2_pilot/rnaseq/SRR12285182_2.fastq.gz
P=$NEW/proteogenomics/data/mm_proxy; mkdir -p "$P"
T=${SLURM_CPUS_PER_TASK:-16}; SCR=/data/tmp/emalekos/proxy_mm${MULTIMAP}_${SLURM_JOB_ID}; mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT
echo "host=$(hostname) MULTIMAP=$MULTIMAP start=$(date -Iseconds)"
# same 30M subsample for both postures (deterministic head, R1/R2 in sync)
set +o pipefail
zcat "$R1" | head -n 120000000 | gzip > "$SCR/s1.fq.gz"
zcat "$R2" | head -n 120000000 | gzip > "$SCR/s2.fq.gz"
set -o pipefail
"$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$SCR/s1.fq.gz" "$SCR/s2.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax "$MULTIMAP" --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/DoHH2."
grep -E "Uniquely mapped reads number|Uniquely mapped reads %|multiple loci" "$SCR/DoHH2.Log.final.out" | head
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$SCR/DoHH2.Aligned.toTranscriptome.out.bam" \
  "$P/DoHH2_mm${MULTIMAP}_coverage.hd5" "DoHH2_mm${MULTIMAP}"
echo "done=$(date -Iseconds) -> $P/DoHH2_mm${MULTIMAP}_coverage.hd5"