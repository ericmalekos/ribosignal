#!/usr/bin/env bash
# Align one Fibroblast RNAseq PE sample and write per-nt transcriptome coverage via
# deepTools bamCoverage (multithreaded). Faster + full core use vs the pysam counter.
#   STAR (TranscriptomeSAM, no genomic BAM) -> samtools sort+index -> bamCoverage
#   (--binSize 1 --filterRNAstrand forward [ISR sense], secondary counted [posture A],
#   supplementary excluded) -> pyBigWig -> per-sample coverage hd5. All in node-local
#   scratch; the BAMs never touch the group filesystem. Resumable (skip if hd5 exists).
set -euo pipefail
SRR="$1"

RIBOSEQ=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
DEEPT=/private/groups/carpenterlab/emalekos/conda_envs/perturb-atac-py
CAS12A=/private/groups/carpenterlab/emalekos/conda_envs/cas12a
IDX=/private/groups/carpenterlab/emalekos/genomes/star_index_grch38_v49
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
FQ1=$ECH/data/rnaseq_input/${SRR}_1.fastq.gz
FQ2=$ECH/data/rnaseq_input/${SRR}_2.fastq.gz
COVDIR=$NEW/data/rnaseq_coverage/per_sample
OUTCOV=${2:-$COVDIR/${SRR}_coverage.hd5}   # optional 2nd arg overrides output path (validation smoke)
LOGDIR=$NEW/logs/rnaseq_star
NT=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$COVDIR" "$LOGDIR"

if [ -s "$OUTCOV" ]; then echo "[$SRR] coverage exists, skipping"; exit 0; fi
for f in "$FQ1" "$FQ2"; do [ -s "$f" ] || { echo "MISSING $f" >&2; exit 1; }; done

SCRATCH=/data/tmp/emalekos/rnaseq_star/$SRR
rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"
echo "[$SRR] host=$(hostname) nt=$NT start=$(date -Iseconds)"

# 1) STAR -> transcriptome BAM (no genomic BAM)
"$RIBOSEQ/bin/STAR" --genomeDir "$IDX" \
  --readFilesIn "$FQ1" "$FQ2" --readFilesCommand zcat --runThreadN "$NT" \
  --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCRATCH/${SRR}."
TX="$SCRATCH/${SRR}.Aligned.toTranscriptome.out.bam"
[ -s "$TX" ] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCRATCH/${SRR}.Log.final.out" "$LOGDIR/${SRR}.Log.final.out" 2>/dev/null || true

# 2) coordinate sort + index (bamCoverage needs an indexed BAM)
"$RIBOSEQ/bin/samtools" sort -@ "$NT" -m 2G -T "$SCRATCH/st" -o "$SCRATCH/${SRR}.sorted.bam" "$TX"
rm -f "$TX"
"$RIBOSEQ/bin/samtools" index -@ "$NT" "$SCRATCH/${SRR}.sorted.bam"

# 3) per-nt coverage: ISR sense (forward), count secondary (posture A), drop supplementary (0x800)
"$DEEPT/bin/bamCoverage" -b "$SCRATCH/${SRR}.sorted.bam" -o "$SCRATCH/${SRR}.cov.bw" \
  --binSize 1 --filterRNAstrand forward --normalizeUsing None \
  --samFlagExclude 2048 --numberOfProcessors "$NT" --outFileFormat bigwig

# 4) bigWig -> hd5
"$CAS12A/bin/python3" "$NEW/scripts/bw_to_coverage_hd5.py" \
  "$SCRATCH/${SRR}.cov.bw" "$OUTCOV" "$SRR"

rm -rf "$SCRATCH"
echo "[$SRR] done=$(date -Iseconds)  -> $OUTCOV"
