# Align and build coverage

The two assays use opposite settings, and the differences matter more than they look.

RNA-seq keeps multimappers up to 10 and must **not** get `--maximum-length` or
`--discard-untrimmed`, both of which are correct for footprints and wrong here.

```bash
mkdir -p rna logs
cutadapt --trim-n -m 20 -j 16 -o rna/t1.fq.gz -p rna/t2.fq.gz \
  fastq/SRR15513269_1.fastq.gz fastq/SRR15513269_2.fastq.gz > logs/cutadapt_rna.log 2>&1

STAR --genomeDir ref/star_index --readFilesIn rna/t1.fq.gz rna/t2.fq.gz --readFilesCommand zcat \
     --runThreadN 16 --outSAMtype None --quantMode TranscriptomeSAM \
     --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM \
     --outFileNamePrefix rna/SRR15513269.

grep -E "Number of input reads|Uniquely mapped reads %" rna/SRR15513269.Log.final.out
samtools quickcheck rna/SRR15513269.Aligned.toTranscriptome.out.bam
python $RIBO_SCRIPTS/rnaseq_coverage.py --bam rna/SRR15513269.Aligned.toTranscriptome.out.bam \
       --out pack/coverage.hd5 --sample SRR15513269 --min-covered-tx 1000
```

Ribo-seq is length-selected, keeps unique alignments only, and uses `EndToEnd` because
soft-clipping shifts the inferred P-site.

```bash
mkdir -p ribo pack
cutadapt -a "$RIBO_ADAPTER" -m 20 -M 40 --discard-untrimmed -j 16 \
  -o ribo/t.fq.gz fastq/SRR15513208.fastq.gz > logs/cutadapt_ribo.log 2>&1

STAR --genomeDir ref/star_index --readFilesIn ribo/t.fq.gz --readFilesCommand zcat \
     --runThreadN 16 --outSAMtype None --quantMode TranscriptomeSAM \
     --outFilterMultimapNmax 1 --alignEndsType EndToEnd \
     --outSAMattributes NH HI AS nM --outFileNamePrefix ribo/SRR15513208.

python $RIBO_SCRIPTS/ribo_psites.py --bam ribo/SRR15513208.Aligned.toTranscriptome.out.bam \
       --gtf ref/chr22.gtf --out pack/psites.hd5 --sample SRR15513208 --lengths 25:35 \
       2>&1 | tee logs/psites.log
```

The frame-0 fraction printed by the last command is the whole QC: 1/3 is noise, a real library
reaches 0.5 to 0.7.
