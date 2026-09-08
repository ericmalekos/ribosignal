# Fetch and filter the annotation

Download GENCODE v49 and the GRCh38 primary assembly, then cut both down to chromosome 22. The
whole tutorial stays on chr22: about 1,750 genes and 11,600 transcripts, roughly 1 GB of RAM.

```bash
mkdir -p ref && cd ref
curl -sSL -o gencode.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
curl -sSL -o genome.fa.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz

gunzip -c genome.fa.gz | awk '/^>/{p=($1==">chr22")} p' > chr22.fa
gunzip -c gencode.gtf.gz | awk -F'\t' '$1=="chr22"' > chr22.gtf
samtools faidx chr22.fa
```

Extract the transcript sequences. The header is trimmed to the bare versioned transcript id,
because everything downstream joins on that.

```bash
gffread -w tx.raw -g chr22.fa chr22.gtf
awk '/^>/{print $1; next}{print}' tx.raw > transcripts.fa && rm tx.raw
grep -c '^>' transcripts.fa
```

Build the STAR index. `--sjdbOverhang` is the read length minus one, and
`--genomeSAindexNbases` must be scaled down for a single chromosome or STAR builds a poor index
without complaining.

```bash
GLEN=$(awk '{s+=$2} END{print s}' chr22.fa.fai)
NB=$(python -c "import math; print(min(14, int(math.log2($GLEN)/2 - 1)))")
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir star_index \
     --genomeFastaFiles chr22.fa --sjdbGTFfile chr22.gtf \
     --sjdbOverhang 75 --genomeSAindexNbases "$NB" --outFileNamePrefix star_index/
cd ..
```
