# Get the sequencing data

One donor, both assays, from Chothani et al. (GEO GSE182371 Ribo-seq, GSE182372 RNA-seq),
Hepatocytes_1. Hepatocytes is the tissue held out of training for both released checkpoints, so
this is a genuine held-out test.

| run | assay | layout | role |
|---|---|---|---|
| `SRR15513269` | RNA-seq | paired, 76 nt | the model's input |
| `SRR15513208` | Ribo-seq | single | the ground truth |

Fetch from ENA directly, which needs no SRA toolkit.

```bash
mkdir -p fastq && cd fastq
for f in SRR15513269_1 SRR15513269_2; do
  curl -sSL -O "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/069/SRR15513269/$f.fastq.gz"
done
curl -sSL -O "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/008/SRR15513208/SRR15513208.fastq.gz"
ls -la && cd ..
```
