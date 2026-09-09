# The same run as a Nextflow pipeline

Everything in this tutorial is also a Nextflow pipeline, so a whole-genome run does not need the
intermediates managed by hand. It mirrors `scripts/demo/run_demo.sh` step for step.

The quick start fetches its own reference and data, and runs chr22 end to end.

```bash
nextflow run . -profile test,docker,cpu
```

With your own FASTQ, give it the reference and the reads. The adapter is measured unless you set
`--ribo_adapter`.

```bash
nextflow run . -profile docker,gpu \
  --rna_fastq  'fastq/*_{1,2}.fastq.gz' \
  --ribo_fastq 'fastq/ribo.fastq.gz' \
  --gtf ref/annotation.gtf --fasta ref/genome.fa \
  --outdir results
```

If you already have STAR transcriptome BAMs, pass them instead and the pipeline skips both the
alignment and the STAR index.

```bash
nextflow run . -profile docker,cpu \
  --rna_bam rna.toTranscriptome.bam --ribo_bam ribo.toTranscriptome.bam \
  --gtf ref/annotation.gtf --fasta ref/genome.fa
```

On SLURM, add the `slurm` profile. Nextflow submits one job per process.

```bash
nextflow run . -profile slurm,singularity,gpu --outdir results -resume
```

## Profiles

| profile | effect |
|---|---|
| `docker` / `singularity` | run every process in the container |
| `local` | no container: use tools already on `PATH`, plus this checkout's `scripts/` |
| `cpu` | `--device cpu`, and the CPU image, which has no `mamba_ssm` |
| `gpu` | `--device cuda`, the GPU image, one accelerator per predict task |
| `slurm` | submit each process to the scheduler |
| `test` | chr22 with the two public accessions used in this tutorial |

## Parameters worth setting

| parameter | default | note |
|---|---|---|
| `--arch` | `attn,mamba4` | on GPU mamba4 is the faster of the two; on CPU it is 6.3x slower |
| `--ribo_adapter` | `null` | `null` measures it, `''` asserts already-trimmed, a string forces it |
| `--chrom` | `chr22` | `all` for the primary assembly, which needs about 32 GB to index |
| `--max_tx_length` | `10000` | the released checkpoints never saw a longer transcript |
| `--pred_scale` | `0.05` | the Poisson dial on the second calling arm |
| `--star_index` | `null` | reuse an existing index instead of building one |

`nextflow run . --help` prints the same list.

## What it writes

```
results/
  reference/     GTF, genome FASTA, transcript FASTA, STAR index
  align/         transcriptome BAMs and STAR logs
  pack/          coverage.hd5, psites.hd5
  pack/          tx_order.txt, offsets, lengths, target_counts, coverage, ORF track
  pred/          <arch>/pred_profiles.npz and the per-arch timing
  calls/         RiboCode real and pred_preddepth arms, the Poisson arm, and Ribo-TISH
  RESULTS.txt    predicted calls scored against the observed Ribo-seq
```

Both calling arms are always written. Quoting the deterministic arm without the Poisson one, or
the reverse, gives a number better than the method is.
