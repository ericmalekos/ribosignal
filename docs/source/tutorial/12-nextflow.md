# Run it as a Nextflow pipeline

The pipeline runs the prediction half of this tutorial on RNA-seq alone: trim, align, build the
pack and ORF track, predict, and call ORFs from the prediction. It reads no Ribo-seq. Scoring the
calls against real Ribo-seq, as the score page does, stays a tutorial step.

The `test` profile runs chr22 end to end. It fetches the reference itself but not the reads: put
the two RNA-seq FASTQs from the sequencing-data page in `./fastq/` first.

```bash
nextflow run . -profile test,docker,cpu
```

With your own data, give it the reference and the reads.

```bash
nextflow run . -profile docker,gpu \
  --rna_fastq 'fastq/*_{1,2}.fastq.gz' \
  --gtf ref/annotation.gtf --fasta ref/genome.fa \
  --outdir results
```

- RNA-seq must be paired-end. For single-end RNA-seq, align it yourself and pass `--rna_bam`.
- `--fasta` must be uncompressed, with a `samtools faidx` index beside it.
- It predicts every isoform of the genes at 5 TPM or more. salmon quantifies the same alignment,
  and `--min_tpm` sets the cutoff; `--tpm_level transcript` applies it to each isoform's own TPM
  instead. The default comes from a sweep on the tutorial sample, on the score page. TPM is
  normalised over the reference, so on the chr22 `test` profile it runs about 47 times higher
  than genome-wide and the cutoff admits more transcripts than it would on a whole genome.
- On 16 CPU threads `attn` predicts about 8,500 chr22 transcripts in 32 min, and `mamba4` is
  6.3x slower on CPU, so a whole-genome run belongs on a GPU, or on `attn`.

If you already have a STAR transcriptome BAM, pass it instead and the pipeline skips the
alignment and the STAR index.

```bash
nextflow run . -profile docker,cpu \
  --rna_bam rna.toTranscriptome.bam --gtf ref/annotation.gtf --fasta ref/genome.fa
```

On SLURM, add the `slurm` profile. Nextflow submits one job per process.

```bash
nextflow run . -profile slurm,singularity,gpu \
  --rna_fastq 'fastq/*_{1,2}.fastq.gz' --gtf ref/annotation.gtf --fasta ref/genome.fa \
  --outdir results -resume
```

## Profiles

| profile | effect |
|---|---|
| `docker` / `singularity` | run every process in the container |
| `local` | no container: use tools already on `PATH`, plus this checkout's `scripts/` |
| `cpu` | `--device cpu`, and the CPU image, which has no `mamba_ssm` |
| `gpu` | `--device cuda`, the GPU image, one accelerator per predict task |
| `slurm` | submit each process to the scheduler |
| `test` | chr22 with the tutorial's RNA-seq accession |

## Parameters worth setting

| parameter | default | note |
|---|---|---|
| `--arch` | `attn,mamba4` | on GPU mamba4 is the faster of the two; on CPU it is 6.3x slower |
| `--min_tpm` | `5` | salmon TPM cutoff for the transcripts to predict |
| `--tpm_level` | `gene` | `gene`: every isoform of a gene at the cutoff; `transcript`: each isoform's own TPM |
| `--chrom` | `chr22` | only when the pipeline fetches the reference; `all` for the primary assembly, which needs about 32 GB to index |
| `--max_tx_length` | `10000` | the released checkpoints never saw a longer transcript |
| `--pred_scale` | `0.05` | the Poisson dial on the second calling arm |
| `--star_index` | `null` | reuse an existing index instead of building one |

`nextflow run . --help` prints a shorter summary; `nextflow.config` lists every parameter.

## What it writes

```
results/
  reference/     GTF, genome FASTA, transcript FASTA, STAR index
  align/         transcriptome BAM and STAR log
  logs/          cutadapt log
  coverage/      coverage.hd5, per-nucleotide RNA-seq coverage
  expression/    salmon quant.sf, and predicted_tx.txt, the transcripts passing the TPM cutoff
  pack/          tx_order.txt, offsets, lengths, coverage, expressed_tx.txt, ORF track
  weights/       the checkpoints, when fetched from Hugging Face
  pred/          <arch>/pred_profiles.npz and the per-arch timing
  calls/         RiboCode pred_preddepth and Poisson arms, and Ribo-TISH
```

Both calling arms are always written. Quoting the deterministic arm without the Poisson one, or
the reverse, gives a number better than the method is.
