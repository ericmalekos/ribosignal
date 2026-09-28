#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RIBOSIGNAL } from './workflows/ribosignal'

def helpMessage() {
    log.info """
    ribosignal ${workflow.manifest.version}

    Predict per-nucleotide Ribo-seq P-site profiles from transcript sequence and RNA-seq
    coverage, then call ORFs from the prediction. Takes RNA-seq only; no Ribo-seq is read.

    Quick start (chr22; fetches the reference, reads the tutorial's RNA-seq FASTQs from ./fastq/):

      nextflow run . -profile test,docker,cpu

    Bring your own data (paired-end RNA-seq, and an uncompressed FASTA with a samtools faidx
    index beside it):

      nextflow run . -profile docker,gpu \\
        --rna_fastq 'fastq/*_{1,2}.fastq.gz' \\
        --gtf ref/annotation.gtf --fasta ref/genome.fa \\
        --outdir results

    Already have a STAR transcriptome BAM:

      nextflow run . -profile docker,cpu \\
        --rna_bam rna.toTranscriptome.bam --gtf ref/annotation.gtf --fasta ref/genome.fa

    Key parameters (see nextflow.config for the rest):

      --arch          attn, mamba4, or both (default '${params.arch}').
                      On GPU mamba4 is faster than attn; on CPU it is 6.3x slower.
      --device        cpu or cuda. The cpu and gpu profiles set this.
      --min_coverage  predict transcripts whose total RNA coverage reaches this
                      (default ${params.min_coverage}).
      --chrom         chr22 (default) or all, when the reference is fetched. 'all' needs
                      ~32 GB to index.
      --max_tx_length ${params.max_tx_length}. The checkpoints never saw a longer transcript.
      --pred_scale    Poisson dial for the second calling arm (default ${params.pred_scale}).

    Profiles: docker, singularity, local, cpu, gpu, slurm, test
    """.stripIndent()
}

// Nextflow 26 forbids top-level statements, so the help gate lives inside the workflow.
workflow {
    if (params.help) {
        helpMessage()
    } else {
        RIBOSIGNAL()
    }
}
