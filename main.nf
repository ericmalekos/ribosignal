#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RIBOSIGNAL } from './workflows/ribosignal'

def helpMessage() {
    log.info """
    ribosignal ${workflow.manifest.version}

    Predict per-nucleotide Ribo-seq P-site profiles from RNA-seq coverage, then call ORFs
    from the prediction. Mirrors scripts/demo/run_demo.sh.

    Quick start (chr22, downloads its own reference and data):

      nextflow run . -profile test,docker,cpu

    Bring your own data:

      nextflow run . -profile docker,gpu \\
        --rna_fastq  'fastq/*_{1,2}.fastq.gz' \\
        --ribo_fastq 'fastq/ribo.fastq.gz' \\
        --gtf ref/annotation.gtf --fasta ref/genome.fa \\
        --outdir results

    Already have STAR transcriptome BAMs:

      nextflow run . -profile docker,cpu \\
        --rna_bam rna.toTranscriptome.bam --ribo_bam ribo.toTranscriptome.bam \\
        --gtf ref/annotation.gtf --fasta ref/genome.fa

    Key parameters (see nextflow.config for the rest):

      --arch          attn, mamba4, or both (default '${params.arch}').
                      On GPU mamba4 is faster than attn; on CPU it is about 11x slower.
      --device        cpu or cuda. The cpu and gpu profiles set this.
      --ribo_adapter  null MEASURES it (default), '' asserts already-trimmed, a
                      string forces it. Never default an adapter you have not measured.
      --chrom         chr22 (default) or all. 'all' needs ~32 GB to index.
      --max_tx_length ${params.max_tx_length}. The checkpoints never saw a longer transcript.
      --pred_scale    Poisson dial for the second calling arm (default ${params.pred_scale}).

    Profiles: docker, singularity, cpu, gpu, slurm, test
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
