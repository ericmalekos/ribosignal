// RNA-seq trimming and alignment. The settings are the RNA-seq ones on purpose: footprint
// settings (--maximum-length, --discard-untrimmed, unique mappers only) are wrong here, and
// getting them backwards is silent.

process TRIM_RNA {
    label 'alignment'
    publishDir "${params.outdir}/logs", mode: params.publish_mode, pattern: '*.log'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('rna_t{1,2}.fq.gz'), emit: trimmed
    path 'cutadapt_rna.log',                     emit: log

    script:
    """
    # -m 20 only. NOT --maximum-length and NOT --discard-untrimmed: both are right for
    # footprints and wrong for RNA. On one RNA-seq library the footprint settings kept 16.4%
    # of reads.
    cutadapt --trim-n -m 20 -j ${task.cpus} -o rna_t1.fq.gz -p rna_t2.fq.gz \\
        ${reads[0]} ${reads[1]} > cutadapt_rna.log 2>&1
    """
}

process ALIGN_RNA {
    label 'alignment'
    publishDir "${params.outdir}/align", mode: params.publish_mode

    input:
    tuple val(sample), path(reads)
    path index

    output:
    tuple val(sample), path("${sample}.Aligned.toTranscriptome.out.bam"), emit: bam
    path "${sample}.Log.final.out",                                       emit: log

    script:
    """
    # RNA keeps multimappers up to 10: the coverage channel wants the depth, and dropping
    # them biases against paralogues and repeat-adjacent transcripts.
    STAR --genomeDir ${index} --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand zcat \\
         --runThreadN ${task.cpus} --outSAMtype None --quantMode TranscriptomeSAM \\
         --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM \\
         --outFileNamePrefix ${sample}.
    # STAR can exit 0 having read nothing, leaving a valid-looking empty BAM.
    n=\$(awk -F'\\t' '/Number of input reads/{gsub(/ /,"",\$2);print \$2}' ${sample}.Log.final.out)
    [ "\${n:-0}" -gt 0 ] || { echo "ABORT: STAR read 0 reads" >&2; exit 1; }
    samtools quickcheck ${sample}.Aligned.toTranscriptome.out.bam
    """
}
