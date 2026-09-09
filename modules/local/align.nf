// Adapter measurement, trimming, and the two alignments. RNA-seq and Ribo-seq take opposite
// settings on purpose; the comments say which, because getting them backwards is silent.

process MEASURE_ADAPTER {
    publishDir "${params.outdir}/logs", mode: params.publish_mode

    input:
    tuple val(assay), path(reads)

    output:
    tuple val(assay), path("adapter_${assay}.txt"), emit: report
    tuple val(assay), stdout,                       emit: adapter

    script:
    def first = reads instanceof List ? reads[0] : reads
    """
    measure_adapter.py ${first} > adapter_${assay}.txt
    cat adapter_${assay}.txt >&2
    # Empty stdout is the answer for an already-trimmed library, and is not an error.
    awk '/Use `cutadapt -a /{print \$4; found=1} END{if(!found) printf ""}' adapter_${assay}.txt | tr -d '\\n'
    """
}

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
    # footprints and wrong for RNA. On this library the footprint settings kept 16.4% of
    # reads where the correct call kept 99.0%.
    cutadapt --trim-n -m 20 -j ${task.cpus} -o rna_t1.fq.gz -p rna_t2.fq.gz \\
        ${reads[0]} ${reads[1]} > cutadapt_rna.log 2>&1
    """
}

process TRIM_RIBO {
    label 'alignment'
    publishDir "${params.outdir}/logs", mode: params.publish_mode, pattern: '*.log'

    input:
    tuple val(sample), path(reads)
    val adapter

    output:
    tuple val(sample), path('ribo_t.fq.gz'), emit: trimmed
    path 'cutadapt_ribo.log',                emit: log

    script:
    // An empty adapter is a legitimate measured answer. Passing -a --discard-untrimmed to an
    // already-trimmed library keeps only the reads carrying an adapter it never had.
    def trim = adapter ? "-a ${adapter} -m 20 -M 40 --discard-untrimmed" : "--trim-n -m 20 -M 40"
    """
    echo "ribo adapter: '${adapter ?: '<none detected>'}'"
    cutadapt ${trim} -j ${task.cpus} -o ribo_t.fq.gz ${reads} > cutadapt_ribo.log 2>&1
    kept=\$(awk -F'[(),]' '/Reads written/{print \$2}' cutadapt_ribo.log)
    echo "kept: \$kept"
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

process ALIGN_RIBO {
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
    # Unique alignments only: RiboCode and Ribo-TISH do not drop multimappers themselves, so
    # keeping them inflates P-site counts and breaks the periodicity test.
    # EndToEnd because soft-clipping shifts the inferred P-site.
    STAR --genomeDir ${index} --readFilesIn ${reads} --readFilesCommand zcat \\
         --runThreadN ${task.cpus} --outSAMtype None --quantMode TranscriptomeSAM \\
         --outFilterMultimapNmax 1 --alignEndsType EndToEnd \\
         --outSAMattributes NH HI AS nM --outFileNamePrefix ${sample}.
    n=\$(awk -F'\\t' '/Number of input reads/{gsub(/ /,"",\$2);print \$2}' ${sample}.Log.final.out)
    [ "\${n:-0}" -gt 0 ] || { echo "ABORT: STAR read 0 reads" >&2; exit 1; }
    samtools quickcheck ${sample}.Aligned.toTranscriptome.out.bam
    """
}
