// Fetch GENCODE, subset to one chromosome, derive the transcript FASTA, build the STAR index.

process FETCH_REFERENCE {
    label 'big_mem'
    publishDir "${params.outdir}/reference", mode: params.publish_mode

    input:
    val release
    val chrom

    output:
    tuple path("${chrom}.gtf"), path("${chrom}.fa"), path("${chrom}.fa.fai"), emit: ref

    script:
    def base = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${release}"
    """
    curl -sSL -o gencode.gtf.gz "${base}/gencode.v${release}.annotation.gtf.gz"
    curl -sSL -o genome.fa.gz   "${base}/GRCh38.primary_assembly.genome.fa.gz"

    if [ "${chrom}" = all ]; then
        gunzip -c genome.fa.gz > ${chrom}.fa
        gunzip -c gencode.gtf.gz > ${chrom}.gtf
    else
        # Match on \$1 of the header so chr2 does not also pull chr22.
        gunzip -c genome.fa.gz | awk -v c=">${chrom}" '/^>/{p=(\$1==c)} p' > ${chrom}.fa
        gunzip -c gencode.gtf.gz | awk -F'\\t' -v c="${chrom}" '\$1==c' > ${chrom}.gtf
    fi
    samtools faidx ${chrom}.fa
    [ -s ${chrom}.fa ] || { echo "ABORT: no sequence for ${chrom}" >&2; exit 1; }
    """
}

process TRANSCRIPT_FASTA {
    publishDir "${params.outdir}/reference", mode: params.publish_mode

    input:
    tuple path(gtf), path(fa), path(fai)

    output:
    path 'transcripts.fa', emit: tx

    script:
    """
    gffread -w tx.raw -g ${fa} ${gtf}
    # gffread appends its own description after the id; the pack keys on the id alone.
    awk '/^>/{print \$1; next}{print}' tx.raw > transcripts.fa
    n=\$(grep -c '^>' transcripts.fa)
    echo "transcripts: \$n"
    [ "\$n" -gt 100 ] || { echo "ABORT: only \$n transcripts" >&2; exit 1; }
    """
}

process STAR_INDEX {
    label 'alignment'
    publishDir "${params.outdir}/reference", mode: params.publish_mode

    input:
    tuple path(gtf), path(fa), path(fai)

    output:
    path 'star_index', emit: index

    script:
    """
    # genomeSAindexNbases must scale with genome size; the default is right for a whole
    # genome and wastes memory (or fails) on a single chromosome.
    GLEN=\$(awk '{s+=\$2} END{print s}' ${fai})
    NB=\$(python -c "import math;print(min(14,int(math.log2(\$GLEN)/2-1)))")
    echo "genome \$GLEN bp, --genomeSAindexNbases \$NB"
    mkdir -p star_index
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star_index \\
         --genomeFastaFiles ${fa} --sjdbGTFfile ${gtf} --sjdbOverhang 75 \\
         --genomeSAindexNbases \$NB --outFileNamePrefix star_index/
    [ -s star_index/geneInfo.tab ] || { echo "ABORT: index has no transcriptome tables" >&2; exit 1; }
    """
}
