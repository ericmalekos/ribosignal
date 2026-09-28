// RNA-seq BAM -> per-nucleotide coverage HDF5 -> the pack the model reads, plus the ORF track.

process RNA_COVERAGE {
    label 'big_mem'
    // NOT outdir/pack: BUILD_PACK publishes a directory of that name, and publishing a
    // directory replaces what is there, so coverage.hd5 was deleted on every fresh run.
    publishDir "${params.outdir}/coverage", mode: params.publish_mode

    input:
    tuple val(sample), path(bam)
    val floor

    output:
    path 'coverage.hd5', emit: coverage

    script:
    """
    rnaseq_coverage.py --bam ${bam} --out coverage.hd5 --sample ${sample} \\
        --min-covered-tx ${floor}
    """
}

process SALMON_QUANT {
    label 'big_mem'
    publishDir "${params.outdir}/expression", mode: params.publish_mode

    input:
    tuple val(sample), path(bam)
    path tx_fasta
    tuple path(gtf), path(fa), path(fai)
    val min_tpm
    val tpm_level

    output:
    path 'quant.sf',         emit: quant
    path 'predicted_tx.txt', emit: tx_list

    script:
    """
    # Which transcripts to predict. With no Ribo-seq there is no observed-P-site floor, so the
    # selection is by expression: salmon alignment mode on STAR's transcriptome BAM, whose EM
    # splits reads shared between isoforms, which per-transcript coverage cannot. TPM is
    # normalised over the transcripts in the reference, so on a one-chromosome reference (the
    # chr22 test) it runs far higher than a whole-transcriptome TPM.
    salmon quant -t ${tx_fasta} -l A -a ${bam} -p ${task.cpus} -o salmon > salmon.log 2>&1 \\
        || { cat salmon.log >&2; exit 1; }
    cp salmon/quant.sf quant.sf

    # gene: every isoform of a gene whose summed TPM reaches the cutoff (the default; keeps the
    # minor isoforms that carry many non-canonical ORFs). transcript: the isoform's own TPM.
    awk -F'\\t' '\$3 == "transcript" {
        match(\$9, /gene_id "[^"]+"/);       g = substr(\$9, RSTART + 9,  RLENGTH - 10)
        match(\$9, /transcript_id "[^"]+"/); t = substr(\$9, RSTART + 15, RLENGTH - 16)
        print t "\\t" g }' ${gtf} > tx2gene.tsv
    awk -F'\\t' -v m=${min_tpm} -v lvl=${tpm_level} '
        FNR == NR { gene[\$1] = \$2; next }
        FNR > 1   { tpm[\$1] = \$4; if (\$1 in gene) gsum[gene[\$1]] += \$4 }
        END { for (t in tpm) {
                  v = (lvl == "gene" && (t in gene)) ? gsum[gene[t]] : tpm[t]
                  if (v >= m) print t } }' tx2gene.tsv quant.sf | sort > predicted_tx.txt
    n=\$(wc -l < predicted_tx.txt)
    echo "transcripts selected at ${tpm_level} TPM >= ${min_tpm}: \$n of \$(( \$(wc -l < quant.sf) - 1 ))"
    [ "\$n" -gt 0 ] || { echo "ABORT: no transcript reaches ${tpm_level} TPM ${min_tpm}" >&2; exit 1; }
    """
}

process BUILD_PACK {
    label 'big_mem'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path tx_fasta
    path coverage
    val max_length

    output:
    path 'pack', emit: pack

    script:
    """
    # --max-length is not cosmetic. The released checkpoints were trained on a universe
    # capped at 10,000 nt, and attention is O(L^2): one 37,852 nt chr22 transcript asks for
    # 42.7 GiB across 8 heads and takes the whole run down.
    # No --psites: the pipeline reads no Ribo-seq, so the pack's P-site target is all zeros.
    # Which transcripts get predicted is SALMON_QUANT's TPM cutoff, not this pack.
    build_pack.py --fasta ${tx_fasta} --coverage ${coverage} \\
        --out pack --max-length ${max_length}
    """
}

process ORF_TRACK {
    label 'big_mem'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path pack
    path tx_fasta

    output:
    path "${pack}/orf_track_v2_nokozak.npy", emit: track

    script:
    """
    # --kozak none on purpose: the four-arm ablation tied in distribution (0.642 to 0.646)
    # and the explicit Kozak gate was redundant or harmful cross-tissue. Both released
    # checkpoints are trained without it, so inference must match.
    build_orf_track.py --pack ${pack} --fasta ${tx_fasta} --mode ext --kozak none
    """
}
