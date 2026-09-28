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

process BUILD_PACK {
    label 'big_mem'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path tx_fasta
    path coverage
    val max_length
    val min_coverage

    output:
    path 'pack', emit: pack

    script:
    """
    # --max-length is not cosmetic. The released checkpoints were trained on a universe
    # capped at 10,000 nt, and attention is O(L^2): one 37,852 nt chr22 transcript asks for
    # 42.7 GiB across 8 heads and takes the whole run down.
    # No --psites: the pipeline reads no Ribo-seq. build_pack.py then writes expressed_tx.txt,
    # the transcripts with total RNA coverage >= --min-coverage, which PREDICT scores.
    build_pack.py --fasta ${tx_fasta} --coverage ${coverage} \\
        --out pack --max-length ${max_length} --min-coverage ${min_coverage}
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
