// BAM -> per-nucleotide HDF5 -> the pack the model reads, plus the ORF track.

process RNA_COVERAGE {
    label 'big_mem'
    publishDir "${params.outdir}/pack", mode: params.publish_mode

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

process RIBO_PSITES {
    label 'big_mem'
    publishDir "${params.outdir}/pack", mode: params.publish_mode, pattern: '*.hd5'
    publishDir "${params.outdir}/logs", mode: params.publish_mode, pattern: '*.log'

    input:
    tuple val(sample), path(bam)
    tuple path(gtf), path(fa), path(fai)

    output:
    path 'psites.hd5',  emit: psites
    path 'psites.log',  emit: log

    script:
    """
    # The frame-0 fraction is the whole QC: 1/3 is noise, a real library reaches 0.5 to 0.7.
    # ribo_psites.py exits non-zero below --min-frame0, so a bad library fails here, loudly.
    ribo_psites.py --bam ${bam} --gtf ${gtf} --out psites.hd5 \\
        --sample ${sample} --lengths 25:35 2>&1 | tee psites.log
    grep -E "frame-0 fraction" psites.log || true
    """
}

process BUILD_PACK {
    label 'big_mem'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path tx_fasta
    path coverage
    path psites
    val max_length

    output:
    path 'pack', emit: pack

    script:
    """
    # --max-length is not cosmetic. The released checkpoints were trained on a universe
    # capped at 10,000 nt, and attention is O(L^2): one 37,852 nt chr22 transcript asks for
    # 42.7 GiB across 8 heads and takes the whole run down.
    build_pack.py --fasta ${tx_fasta} --coverage ${coverage} --psites ${psites} \\
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
