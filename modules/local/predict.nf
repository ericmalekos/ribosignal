process FETCH_WEIGHTS {
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    val repo

    output:
    path 'weights', emit: weights

    script:
    """
    python -c "from huggingface_hub import snapshot_download; snapshot_download('${repo}', local_dir='weights')"
    cd weights && sha256sum -c SHA256SUMS
    """
}

process PREDICT {
    label 'predict'
    publishDir "${params.outdir}/pred", mode: params.publish_mode

    input:
    val arch
    path pack
    path track
    path tx_fasta
    path weights
    path tx_list
    val device

    output:
    tuple val(arch), path("${arch}/pred_profiles.npz"), emit: profiles
    path "${arch}_seconds.tsv",                         emit: timing

    script:
    """
    export RIBO_ONEHOT_FASTA=\$(readlink -f ${tx_fasta})
    export RIBO_ORF_TRACK=\$(readlink -f ${track})
    # Pin the thread count. torch otherwise takes every core it can see, which makes any
    # reported CPU time unreproducible; on a 160-core box it ran 84 threads unasked.
    export OMP_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus} NUMEXPR_NUM_THREADS=${task.cpus}

    T0=\$SECONDS
    # --tx_list: the transcripts at or above the TPM cutoff (SALMON_QUANT). Without it the
    # script falls back to an observed-P-site floor, which a pack with no Ribo-seq never meets.
    dump_pred_profiles.py --run ${weights}/ --arch ${arch} --pack ${pack} \\
        --tx_list ${tx_list} --out ${arch} --device ${device}
    printf '%s\\t%s\\t%s\\t%s\\n' "${arch}" "${device}" "${task.cpus}" "\$((SECONDS-T0))" \\
        > ${arch}_seconds.tsv
    cat ${arch}_seconds.tsv
    """
}
