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
    val device

    output:
    tuple val(arch), path("${arch}/pred_profiles.npz"), emit: profiles
    // The directory too, not just the file: score_demo.py looks for
    // pred/<arch>/pred_profiles.npz, so staging the bare npz would both lose the arch and
    // collide on the second architecture.
    path "${arch}",                                     emit: dir
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
    dump_pred_profiles.py --run ${weights}/ --arch ${arch} --pack ${pack} \\
        --out ${arch} --device ${device}
    printf '%s\\t%s\\t%s\\t%s\\n' "${arch}" "${device}" "${task.cpus}" "\$((SECONDS-T0))" \\
        > ${arch}_seconds.tsv
    cat ${arch}_seconds.tsv
    """
}
