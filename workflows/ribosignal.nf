include { FETCH_REFERENCE; TRANSCRIPT_FASTA; STAR_INDEX } from '../modules/local/reference'
include { MEASURE_ADAPTER; TRIM_RNA; TRIM_RIBO; ALIGN_RNA; ALIGN_RIBO } from '../modules/local/align'
include { RNA_COVERAGE; RIBO_PSITES; BUILD_PACK; ORF_TRACK } from '../modules/local/pack'
include { FETCH_WEIGHTS; PREDICT } from '../modules/local/predict'
include { RIBOCODE_ANNOT; RIBOCODE_CALL; RIBOTISH_CALL; SCORE } from '../modules/local/orfcall'

workflow RIBOSIGNAL {

    // ---- reference -------------------------------------------------------------------
    if (params.gtf && params.fasta) {
        ch_ref = Channel.fromPath(params.fasta, checkIfExists: true)
                    .map { fa -> tuple(file(params.gtf, checkIfExists: true), fa, file("${fa}.fai")) }
    } else {
        ch_ref = FETCH_REFERENCE(params.gencode_release, params.chrom).ref
    }

    ch_tx = TRANSCRIPT_FASTA(ch_ref)

    // Build the index only if something will actually align against it. Supplying both BAMs
    // and still indexing costs an hour and ~32 GB on a whole genome, for nothing.
    def needs_index = !(params.rna_bam && params.ribo_bam)
    ch_index = params.star_index ? Channel.fromPath(params.star_index, checkIfExists: true)
             : needs_index       ? STAR_INDEX(ch_ref).index
             :                     Channel.empty()

    // ---- RNA-seq: FASTQ or a transcriptome BAM straight in ---------------------------
    if (params.rna_bam) {
        ch_rna_bam = Channel.fromPath(params.rna_bam, checkIfExists: true)
                        .map { b -> tuple(b.simpleName, b) }
    } else {
        ch_rna = Channel.fromFilePairs(params.rna_fastq, checkIfExists: true)
        ch_rna_bam = ALIGN_RNA(TRIM_RNA(ch_rna).trimmed, ch_index.first()).bam
    }

    // ---- Ribo-seq --------------------------------------------------------------------
    if (params.ribo_bam) {
        ch_ribo_bam = Channel.fromPath(params.ribo_bam, checkIfExists: true)
                        .map { b -> tuple(b.simpleName, b) }
    } else {
        ch_ribo = Channel.fromPath(params.ribo_fastq, checkIfExists: true)
                     .map { f -> tuple(f.simpleName.replaceAll(/\.fastq$/, ''), f) }

        // ribo_adapter: null MEASURES it, '' asserts already-trimmed, a string forces it.
        // An empty measured answer is a real answer, not a failure to measure.
        if (params.ribo_adapter == null) {
            ch_adapter = MEASURE_ADAPTER(ch_ribo.map { s, f -> tuple('ribo', f) })
                            .adapter.map { assay, a -> a.trim() }
        } else {
            ch_adapter = Channel.value(params.ribo_adapter)
        }
        ch_ribo_bam = ALIGN_RIBO(TRIM_RIBO(ch_ribo, ch_adapter).trimmed, ch_index.first()).bam
    }

    // ---- pack ------------------------------------------------------------------------
    // The coverage guard's whole-transcriptome floor of 5,000 fires spuriously on a
    // chromosome subset, which holds far fewer transcripts. Scale it, do not remove it.
    ch_floor = ch_tx.map { fa -> Math.max(100, fa.countFasta().intdiv(10)) }

    // One sample per assay. With several, the two per-nt channels would pair positionally and
    // BUILD_PACK would quietly make a pack from one coverage file and one psites file, which
    // looks like a normal result. Pooling replicates is a separate step (scripts/prepare/).
    ch_rna_bam.count().subscribe { n ->
        if (n > 1) error "ribosignal takes ONE RNA-seq sample; got ${n}. Pool replicates first."
    }
    ch_ribo_bam.count().subscribe { n ->
        if (n > 1) error "ribosignal takes ONE Ribo-seq sample; got ${n}. Pool replicates first."
    }

    ch_cov    = RNA_COVERAGE(ch_rna_bam, ch_floor.first())
    ch_psites = RIBO_PSITES(ch_ribo_bam, ch_ref.first()).psites
    ch_pack   = BUILD_PACK(ch_tx.first(), ch_cov, ch_psites, params.max_tx_length)
    ch_track  = ORF_TRACK(ch_pack, ch_tx.first())

    // ---- predict ---------------------------------------------------------------------
    ch_weights = file(params.weights).isDirectory()
                    ? Channel.fromPath(params.weights, checkIfExists: true)
                    : FETCH_WEIGHTS(params.weights).weights

    ch_arch = Channel.fromList(params.arch.tokenize(',')*.trim())
    // .first() turns each single-item QUEUE channel into a VALUE channel. Without it the
    // first architecture consumes the pack, track, FASTA and weights and the second never
    // runs -- and with --arch attn alone the pipeline looks correct.
    ch_pred = PREDICT(ch_arch, ch_pack.first(), ch_track.first(), ch_tx.first(),
                      ch_weights.first(), params.device)

    // ---- call ------------------------------------------------------------------------
    ch_annot = RIBOCODE_ANNOT(ch_ref)
    ch_calls = RIBOCODE_CALL(ch_pred.profiles, ch_annot.first(),
                             params.min_aa, params.pval, params.pred_scale)

    if (!params.skip_ribotish) {
        RIBOTISH_CALL(ch_pred.profiles, ch_ref.first())
    }

    // ---- score -----------------------------------------------------------------------
    // Both sides stage DIRECTORIES named for the architecture, which is the layout
    // score_demo.py expects: pred/<arch>/pred_profiles.npz, calls/<arch>/real_collapsed.txt
    // and calls/<arch>_poisson/pred_preddepth_collapsed.txt.
    SCORE(
        ch_pred.dir.collect(),
        ch_calls.calls.mix(ch_calls.poisson).map { arch, d -> d }.collect()
    )
}
