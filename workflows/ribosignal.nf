include { FETCH_REFERENCE; TRANSCRIPT_FASTA; STAR_INDEX } from '../modules/local/reference'
include { TRIM_RNA; ALIGN_RNA } from '../modules/local/align'
include { RNA_COVERAGE; BUILD_PACK; ORF_TRACK } from '../modules/local/pack'
include { FETCH_WEIGHTS; PREDICT } from '../modules/local/predict'
include { RIBOCODE_ANNOT; RIBOCODE_CALL; RIBOTISH_CALL } from '../modules/local/orfcall'

// RNA-seq in, predicted P-site profiles and ORF calls out. No Ribo-seq is read at any step:
// scoring a prediction against real Ribo-seq is the tutorial's job (scripts/demo/score_demo.py),
// not the pipeline's.

workflow RIBOSIGNAL {

    // ---- reference -------------------------------------------------------------------
    if (params.gtf && params.fasta) {
        ch_ref = Channel.fromPath(params.fasta, checkIfExists: true)
                    .map { fa -> tuple(file(params.gtf, checkIfExists: true), fa,
                                       file("${fa}.fai", checkIfExists: true)) }
    } else {
        ch_ref = FETCH_REFERENCE(params.gencode_release, params.chrom).ref
    }

    ch_tx = TRANSCRIPT_FASTA(ch_ref)

    // ---- RNA-seq: FASTQ, or a transcriptome BAM straight in ----------------------------
    // Build the index only if something will align against it: a whole-genome index costs an
    // hour and ~32 GB.
    if (params.rna_bam) {
        ch_rna_bam = Channel.fromPath(params.rna_bam, checkIfExists: true)
                        .map { b -> tuple(b.simpleName, b) }
    } else {
        if (!params.rna_fastq) error "Give --rna_fastq (paired FASTQ glob) or --rna_bam."
        ch_index = params.star_index ? Channel.fromPath(params.star_index, checkIfExists: true)
                 :                     STAR_INDEX(ch_ref).index
        ch_rna = Channel.fromFilePairs(params.rna_fastq, checkIfExists: true)
        ch_rna_bam = ALIGN_RNA(TRIM_RNA(ch_rna).trimmed, ch_index.first()).bam
    }

    // One sample. With several, BUILD_PACK would quietly use one of them, which looks like a
    // normal result. Pool replicates before the pipeline.
    ch_rna_bam.count().subscribe { n ->
        if (n > 1) error "ribosignal takes ONE RNA-seq sample; got ${n}. Pool replicates first."
    }

    // ---- pack ------------------------------------------------------------------------
    // The coverage guard's whole-transcriptome floor of 5,000 fires spuriously on a
    // chromosome subset, which holds far fewer transcripts. Scale it, do not remove it.
    ch_floor = ch_tx.map { fa -> Math.max(100, fa.countFasta().intdiv(10)) }

    ch_cov   = RNA_COVERAGE(ch_rna_bam, ch_floor.first())
    ch_pack  = BUILD_PACK(ch_tx.first(), ch_cov, params.max_tx_length, params.min_coverage)
    ch_track = ORF_TRACK(ch_pack, ch_tx.first())

    // ---- predict ---------------------------------------------------------------------
    // Both branches are VALUE channels, so every architecture reads the same weights.
    ch_weights = file(params.weights).isDirectory()
                    ? Channel.value(file(params.weights, checkIfExists: true))
                    : FETCH_WEIGHTS(params.weights).weights

    ch_arch = Channel.fromList(params.arch.tokenize(',')*.trim())
    // .first() turns each single-item QUEUE channel into a VALUE channel. Without it the
    // first architecture consumes the pack, track and FASTA and the second never runs -- and
    // with --arch attn alone the pipeline looks correct.
    ch_pred = PREDICT(ch_arch, ch_pack.first(), ch_track.first(), ch_tx.first(),
                      ch_weights, params.device)

    // ---- call ------------------------------------------------------------------------
    ch_annot = RIBOCODE_ANNOT(ch_ref)
    RIBOCODE_CALL(ch_pred.profiles, ch_annot.first(),
                  params.min_aa, params.pval, params.pred_scale)

    if (!params.skip_ribotish) {
        RIBOTISH_CALL(ch_pred.profiles, ch_ref.first())
    }
}
