// ORF calling, two callers and two arms.
//
//   real            the caller on OBSERVED P-sites; the reference the others are judged against
//   pred_preddepth  fully de novo: predicted shape AND predicted depth, no Ribo-seq used
//   + Poisson       the same arm through a Poisson dial, trading CDS recall for non-canonical
//                   precision (about 6 points for +0.10)
//
// Both arms are always reported. Quoting one without the other has repeatedly produced a
// number that looked better than the method is.

process RIBOCODE_ANNOT {
    publishDir "${params.outdir}/calls", mode: params.publish_mode

    input:
    tuple path(gtf), path(fa), path(fai)

    output:
    path 'annot', emit: annot

    script:
    """
    prepare_transcripts -g ${gtf} -f ${fa} -o annot
    [ -s annot/transcripts.pickle ] || { echo "ABORT: no transcripts.pickle" >&2; exit 1; }
    """
}

process RIBOCODE_CALL {
    label 'big_mem'
    publishDir "${params.outdir}/calls", mode: params.publish_mode

    input:
    tuple val(arch), path(profiles)
    path annot
    val min_aa
    val pval
    val pred_scale

    output:
    tuple val(arch), path("${arch}"),         emit: calls
    tuple val(arch), path("${arch}_poisson"), emit: poisson

    script:
    """
    ribocode_dropin.py --profiles ${profiles} --annot ${annot} --variant real \\
        --out ${arch} --min_aa ${min_aa} --pval ${pval}
    ribocode_dropin.py --profiles ${profiles} --annot ${annot} --variant pred_preddepth \\
        --out ${arch} --min_aa ${min_aa} --pval ${pval}
    ribocode_dropin.py --profiles ${profiles} --annot ${annot} --variant pred_preddepth \\
        --pred_poisson --pred_scale ${pred_scale} \\
        --out ${arch}_poisson --min_aa ${min_aa} --pval ${pval}
    """
}

process RIBOTISH_CALL {
    label 'big_mem'
    publishDir "${params.outdir}/calls", mode: params.publish_mode

    input:
    tuple val(arch), path(profiles)
    tuple path(gtf), path(fa), path(fai)

    output:
    tuple val(arch), path("ribotish_${arch}"), emit: calls
    path "scored_${arch}.gtf",                 emit: gtf

    script:
    """
    # --emit-gtf makes the wrapper restrict the annotation to exactly the transcripts it
    # profiled, keeping only genes all of whose transcripts survived. Ribo-TISH skips its bam
    # path only when a gene is fully covered, so filtering on the pack instead is wrong: the
    # pack holds every transcript while the prediction covers only those with signal, and
    # Ribo-TISH then goes looking for a bam that is not there.
    ribotish_dropin.py --profiles ${profiles} --gtf ${gtf} --genome ${fa} \\
        --emit-gtf scored_${arch}.gtf --variant pred_preddepth --out ribotish_${arch} \\
        --longest --minaalen 5 --fpth ${params.pval} --numproc ${task.cpus}
    """
}

process SCORE {
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path 'pred/*'
    path 'calls/*'

    output:
    path 'RESULTS.json', emit: json
    path 'RESULTS.txt',  emit: txt

    script:
    """
    score_demo.py --pred-dir pred --calls-dir calls --out RESULTS.json | tee RESULTS.txt

    # score_demo.py prints "no pred_profiles.npz found" and still exits 0, so a stage-in
    # mistake produced a COMPLETED run whose profile half was silently empty. Half a result
    # that reports success is worse than a failure, so assert both halves are populated.
    python -c "\\
import json, sys; \\
r = json.load(open('RESULTS.json')); \\
m = [k for k in ('profile', 'calls') if not r.get(k)]; \\
sys.exit('ABORT: RESULTS.json has no ' + ' and no '.join(m) + ' section. Inputs did not \\
stage as score_demo.py expects: pred/<arch>/pred_profiles.npz and \\
calls/<arch>/real_collapsed.txt.') if m else \\
print('scored %d profile arm(s), %d call arm(s)' % (len(r['profile']), len(r['calls'])))"
    """
}
