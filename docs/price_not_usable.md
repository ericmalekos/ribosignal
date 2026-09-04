# PRICE cannot be run on predicted Ribo-seq profiles

Recorded 2026-09-04. PRICE was tested alongside RiboCode, Ribo-TISH, RiboTaper and ribotricer and
was **dropped**. This file exists so the next person does not repeat the attempt.

## The finding

**PRICE fits its cleavage model to the read MISMATCH distribution, so it needs real read
sequences. A predicted profile has none, and no amount of packaging supplies them.**

PRICE 1.0.4 (GEDI 1.0.6a) runs correctly on real reads. It was run on the same Hepatocytes
Ribo-seq the other four callers saw, projected from the transcriptome bams to the genome by
`scripts/orfcallers/tx_bam_to_genome.py`, and produced `chr1.orfs.tsv` and `chr1.orfs.cit` in
about 15 minutes. It fails on every arm built from a density array.

## The evidence

PRICE writes `<prefix>.estimateData`, a table of read counts stratified by frame and by number of
mismatches: `F0L0/F1L0/F2L0` are zero-mismatch reads, `F0L1/F1L1/F2L1` one-mismatch.

| arm | read lengths | L1 (one-mismatch) columns |
|---|---|---|
| real projected reads | 21 to 30 | populated, e.g. 1044 / 1210 / 988 at length 22 |
| synthesized from a density array | 29 only | **every value 0.00** |

With no mismatch signal the estimator aborts:

    Could not estimate model, there is something wrong with mismatch information
    (did you forget to include the MD attribute?)
      at gedi.riboseq.cleavage.CleavageModelEstimator.estimateBoth(CleavageModelEstimator.java:609)

## Two fixes that were tried and did not work

**Adding the MD tag.** The error message names `MD`, so `samtools calmd` was run against the
reference. Confirmed present afterwards. PRICE failed identically. The message is misleading: the
tag was missing, but supplying it does not supply mismatches.

**Giving the synthetic reads real sequence.** `synth_ribo_bam.py --ref-seq` was added, so each
emitted read carries the reference sequence at its projected position rather than a run of `N`.
Verified: `MD:Z:29` on 20,000 of 20,000 reads, with the P-site self-check still at 0 mismatches.
PRICE failed identically, because reference-matching reads produce `L1 = 0` by construction. This
is the step that makes the diagnosis certain: the obstacle is the absence of mismatches, not the
absence of the tag.

## What was NOT done, deliberately

**Mismatches were not injected.** Adding synthetic sequencing errors would let PRICE fit a model,
but that model would be an artifact of whatever error rate was chosen rather than a property of
the data, and every downstream ORF call would inherit it.

**Model reuse was attempted and abandoned.** PRICE has no option to load a pre-estimated model,
even under `-hhh`. Seeding `chr1.model`, `chr1.estimateData` and `chr1.maxPos` from the real
projected arm into the predicted arm's output directory, in the hope that GEDI's pipeline would
skip the estimation step, was started (`scripts/orfcallers/price_modelreuse.sbatch`) and cancelled
before completing. Anyone resuming this should know it is the only remaining avenue, and that even
if it runs, it applies a cleavage model fitted to a 21-to-30 nt library to a uniform 29 nt one.

## Where this leaves the caller set

Four of five callers work on predicted profiles:

| caller | mechanism | works on predictions |
|---|---|---|
| RiboCode | native density injection | yes |
| Ribo-TISH | `predict --inprofile` | yes |
| RiboTaper | synthesized genome bam | yes |
| ribotricer | synthesized genome bam | yes |
| **PRICE** | genome bam, **plus real mismatch structure** | **no** |

PRICE remains usable on real Ribo-seq, including the transcriptome bams this project keeps, via
the projection tool. It is only the predicted arms it cannot score, which is what the comparison
needed it for.
