# Score against the real Ribo-seq

The run ends with a number rather than an assertion that it completed. `score_demo.py` reports
per-nucleotide profile correlation and ORF-call precision, recall and F1, predicted against
observed.

```bash
python $RIBO_SCRIPTS/demo/score_demo.py --pred-dir pred --calls-dir calls \
       --out RESULTS.json | tee RESULTS.txt
```

## What this run produces

chr22, 6,583 scored transcripts (those with at least 50 observed P-sites), Hepatocytes_1 (the
tissue held out of training). Check a run against this table; the next section explains why its
non-canonical columns overstate the model.

| architecture | arm | overall P | R | F1 | annotated F1 | novel P | novel R |
|---|---|--:|--:|--:|--:|--:|--:|
| `attn` | `pred_preddepth` | 0.723 | 0.893 | 0.799 | 0.980 | 0.576 | 0.779 |
| `attn` | `+ poisson 0.05` | 0.897 | 0.704 | 0.789 | 0.928 | 0.800 | 0.508 |
| `mamba4` | `pred_preddepth` | 0.745 | 0.889 | 0.811 | 0.980 | 0.595 | 0.770 |
| `mamba4` | `+ poisson 0.05` | 0.903 | 0.719 | 0.801 | 0.948 | 0.791 | 0.512 |

**Annotated-CDS F1** should be high; this is the class the annotation makes trustworthy, and 0.980
is what a working run looks like. **The Poisson arm is the dial**: it trades recall for
precision. Report both arms; quoting either alone gives a number better than the method is.

## The chr22-only alignment inflates the non-canonical numbers

Aligning to chr22 alone leaves reads from genes on other chromosomes nowhere to go but their chr22
look-alikes, mostly pseudogenes. Those then look expressed in the RNA-seq the model reads and
translated in the Ribo-seq it is scored against, so model and reference agree on signal that is
not there. 178 of the 729 real calls sit on genes whose genome-wide TPM is below 0.1, and 236 of
the 244 novel calls are gone when both assays are aligned to the whole genome: the transcripts
carrying them drop from a median of 814 P-sites to none.

Scored the clean way, with both assays aligned to the whole genome and then restricted to chr22
(the Ribo-seq also through the project's standard ncRNA and cross-gene read filters), on the 5,116
transcripts with at least 50 observed P-sites:

| architecture | arm | reference | P | R | F1 | annotated P / R | non-canonical P / R |
|---|---|--:|--:|--:|--:|--:|--:|
| `attn` | `pred_preddepth`, chr22-only | 729 | 0.723 | 0.893 | 0.799 | 0.968 / 0.992 | 0.518 / 0.773 |
| `attn` | `pred_preddepth`, whole genome | 381 | 0.717 | 0.945 | 0.815 | 0.957 / 0.997 | 0.281 / 0.714 |
| `attn` | `+ poisson 0.05`, chr22-only | 729 | 0.897 | 0.704 | 0.789 | 0.973 / 0.887 | 0.760 / 0.479 |
| `attn` | `+ poisson 0.05`, whole genome | 381 | 0.879 | 0.843 | 0.861 | 0.968 / 0.958 | 0.404 / 0.329 |
| `mamba4` | `pred_preddepth`, whole genome | 381 | 0.765 | 0.950 | 0.848 | 0.960 / 1.000 | 0.342 / 0.729 |
| `mamba4` | `+ poisson 0.05`, whole genome | 381 | 0.904 | 0.861 | 0.882 | 0.965 / 0.968 | 0.529 / 0.386 |

Annotated CDS calls hold up. Non-canonical precision falls by a third to a half: with a clean
reference, roughly one non-canonical call in three or four is supported on the deterministic arm,
and two in five to one in two on the Poisson arm. The clean reference holds 70 non-canonical calls,
only 9 of them novel, so these rates rest on small counts. The tutorial keeps the chr22-only
alignment because a whole-genome STAR index needs about 32 GB of RAM; treat its non-canonical
numbers as an upper bound.

## Profile correlation, and what it is not

```
attn     n=6583   pearson +0.2050   spearman(tie-corrected) +0.2915
         periodicity  pred 0.7918  obs 0.7289   (1/3 = none)
mamba4   n=6583   pearson +0.1923   spearman(tie-corrected) +0.2729
         periodicity  pred 0.7807  obs 0.7289   (1/3 = none)
```

Predicted periodicity exceeding observed is expected: the prediction is a smooth expectation, and
the single observed run it is compared against is itself noisy.

**The Pearson figure is well below the 0.6585 that `release/orf_v2_*/test_metrics.json` reports**,
and the report says so rather than glossing it. It is the same quantity, since `pred_flat` is
`softmax(logits)` and Pearson is scale-invariant, so the gap is real. It is not a transcript-count
effect and only partly a depth effect: binning this run by observed counts gives 0.07 at 50 to 100
counts, rising to 0.33 at 2,000 to 5,000, then flattening near 0.27 in the deepest bin. The demo
differs from the release evaluation in three ways that have not been separated: one RNA-seq run
rather than the pooled Hepatocytes pack, one Ribo-seq run, and RNA aligned at `mm10` where the
checkpoints were trained on `mm1` coverage. Treat the demo number as a working end-to-end check,
not as a reproduction of the released metric.

**Do not compare the periodicity line to `test_metrics.json`'s `period_*`.** They are different
quantities: this page reports the max-frame fraction, floored at 1/3, while training reports an
autocorrelation contrast (frame lags 3/6/9/12 minus off-frame lags). That is why the released
`period_obs_median` of 0.1457 sits below 1/3, which no fraction can.

## Without Ribo-seq to choose the transcripts

The Nextflow pipeline takes RNA-seq only, so it chooses transcripts by expression: salmon TPM on
the same alignment, either per gene (every isoform of a gene at the cutoff) or per transcript.
Same whole-genome-aligned sample and clean reference as above, each selection re-run through
RiboCode, F1 over all ORF classes:

| selection | transcripts | `attn` | `attn` + Poisson | `mamba4` | `mamba4` + Poisson |
|---|--:|--:|--:|--:|--:|
| >= 50 observed P-sites (needs Ribo-seq) | 5,116 | 0.815 | 0.861 | 0.848 | 0.882 |
| every expressed transcript | 8,194 | 0.460 | 0.663 | 0.467 | 0.656 |
| gene TPM >= 0.5 | 7,342 | 0.559 | 0.703 | 0.577 | 0.723 |
| gene TPM >= 1 | 7,035 | 0.602 | 0.732 | 0.626 | 0.739 |
| gene TPM >= 2 | 6,505 | 0.643 | 0.752 | 0.669 | 0.778 |
| **gene TPM >= 5 (default)** | 5,879 | 0.674 | 0.747 | 0.699 | 0.775 |
| gene TPM >= 10 | 5,275 | 0.660 | 0.715 | 0.676 | 0.737 |
| gene TPM >= 20 | 3,702 | 0.559 | 0.585 | 0.573 | 0.615 |
| transcript TPM >= 0.5 | 1,651 | 0.638 | 0.723 | 0.645 | 0.743 |
| transcript TPM >= 1 | 1,242 | 0.686 | 0.738 | 0.692 | 0.741 |
| transcript TPM >= 2 | 855 | 0.681 | 0.710 | 0.691 | 0.717 |
| transcript TPM >= 5 | 455 | 0.617 | 0.627 | 0.622 | 0.633 |

Gene TPM 5 is the default: the best or near-best F1 on all four arms. Transcript-level cutoffs
reach similar F1 but keep only the isoforms salmon credits, and lose about half the
non-canonical recall (for `attn`, 0.714 at gene TPM 5 against 0.343 at transcript TPM 1). Without
Ribo-seq, non-canonical precision is low: at the default it is 0.178 and 0.229 for the two `attn`
arms and 0.202 and 0.293 for `mamba4`, against 0.281 to 0.529 when the Ribo-seq chooses the
transcripts. Some of those calls may be real ORFs this Ribo-seq is too shallow to confirm.
