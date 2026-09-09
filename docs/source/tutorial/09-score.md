# Score against the real Ribo-seq

The run ends with a number rather than an assertion that it completed. `score_demo.py` reports
per-nucleotide profile correlation and ORF-call precision, recall and F1, predicted against
observed.

```bash
python $RIBO_SCRIPTS/demo/score_demo.py --pred-dir pred --calls-dir calls \
       --out RESULTS.json | tee RESULTS.txt
```

## What this run produces

chr22, 6,583 scored transcripts, Hepatocytes_1 (the tissue held out of training).

| architecture | arm | overall P | R | F1 | annotated F1 | novel P | novel R |
|---|---|--:|--:|--:|--:|--:|--:|
| `attn` | `pred_preddepth` | 0.723 | 0.893 | 0.799 | 0.980 | 0.576 | 0.779 |
| `attn` | `+ poisson 0.05` | 0.897 | 0.704 | 0.789 | 0.928 | 0.800 | 0.508 |
| `mamba4` | `pred_preddepth` | 0.745 | 0.889 | 0.811 | 0.980 | 0.595 | 0.770 |
| `mamba4` | `+ poisson 0.05` | 0.903 | 0.719 | 0.801 | 0.948 | 0.791 | 0.512 |

The two architectures land within about 0.01 F1 of each other on every arm. On a single
chromosome with 729 reference calls that is not a basis for preferring one over the other.

Read three things. **Annotated-CDS F1** should be high; this is the class the annotation makes
trustworthy, and 0.980 is what a working run looks like. **Non-canonical classes are the hard
case**, and the uncalibrated arm over-calls them: `dORF` precision is 0.222, about four false
calls per true one. **The Poisson arm is the dial**, and comparing the two rows shows what it
buys: novel precision 0.576 to 0.800, for novel recall 0.779 down to 0.508. Report both arms.
Quoting either alone gives a number better than the method is.

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
