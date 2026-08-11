# Sequence-attribution experiments: which sequence features drive the ORF signal?

Status: PLAN, not yet run. Written 2026-08-11.

## Why this is worth doing now

The per-class results say the model is at the between-experiment ceiling on annotated CDS
(0.988 vs 0.992-0.995) and far below it on everything else, most starkly `internal` (0.022 against a
ceiling of 0.539-0.620) and `dORF` (0.134 against 0.516-0.635). Two real experiments reproduce those
classes; the model does not attempt them. **That gap is the question these experiments exist to
answer**, and no amount of further benchmarking will answer it.

There is already one suggestive result, from figure C10, but on only three transcripts:

| transcript | class | start-codon saliency rank | in-ORF saliency / mean |
|---|---|--:|--:|
| POLR1D | uORF | 3 | 8.69x |
| LINC02693 | lncRNA ORF | 10 | 4.62x |
| RPL32 | CDS | **292** | 3.40x |

For the non-canonical ORFs the start codon is among the highest-attribution positions; for the
canonical CDS it ranks 292nd. That hints the model finds uORFs by reading start context and finds
canonical CDS some other way (coverage, length, frame continuity). With n=3 it is an anecdote.

Separately, the model contains a directly readable learned feature: `start_ctx_conv`
(`scripts/model.py:115`), a 4x10 nucleotide-by-position kernel that **rediscovered the Kozak motif
unsupervised at r=0.807** without ever being told it exists.

## The four tests

### 1. Aggregate saliency (cheap, ~1 h)

Generalise C10 from 3 transcripts to thousands. Target = predicted profile log-probability summed
over an ORF's in-frame positions; backprop to the one-hot input; aggregate |d target / d one-hot| by
position relative to start and stop, by frame, and by ORF class.

Answers: does the start-context dependence really differ between canonical and non-canonical ORFs?
Is there 3-nt periodic structure in the attribution itself (i.e. has the model learned frame, or is
it copying coverage)?

**Does NOT need reference calls** -- saliency can run over all candidate ORFs in the track, so the
scarce classes are not sample-size limited here.

### 2. ORF-track channel ablation (cheap, ~2 h)

The ORF track is `(sum_L, 5)`. Zero each channel in turn, re-dump, re-run the drop-in, and record the
F1 change per ORF class.

Answers: which ENGINEERED features carry which class. Specifically, whether `internal` fails for want
of a sequence signal or for want of a coverage/expression signal -- the single most decision-relevant
unknown, because it determines whether the fix is a feature or an architecture.

Needs reference calls, so it is sample-size limited on `internal` (see below).

### 3. In-silico mutagenesis (expensive, ~1 day, strongest evidence)

Mutate every position to each of the three alternative bases, measure the change in predicted
in-frame signal. Produces a per-nucleotide x per-base effect matrix.

Answers: lets a PWM be read directly OUT of the model rather than correlated against a prior. The
Kozak rediscovery becomes a measurement instead of an inference, and any second motif the model uses
becomes visible.

Tractable only if scoped: +/- 30 nt around start codons across a few thousand ORFs, not whole
transcripts. Whole-transcript ISM on 36,668 transcripts is not affordable.

### 4. Codon-level attribution (moderate)

Aggregate test 1 or 3 by codon identity.

Answers: has the model learned codon optimality or tRNA-abundance effects? Genuinely novel, because
nothing in the input encodes them -- a positive result would be a finding about what is learnable
from sequence plus coverage alone. Requires many CDS, which is the one thing never in short supply.

## Which datasets, and why

### Do not choose by performance

The instinct to run on the best-performing dataset is the wrong one here, and worth stating plainly:
the best-performing dataset is best because its ANNOTATED CDS calls are easy, and annotated CDS is
already at the ceiling. Attribution there would characterise the part that already works. The classes
that need explaining are the ones where performance is worst. Dataset choice should follow the
question, not the leaderboard.

### CORRECTION (2026-08-11): "held out" does not mean the transcripts are unseen

An earlier draft of this document recommended held-out Hepatocytes on the grounds that "neither
released model has trained on it". **That was wrong**, and the error matters more for attribution
than it does for benchmarking.

LOTO holds out a TISSUE, not transcripts. From `dataset.py::loto_split`: train and val are
chromosome-split against each other for honest early stopping, but `test_tx` is "the held-out
tissue's scorable tx (>= min_signal, not chrM), **all chromosomes**". The nine Chothani tissues share
one fixed 36,668-tx universe; per the LOTO header, "per-tissue expression enters via the coverage
input, not the transcript set". Measured:

| Hepatocytes test transcripts | 33,918 | |
|---|--:|--:|
| also seen in TRAIN gradients (some held-in tissue) | 27,007 | **79.6%** |
| also seen in VAL (early stopping) | 6,813 | 20.1% |
| never seen in any held-in tissue | **98** | **0.3%** |

For a given transcript the SEQUENCE and the ORF TRACK are byte-identical between its train and test
appearances; only coverage and target differ. So the model has seen the exact input sequence of
99.7% of its "held-out" test set.

**How much this inflates benchmarks: probably little.** Cross-species mouse shares NO transcripts
with training and still scores close:

| test set | transcript overlap | standalone F1 | annotated recall |
|---|---|--:|--:|
| Hepatocytes (LOTO) | 99.7% | 0.909-0.913 | 0.997 |
| Janich mouse liver | 0% | 0.893-0.898 | 0.980-0.998 |
| Wang mouse liver | 0% | 0.867 | 0.999 |

A model relying on memorised transcript-to-shape lookups should collapse where no transcript was ever
seen. It loses ~0.02-0.04 F1. Caveats: mouse also differs in species, tissue and library, so this is
not a clean isolation; and the model has no transcript-ID input, so it can only memorise sequence,
which ~5M parameters over ~100 Mnt makes hard but not impossible for high-signal transcripts.

**How much it matters for ATTRIBUTION: a lot.** An F1 that is robust to memorisation is not the same
as a saliency map that is. Attribution on a transcript whose profile the model may have partly
memorised can reveal what was memorised rather than what rule was learned, and aggregating over
33,918 such transcripts propagates that contamination into the summary.

### Revised substrate choice

**Primary: mouse liver GSE243134.** Promoted from secondary. Its transcripts were never seen in any
form -- different species, disjoint transcript IDs -- so attribution cannot be reading back a
memorised sequence. It is also the deepest liver set (150.5 M P-sites, 21 Ribo samples) with its own
RNA, and it carries the failing classes (`internal` n=223, `dORF` n=286) that motivate the work.

**Secondary: held-out Hepatocytes (human), used deliberately as the CONTAMINATED arm.** Not dropped:
comparing attribution between a substrate the model has seen and one it has not is itself the
memorisation test. Features that appear in both are learned rules; features that appear only in the
human holdout are candidates for memorised sequence.

**Control worth running first (cheap): the 98 never-seen transcripts.** Score and attribute the 98
Hepatocytes transcripts absent from every held-in tissue, against a size- and signal-matched sample of
seen transcripts. n=98 is small and will not support per-class claims, but it is the only clean
within-human, within-tissue isolation of the memorisation variable that exists, and it costs almost
nothing.

### Explicitly NOT used

- **GSE39561** -- no usable P-sites, and its predictions duplicate GSE208041's.
- **CAR-T** as a primary -- smallest universe (23,887 tx) and fewest non-canonical reference calls
  (internal 87). Useful as a third replication point, not as a main substrate.
- **Any Chothani training tissue** -- the released models trained on these. Attribution over training
  data confounds "learned feature" with "memorised example", which is the entire failure mode the
  holdout exists to avoid.

## The real constraint is class sample size, not dataset choice

Reference calls available per class:

| dataset | annotated | uORF | novel | dORF | internal |
|---|--:|--:|--:|--:|--:|
| mouse janich | 10,515 | 720 | 769 | 234 | **187** |
| human GSE208041 | 10,654 | 1,191 | 1,347 | - | **138** |
| human CAR-T | 7,399 | 579 | 297 | - | **87** |

`internal` is the limiting class everywhere, and it is the class whose failure most needs explaining.
Consequences:

- **Test 1 is unaffected** -- saliency needs candidate ORFs, not reference calls.
- **Test 2 is underpowered on `internal`** in any single dataset. Report confidence intervals, and
  pool across datasets before drawing a conclusion about that class specifically.
- Do not report a per-class F1 delta on n=87 without an interval. That is how a noise fluctuation
  becomes a mechanism in a writeup.

## Both models, not one

Run tests 1 and 2 on `attn` AND `mamba4`. A feature that matters to only one architecture is a fact
about that architecture, not about translation. The two mixers already trade periodicity against
precision differently, so agreement between them is meaningful evidence and disagreement is itself a
result.

## Suggested order

1. Test 1 on Hepatocytes, both models. Cheapest, and settles whether the C10 n=3 asymmetry is real.
2. Test 2 on Hepatocytes, both models. Determines whether `internal` is a sequence or a coverage
   problem, which decides what any fix would even look like.
3. Whichever of tests 3 and 4 the first two motivate. Do not pre-commit: if test 2 shows `internal`
   fails on coverage rather than sequence, ISM around start codons is the wrong next experiment.
4. Replicate the surviving findings on mouse GSE243134.
