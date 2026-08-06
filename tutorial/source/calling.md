# Calling ORFs from a predicted profile

The model emits a per-nucleotide P-site distribution and a predicted total. Turning that into ORF calls uses
the same caller the training data came from -- RiboCode -- with the predicted density substituted for the
experimental one. This page covers the substitution, the calibration dial, and the fact that the calibrated
arm is **stochastic**, which changes how its output must be used.

## The drop-in

`scripts/ribocode_dropin.py` bypasses BAM reading entirely. RiboCode's density object `tpsites_sum` is just a
`{transcript_id: np.ndarray}` dict, so the predicted profile is injected directly and passed to
`detectORF.main` with the **same prebuilt annotation and default parameters** the experimental runs used.
Nothing about the caller changes; only the numbers it reads.

```{list-table}
:header-rows: 1
:widths: 24 76

* - Variant
  - Density
* - `real`
  - observed pooled P-site counts (the reference)
* - `pred_obsdepth`
  - predicted **shape** x real per-transcript total -- isolates profile shape
* - `pred_preddepth`
  - predicted shape x **predicted** total -- fully standalone, no Ribo-seq at all
```

`pred_preddepth` is the deployment case: sequence + RNA-seq in, ORF calls out.

```{admonition} Why the density is integerised
:class: note
A smooth density has almost no exact zeros, so it slips through RiboCode's per-codon coverage and tie logic
far too easily and looks artificially periodic. Scaling to a realistic depth is also mandatory: RiboCode hard-
gates on per-transcript sum >= 5, per-ORF frame-0 sum >= 5, and >= 5 nonzero frame-0 codons.
```

## The stringency dial

$\theta$ scales the predicted effective depth:

```python
lam = p * (pred_total * theta)          # p is the predicted position distribution
```

$\theta$ = 0.05 simulates a Ribo-seq experiment at 5% of the predicted depth. Fewer counts means weaker
per-ORF evidence, so RiboCode's gates reject marginal ORFs: **lower $\theta$ is stricter**, and the stringency
is applied *through the data* rather than by thresholding a score.

Two arms are always reported ({doc}`metrics`):

- **standard** -- $\theta$ = 1, deterministic rounding. The naive standalone prediction, which over-calls.
- **poisson** -- $\theta$ = $\theta$*, density drawn as `Poisson(lam)`. The calibrated recipe.

## Choosing $\theta$*

Sweep $\theta$, run the full caller at each, and count how many annotated-CDS transcripts get a significant
`annotated` call. Annotated CDS are the ORFs we are most confident really are translated -- an internal
positive control present in every sample, requiring no labelled non-canonical truth.

$$\text{cds\_recall\_rel}(\theta) = \frac{\text{cds\_hits}(\theta)}{\max_\theta \text{cds\_hits}(\theta)}$$

**$\theta$\* = the smallest $\theta$ with `cds_recall_rel` >= target** (default 0.90). Smallest, because
tightening trades non-canonical yield for precision monotonically, so you want the most stringent point that
still retains 90% of the CDS the sample can support.

```{admonition} The denominator is relative, and it has to be
:class: important
Absolute recall against the full annotation **plateaus at 0.533** on BMDM (10,236 of 19,190 CDS-bearing
transcripts), so an absolute 0.90 target is unreachable *by construction* and the dial would always
saturate at $\theta$=1.

The cause is **annotation redundancy, not detection failure**. The universe carries **1.90** CDS-bearing
isoforms per gene; Ribo-seq signal concentrates on essentially one of them and the caller finds it. At the
ceiling those 10,236 transcripts span **9,777 of 10,080 genes**, so transcript-level recall of 53.3% is
**gene-level recall of 97.0%**, at 1.05 called isoforms per gene. The ceiling measures how many isoforms
the annotation lists, not what the sample translates.

The relative anchor asks the question the dial is actually for: at this stringency, what fraction of the
CDS this sample **can** support is retained. It normalises away that annotation artifact and is
self-calibrating across datasets of differing depth and isoform density.
```

BMDM, `mamba4` union model:

```{list-table}
:header-rows: 1
:widths: 12 22 16 14 14 22

* - $\theta$
  - cds_recall_rel
  - calls
  - uORF
  - dORF
  - lncRNA_orf
* - 0.02
  - 0.8404
  - 8,999
  - 178
  - 28
  - 183
* - **0.05**
  - **0.9504**
  - 10,786
  - 509
  - 58
  - 455
* - 0.10
  - 0.9863
  - 11,996
  - 918
  - 107
  - 822
* - 1.00
  - 1.0000
  - 17,436
  - 3,445
  - 822
  - 2,793
```

A guard checks the top two $\theta$ agree within 2% before trusting the ceiling; otherwise the grid needs
extending upward and `pgx.calibrate` says so.

```{admonition} "Optimal" means most-stringent-meeting-a-floor
:class: warning
$\theta$* is **not** optimal in a loss-minimising sense. There is no held-out non-canonical ground truth to
optimise against -- which is precisely why an anchor is used instead. The 0.90 floor is a chosen operating
point, and the proteogenomics results suggest it may be set too conservatively ({doc}`proteogenomics`).
```

## The Poisson arm is stochastic

`Poisson(lam)` is a **random draw**. The same model, the same $\theta$, and a different seed give a different
call set. This is not a detail to gloss: it changes what a single run's output means.

Five seeds at $\theta$* = 0.05 on BMDM:

```{list-table}
:header-rows: 1
:widths: 12 16 18 14 14 26

* - seed
  - total calls
  - annotated
  - uORF
  - dORF
  - lncRNA_orf
* - 0
  - 10,786
  - 9,730
  - 509
  - 58
  - 455
* - 1
  - 10,846
  - 9,736
  - 528
  - 60
  - 489
* - 2
  - 10,811
  - 9,735
  - 516
  - 55
  - 478
* - 3
  - 10,825
  - 9,746
  - 517
  - 62
  - 474
* - 4
  - 10,846
  - 9,745
  - 530
  - 65
  - 480
```

**Counts are stable; membership is not.**

```{list-table}
:header-rows: 1
:widths: 52 48

* - Quantity
  - Value
* - total calls, mean +- sd
  - 10,823 +- 23 (**CV 0.21%**)
* - pairwise Jaccard, 10 seed pairs
  - **0.822** (0.819 - 0.827)
* - called in **all 5** seeds
  - 8,869 = **66.9%** of the 13,250-call union
* - called in **exactly one** seed
  - 1,673 = **12.6%** of the union
```

So an aggregate count is reproducible to within a fraction of a percent, while roughly **18% of any single
seed's call set is draw-specific**. Two runs that differ only in seed agree on about four calls in five.

The instability is concentrated exactly where it matters least for calibration and most for discovery:

- `annotated` (canonical): 9,730 - 9,746, a spread of **0.16%**. The CDS anchor is rock solid, which is what
  makes $\theta$* itself stable.
- non-canonical: uORF spread ~4%, lncRNA_orf ~7%, dORF ~17% -- the classes the proteogenomic search is for.

```{admonition} Consequences you must not ignore
:class: warning
1. **Always set and record the seed.** `--seed` is a first-class parameter of the result, not a nuisance
   argument. An unrecorded seed makes the run unreproducible.
2. **Never compare two arms or two models across different seeds.** The ~18% membership churn will swamp
   real differences in the non-canonical tail.
3. **A single-seed non-canonical call list is not a stable object.** Treat any downstream artefact built
   from one draw -- including a search database -- as one sample from a distribution.
```

## Mitigation: consensus over draws

The fix mirrors the peptide-level replication logic in {doc}`proteogenomics`: run the calibrated arm at
several seeds and keep the ORFs that recur. `pgx.consensus` does this.

```bash
python -m pgx.consensus --calib-dir <dir>/calib --theta 0.05 \
    --profiles p.npz --species mouse --seeds 5 --min-seeds 3 --out <dir>/consensus
```

Requiring an ORF in >= 3 of 5 draws removes the 12.6% single-draw tail while retaining the 66.9% core, at a
cost of five caller runs (a few minutes each). Use it whenever the call list is an end product -- a search
database, a published table -- rather than an intermediate.

On BMDM the cleanup is almost exactly the size the churn predicts:

```{list-table}
:header-rows: 1
:widths: 34 16 16 16 18

* - Novel calls
  - total
  - uORF
  - dORF
  - lncRNA_orf
* - single draw (seed 0)
  - 1,022
  - 509
  - 58
  - 455
* - **consensus (>= 3 of 5)**
  - **833**
  - 387
  - 43
  - 403
* - removed as draw-specific
  - -189 (**18.5%**)
  - -122
  - -15
  - -52
```

Nearly a fifth of the single-draw novel calls do not survive replication -- and they would otherwise have
gone straight into a search database as though they were findings.

```{admonition} The deterministic arm has no such problem
:class: note
$\theta$ = 1 with rounding is exactly reproducible, which is one reason both arms are always reported. It
over-calls, but it over-calls the *same way* every time.
```
