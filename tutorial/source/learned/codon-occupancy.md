# Has the model learned codon-level ribosome occupancy?

The model never receives codon identity, amino acid, or tRNA abundance as an input -- only one-hot
sequence, an ORF track, and RNA-seq coverage. So any codon structure in its predictions was learned
from sequence. This asks whether there is any, and whether it transfers.

## The same statistic, computed twice

Every `pred_profiles.npz` carries `obs_flat` and `pred_flat` over an identical transcript axis, so the
empirical and model codon statistics are the SAME computation applied to two arrays -- no second
pipeline, nothing that can differ between the two sides except the numbers.

Per transcript with an annotated CDS and >=100 observed CDS counts: take the per-nt profile over the
CDS, divide by that transcript's CDS mean (removing expression, so one highly translated gene cannot
dominate), walk the CDS in codons, and average genome-wide. Both P-site and A-site (P+3, the
conventional dwell proxy) are reported; quoting one without saying which is a common way codon-dwell
numbers become non-comparable between studies.

The signal gate is applied to OBSERVED counts only. `pred_flat` is a per-transcript probability
distribution summing to 1.0, so gating it on a count threshold silently rejects everything; gating on
observed also guarantees both sides average over the identical transcript set.

## The empirical ceiling comes first

Codon dwell is protocol-sensitive -- cycloheximide, digestion, depth -- so a model correlation of
0.26 means very different things depending on what two real experiments achieve. Observed-vs-observed
across independent datasets (same-dataset pairs excluded, since they share one `obs_flat` and
correlate at 1.000 by construction):

| comparison | n pairs | median r | range |
|---|--:|--:|---|
| **within species** | 31 | **0.944** | 0.344 to 0.956 |
| **across species** (mouse vs human) | 36 | **0.269** | -0.006 to 0.435 |

Two findings before the model is involved at all. Codon occupancy is **highly reproducible between
independent experiments of the same species**, so it is a real and stable measurement. And it is
**almost entirely species-specific** -- a human codon-dwell profile predicts a mouse one at r = 0.27,
consistent with different tRNA pools.

```{figure} /img/C12_codon_occupancy/C12_codon_occupancy.png
:width: 100%

**Left:** the empirical ceiling, no model involved -- codon occupancy is highly reproducible within a
species (0.944) and almost entirely species-specific (0.269). **Right:** the model against whichever
ceiling applies to it. Mouse transcripts were never seen and reach 98% of the across-species ceiling,
i.e. the transfer limit. The human number is higher in absolute terms but sits on transcripts ~99.7%
present in training, so the two are not a ranking.
```

## The model, against the right ceiling

| substrate | model vs observed | relevant ceiling | fraction of ceiling |
|---|--:|--:|--:|
| human (THP-1, CAR-T) | **0.548** | 0.944 (within species) | **58%** |
| mouse (Janich, GSE243134, Wang) | **0.263** | 0.269 (across species) | **98%** |

Per dump: human 0.500-0.599 (n=4), mouse 0.145-0.474 (n=9).

**On human the model recovers 58% of the reproducible codon signal** -- substantial structure, learned
from sequence alone, for a quantity it was never given.

**On mouse it sits at the cross-species ceiling.** The right comparison for a human-trained model
predicting mouse codon dwell is not the within-species ceiling but the across-species one, because
that is the ceiling on how much codon structure transfers between the species at all. At 0.263
against 0.269, the model is at the transfer limit: it predicts mouse codon occupancy about as well as
a real human experiment does. The constraint is the biology, not the model.

## These numbers do not depend on the alignment recipe

The mouse side was recomputed from scratch after the Ribo-seq alignments were rebuilt under the
canonical recipe (`docs/PIPELINE_POLICY.md`), which changed total P-sites by 3.5-4.9%. Every median
above came back identical to three decimals.

That is not a skipped rerun -- the two sets of tables are different files, differing by up to 0.0385
in individual codons. They correlate at **r = 0.99952**. The reason is in the method: occupancy is
normalized per transcript against that transcript's own CDS mean, so a roughly uniform change in
read depth divides out.

Worth stating explicitly, because it is the kind of thing a reader should not have to take on trust:
the codon result would be the same had the data been processed with either recipe. Contrast the
profile-correlation metrics in {doc}`memorisation`, which are sensitive to exactly this.

```{admonition} What this does NOT show
:class: caution
- n is small (4 human, 9 mouse dumps) and the across-species ceiling rests on few independent pairs.
- The human datasets here are new (THP-1, CAR-T) but their TRANSCRIPTS largely overlap the training
  universe, so the human number carries the memorisation component quantified in
  {doc}`memorisation`. The clean statement is the mouse one.
- **The CAR-T half of the human number is built on a pack with ~29% PCR duplicates.** That pipeline
  ran `umi_tools dedup` on the genome BAM but built the pack from the undeduplicated transcriptome
  BAM (`docs/PIPELINE_POLICY.md`). Duplicates pile at specific positions rather than spreading
  evenly, so they distort profile shape. The human 0.548 is provisional pending a rebuild; the mouse
  number, which is the one the argument rests on, is unaffected.
- Recovering codon occupancy does not establish that the model represents tRNA abundance or codon
  optimality. It shows only that predicted profiles carry codon-correlated structure.
```
