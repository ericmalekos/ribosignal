# What sequence does the model read?

Figure C10 computed input saliency for THREE transcripts and found the start codon ranked 3rd
for a uORF but 292nd for a canonical CDS. This repeats that gradient over ~1,800 ORFs per
substrate, for both released models, so the asymmetry either survives as a distribution or does
not.

Target is the predicted profile log-probability summed over an ORF's in-frame positions;
saliency is |d target / d one-hot| per nucleotide. The coverage channel's gradient is tracked
separately, which is what separates *reads the sequence* from *reads the RNA-seq*.

**Substrate matters.** Mouse GSE243134 is primary because none of its transcripts were ever in
training. Human Hepatocytes is retained as the deliberately CONTAMINATED arm ({doc}`memorisation`).
A pattern present on BOTH is a learned rule; one present only on human is a memorisation candidate.

The mouse arms below were recomputed on canonically-aligned data (`docs/PIPELINE_POLICY.md`); the
human arms are unchanged because those packs were not rebuilt. Class ordering is identical before and
after, which is what the method predicts: saliency is a gradient with respect to sequence, ORF track
and coverage, and the realignment moved none of the three.

## Per-class saliency, by model and substrate

### mamba4 on mouse (PRIMARY -- transcripts never seen)

| ORF class | n | median start-codon rank | in-ORF saliency ratio | sequence fraction |
|---|--:|--:|--:|--:|
| uORF | 265 | 35 | 9.59 | 0.635 |
| Overlap_uORF | 292 | 36 | 7.13 | 0.613 |
| internal | 285 | 74 | 6.07 | 0.596 |
| novel | 297 | 81 | 2.63 | 0.538 |
| dORF | 295 | 84 | 3.48 | 0.816 |
| Overlap_dORF | 55 | 99 | 3.49 | 0.787 |
| annotated | 299 | 160 | 1.28 | 0.602 |

### attn on mouse

| ORF class | n | median start-codon rank | in-ORF saliency ratio | sequence fraction |
|---|--:|--:|--:|--:|
| uORF | 265 | 53 | 6.20 | 0.619 |
| Overlap_uORF | 292 | 52 | 4.49 | 0.609 |
| novel | 297 | 114 | 2.09 | 0.472 |
| Overlap_dORF | 55 | 116 | 2.82 | 0.766 |
| annotated | 299 | 146 | 1.22 | 0.603 |
| dORF | 295 | 148 | 2.84 | 0.803 |
| internal | 285 | 157 | 3.81 | 0.573 |

### mamba4 on human Hepatocytes (contaminated arm)

| ORF class | n | median start-codon rank | in-ORF saliency ratio | sequence fraction |
|---|--:|--:|--:|--:|
| uORF | 265 | 37 | 9.77 | 0.538 |
| Overlap_uORF | 296 | 42 | 5.86 | 0.542 |
| novel | 295 | 58 | 2.94 | 0.424 |
| internal | 278 | 70 | 6.93 | 0.566 |
| dORF | 210 | 107 | 4.20 | 0.783 |
| Overlap_dORF | 66 | 88 | 3.72 | 0.691 |
| annotated | 300 | 122 | 1.20 | 0.543 |

### attn on human Hepatocytes (contaminated arm)

| ORF class | n | median start-codon rank | in-ORF saliency ratio | sequence fraction |
|---|--:|--:|--:|--:|
| uORF | 265 | 67 | 6.05 | 0.492 |
| Overlap_uORF | 296 | 62 | 3.78 | 0.478 |
| novel | 295 | 103 | 2.19 | 0.343 |
| internal | 278 | 136 | 4.51 | 0.480 |
| dORF | 210 | 108 | 3.12 | 0.730 |
| Overlap_dORF | 66 | 116 | 2.74 | 0.615 |
| annotated | 300 | 164 | 1.19 | 0.471 |

## What replicates

**The C10 asymmetry is real and it is a learned rule.** uORF start codons rank ~35 (mamba4) and
~53 (attn) on mouse; canonical CDS ranks ~160 and ~146. A 3-5x difference, same ordering on both
substrates and both models. Because it holds on mouse, where no transcript was ever seen, it is
not memorised sequence.

**The model reads different inputs for different ORF classes.** dORFs are strongly
sequence-driven (sequence fraction 0.77-0.82); novel ORFs are the most coverage-driven
(0.34-0.54); canonical CDS leans on coverage and has the WORST-ranked start codon of any class.
The model finds annotated ORFs largely from RNA-seq and non-canonical ones largely from sequence.

**Unseen transcripts push it toward sequence.** Sequence fraction is higher on mouse than human in
every class (annotated 0.603 vs 0.471 for attn; uORF 0.619 vs 0.492). With no memorised profile
available the model falls back on sequence -- mechanistically consistent with {doc}`memorisation`.

**mamba4 attends far more sharply than attn** at identical ORFs (internal start rank 74 vs 157;
uORF in-ORF ratio 9.59 vs 6.20), consistent with its better held-out Pearson and its 5/5 win on
cross-species CDS F1.

```{admonition} The internal-ORF puzzle
:class: important
`internal` ORFs score F1 **0.022** against a between-experiment ceiling of 0.539-0.620. Yet mamba4
ranks their start codons at 74 on mouse and attends 6.1x above transcript mean inside them. The
model LOCATES internal ORFs and still does not call them.

Two caveats before reading that as a calling problem. Internal ORFs are nested inside CDS, so a
high in-ORF ratio partly reflects CDS signal the region would carry anyway. And saliency measures
what the model attends to, not what it decides. The ORF-track channel ablation is the clean test,
and it determines whether in-silico mutagenesis around start codons is worth running at all: if
removing the start channel barely moves internal F1, the failure is not a missing sequence feature.
```

```{figure} /img/C11_channel_ablation/C11_channel_ablation.png
:width: 100%

ORF-track channel ablation, 6 arms x 2 models. Red = the model calls that class **better without**
the channel. The `annotated` column is flat (~0.000) because canonical CDS calling does not lean on
this track at all. The `internal` column is red for every frame channel: those channels encode where
the annotated frame lies, which is precisely the prior that argues against an ORF nested inside a
known CDS in a different frame. `ch3` (start propensity) is the mirror image -- worthless on
annotated CDS, decisive for uORFs.
```

## The channel ablation resolves it

Each of the five ORF-track channels was zeroed in turn and prediction re-run, plus an `f012` arm
zeroing all three frame-occupancy channels together. Twelve arms (6 x 2 models), scored by
`scripts/score_channel_ablation.py` against the same observed reference. Positive delta = the model
calls that class **better without** the channel.

| channel | annotated | uORF | internal |
|---|--:|--:|--:|
| ch0 frame-0 occupancy | -0.000 / -0.001 | -0.016 / -0.000 | **+0.053 / +0.039** |
| ch1 frame-1 occupancy | +0.000 / -0.000 | -0.018 / -0.014 | **+0.031 / +0.037** |
| ch2 frame-2 occupancy | -0.001 / -0.001 | -0.022 / -0.006 | **+0.050 / +0.009** |
| ch3 start propensity | +0.000 / +0.001 | **-0.069 / -0.073** | -0.000 / -0.018 |
| ch4 is_stop | -0.001 / -0.001 | -0.024 / +0.007 | +0.031 / +0.006 |
| f012 all frame channels | -0.002 / -0.002 | -0.048 / -0.037 | +0.036 / +0.025 |

*attn / mamba4, `pred_obsdepth`. Full table: `results/orf_channel_ablation_canon/ablation_scores.tsv`.*

**The answer to the puzzle is that the frame channels suppress internal ORFs.** Removing any one of
them roughly triples internal F1 (0.008 to 0.061 for attn; 0.025 to 0.065 for mamba4). Those channels
encode where the *annotated* reading frame lies, which is exactly the prior that argues against a
plausible ORF sitting inside a known CDS in a different frame. The model locates internal ORFs, as
the saliency showed, and is then told by its own input not to call them.

This is not a bug to fix so much as a design consequence to state. The same prior is what makes
annotated CDS calling near-perfect (F1 0.986-0.990, and removing frame channels costs `ALL` -0.024).
Internal ORFs are the price.

**ch3 is the mirror image.** The start-propensity channel is worth nothing on annotated CDS
(delta +0.000) and everything on uORFs (-0.069 / -0.073) -- the one class whose defining evidence is
a start codon upstream of the annotated one. It is also the only channel whose removal *hurts*
internal ORFs, and it hurts them far less than the frame channels help.

```{admonition} Why in-silico mutagenesis was cancelled
:class: note
The ISM plan was to mutate sequence around start codons and watch internal-ORF F1 respond. The
ablation makes that uninformative in advance: removing the start channel entirely moves internal F1
by -0.0002 (attn). A feature worth nothing when deleted will not become informative when perturbed
more gently. The bottleneck is the frame prior, which no amount of sequence mutagenesis reaches.

Preferring the cheap decisive experiment over the expensive suggestive one is the general point;
here it saved a GPU sweep that could only have confirmed a null.
```

```{admonition} Does this mean the frame channels should be removed?
:class: warning
No -- and the heatmap on its own could be read that way.

The `internal` gain is **noise, not recovered ORFs**. In the `f012` arm internal precision is 0.035
(attn) / 0.041 (mamba4), so ~96% of the internal calls the model starts emitting are wrong; F1 rises
only because recall climbs from ~0 to ~0.06 on top of that. The frame channels are not hiding good
internal calls -- they suppress a class the model calls badly either way.

The classes that LOSE are also bigger, and are the ones the application needs: uORF -0.048/-0.037 and
novel -0.038/-0.043 over 697 and 663 reference ORFs, against 207 internal ORFs gaining. Aggregate F1
falls 0.024 on both models.

And the deeper limit: every arm here zeroes a channel on a model **trained with it**, so this measures
what the trained model RELIES on, not what an architecture without the channel would achieve. A model
trained without frame channels could redistribute onto sequence and never acquire the suppressing
prior at all. That question is open; the clean test is a training arm on a frame-zeroed ORF track.
```

All twelve arms were re-run on canonically-aligned data. Sign agreement with the original off-recipe
run is 57/60 across five ORF classes, two models and six arms; the three disagreements all sit within
0.005 of zero, i.e. noise about no-effect rather than reversals.
