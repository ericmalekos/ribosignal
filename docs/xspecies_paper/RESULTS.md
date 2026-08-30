# Results: zero-shot transfer of a per-nucleotide Ribo-seq model from primates to yeast

Eric Malekos. Compiled 2026-08-29. Methods: `docs/xspecies_paper/METHODS.md`.

A per-nucleotide ribosome-profiling signal model trained on human and mouse alone was
applied, without any parameter update, to seven genomes spanning primates, zebrafish,
*C. elegans* and yeast.

## Zero-shot ORF calling matches the in-distribution benchmark

The same ORF caller was run under identical statistics on the model's predicted profile
and on the observed profile, and the two call sets compared (Figure 2a). Using the fully
de novo variant, in which no observed data enters the prediction:

| species | precision | recall | F1 |
|---|---|---|---|
| yeast | 0.997 | 0.999 | 0.998 |
| *C. elegans* | 0.976 | 0.981 | 0.979 |
| human | 0.912 | 0.933 | 0.922 |
| gorilla | 0.920 | 0.923 | 0.921 |
| chimpanzee | 0.920 | 0.922 | 0.921 |
| macaque | 0.921 | 0.909 | 0.915 |
| zebrafish | 0.809 | 0.981 | 0.887 |

The cross-species human arm (0.912 / 0.933 / 0.922) reproduces an established
in-distribution figure of 0.903 / 0.928 / 0.915 for the same variant and filters on
held-out human tissue, and the three non-human primates fall within 0.006 F1 of it. This
is the study's principal control: the model performs on unseen genomes at the level it
performs on its training distribution. A transformer and a Mamba state-space architecture
were evaluated independently and agree within 0.003 F1 on every arm.

## The headline number is carried by annotated coding sequence

Annotated CDS is 65 to 100% of every call set and is the easy case, recovered at recall
0.993 to 1.000 in every species including the four absent from training (Figure 2b). The
discriminating classes are those where the model must place density at a position the
annotation does not mark:

| species | uORF recall (n) | ncORF recall (n) |
|---|---|---|
| human | 0.712 (4,392) | 0.768 (574) |
| chimpanzee | 0.706 (4,611) | 0.757 (222) |
| macaque | 0.665 (5,327) | 0.747 (352) |
| gorilla | 0.580 (2,549) | 0.770 (187) |
| zebrafish | 0.645 (124) | 1.000 (21) |
| *C. elegans* | 0.533 (195) | 1.000 (5) |

ncORFs are defined on transcript biotype rather than on the caller's ORF-type label, which
conflates lncRNA ORFs with novel ORFs on coding transcripts. Median ORF length largely
explains the difficulty ordering: CDS 1,110 to 1,497 nt, ncORF 165 to 255 nt, uORF 135 to
213 nt. The short classes are where the periodicity test is fragile and where the model,
not the annotation, supplies the evidence.

**uORF recall does not track phylogeny, and an earlier draft of this analysis wrongly
claimed it did.** Gorilla (0.580) falls below macaque (0.665) despite being the closer
outgroup to human, and depth does not rescue the ordering: chimpanzee reaches 0.706 on
64.4 M pooled frame-0 P-sites while macaque reaches 0.665 on 123.9 M. Human, the training
species, is highest; beyond that the ordering follows neither divergence nor depth, and no
trend is claimed (Figure 2c).

## Profiles are predicted well above chance but imperfectly

Median per-transcript correlation between predicted and observed per-nucleotide profiles
is 0.49 to 0.54 across primates, 0.33 in *C. elegans*, 0.27 in yeast and 0.14 in zebrafish
(Figure 3a). A position-shuffled null gives |r| < 0.002 everywhere, so the signal is
positional rather than an artifact of predicting overall magnitude.

Predicting a transcript's total output is a much easier task and is reported separately so
the two are never conflated: 0.58 to 0.86. The gap between the two columns is the distance
still to be closed. Zebrafish is weakest on every profile metric and is also the only arm
that over-calls ORFs (precision 0.809, recall 0.981, predicted call set 21% larger than
observed), both consistent with a count head predicting too much density there.

## A tie-handling artifact produced spurious negative rank correlations

The rank correlations reported alongside these Pearson values were negative for primates
and zebrafish, which would imply anti-correlated profiles. They do not.

The shared Spearman implementation assigns ordinal ranks without averaging ties. These
profiles are 74 to 95% zeros, so every tied zero receives a distinct rank in positional
order and the observed rank vector becomes a proxy for position rather than signal.
Because the model correctly predicts low density in the 3'UTR, and those positions are
both observed-zero and last in position, the vectors are driven into apparent
anti-correlation. The effect scales with the zero fraction, and correcting the tie handling
makes every species positive and consistent with the Pearson value (Figure 3c): human
-0.195 to 0.379, zebrafish -0.405 to 0.142. The shared function was flagged rather than
silently changed, because it feeds previously reported columns elsewhere; any rank
correlation on sparse data through it is affected.

## Retaining multimapping reads changes neither the call set nor the model's score

Re-aligning all 67 Ribo-seq runs with up to 25 alignments per read and calling ORFs
independently changes the call count by 1.000x to 1.013x in six of seven arms, with no
class preferentially inflated (Figure 4a). A 1.33x inflation with a 2.5-fold rise in novel
calls, previously seen in a mouse liver dataset, does not reproduce in any of these
genomes.

The apparent exception, chimpanzee at 0.884x, is a depth artifact: five of its eight runs
lost three-nucleotide periodicity under the permissive posture and were dropped by the
caller's own gate, so the two call sets rest on different data. The comparison now counts
contributing runs per posture and flags any mismatch. Scoring the model against the
permissive call set instead of the unique-only one moves F1 by at most 0.001 on every
comparable arm (Figure 4b), so the reported accuracy is not an artifact of alignment
posture.

## Two data-quality findings

**The RefSeq *C. elegans* annotation is defective.** NCBI's GTF for GCF_000002985.6 writes
each isoform's own name into the feature-type column where `CDS` belongs, on 204,530 of
204,542 CDS rows. The file still yields a correct genome index, transcriptome and
transcript universe, and all transcript counts reconcile, so nothing upstream fails. The
caller does not fail either: with no annotated CDS to compare against it typed all 17,745
worm ORFs as novel and none as annotated. After repair the same arm yields 16,594
annotated ORFs (93.5%), from an identical call set of the same 17,745 genomic loci with
none gained or lost. Detection was always correct; only the labelling was wrong. The
defect surfaced by comparing ORF-type composition across species, not from any assertion
in the pipeline.

***Drosophila* ribosome profiling could not be used.** Two independent datasets were
tested and neither shows usable three-nucleotide periodicity: 15.2 / 56.2 / 28.7% across
frames for an RNase T1 library and 21.5 / 19.5 / 59.0% for an RNase I library, against
93.9 / 2.6 / 3.6% for the yeast positive control (Figure 1c). The caller assigns no P-site
offset for any read length in either. The arm was dropped rather than forced through by
lowering the periodicity threshold.

## Limitations

Predicted and observed calls are not independent: the model is trained to produce periodic
profiles and the caller tests frame-0 dominance, so high CDS agreement partly reflects both
sides being keyed to the same annotation. The uORF and ncORF columns are the load-bearing
evidence.

Several arms are shallow. Gorilla contributes three runs and zebrafish four (after four
were discarded for weak periodicity). The ncORF class holds only 5 to 21 loci in
*C. elegans* and zebrafish, so those recall values of 1.000 rest on very few events, and
the zebrafish uORF class (n=124) with precision 0.084 should be read as an over-calling
signal rather than a stable estimate.

Each species contributes a single tissue or developmental context, so within-species
generalisation across cell types is untested. The *C. elegans* Ribo and RNA arms are
replicate timecourses matched by stage rather than by lysate.
