# What the model learned

Benchmarks say how well the model scores. They do not say what it has actually seen, which parts of
the input it reads, or whether it has learned anything about translation beyond locating ORFs. These
three pages answer those separately, and each is built so it can come out negative.

## The three answers

**It has memorised something, and the amount is measurable.** On held-out Hepatocytes, transcripts
the model never saw during training reach profile r = 0.148, while seen transcripts *matched on depth
and length* reach 0.248 -- disjoint confidence intervals. Gradient exposure is worth about +0.10
profile correlation. The naive comparison against all seen transcripts (r = 0.616) is inflated almost
entirely by depth. See {doc}`memorisation`.

**It reads different inputs for different ORF classes.** Saliency ranks uORF start codons ~35th on
mouse and annotated CDS start codons ~160th: the model finds canonical CDS mostly from RNA-seq
coverage and non-canonical ORFs mostly from sequence. On unseen mouse transcripts the sequence
fraction rises in every class, which is what a model falling back on sequence in the absence of a
memorised profile should do. See {doc}`attribution`.

**It has learned codon-level structure it was never given.** Predicted codon occupancy correlates
with measured occupancy at r = 0.263 on mouse against a cross-species empirical ceiling of 0.269 --
**97% of what is achievable**. Nothing in training supplied codon identity or dwell time. See
{doc}`codon-occupancy`.

## One negative result, kept

The ORF-track channel ablation was run to decide whether in-silico mutagenesis around start codons
was worth a GPU sweep. It showed that removing the start channel entirely moves internal-ORF F1 by
-0.0002, so the sweep was cancelled before it ran. The same ablation explained the internal-ORF
puzzle as a side effect: the frame-occupancy channels actively *suppress* internal ORFs, roughly
tripling their F1 when removed. That is reported in {doc}`attribution` rather than filed away,
because it is a design consequence of the input representation and belongs in the record.

```{admonition} These results survived a change of alignment recipe
:class: note
Partway through this work the Ribo-seq alignments were found to deviate from
`docs/PIPELINE_POLICY.md` and were rebuilt: 40 alignments, 3 pools, 9 packs, 18 dumps, then every
analysis re-run. Total P-sites moved by 3.5-4.9%.

Expectations were written down before the comparison rather than after. Saliency was predicted to be
unchanged (it is a gradient with respect to sequence, ORF track and coverage, none of which moved)
and was. Ablation signs were predicted to hold and did, 57/60. Codon occupancy was flagged as
genuinely uncertain, because its empirical side rebuilds from shifted P-sites -- it came back
identical to three decimals, for a reason worth reading in {doc}`codon-occupancy`.

Old and new outputs are both kept, under `results/*_canon` and
`results/_archive_offrecipe_2026_08_13/`.
```

```{toctree}
:maxdepth: 2

memorisation
attribution
codon-occupancy
```
