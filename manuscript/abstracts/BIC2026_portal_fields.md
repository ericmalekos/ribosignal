# BIC 2026 portal fields (entered manually, NOT in the PDF)

## 1. List of Contributing Authors

Required delimiter format: `First and Last Name, Department Name; Institution Name`

```
Eric Malekos, <DEPARTMENT>; University of California, Santa Cruz
<Co-author>, <DEPARTMENT>; University of California, Santa Cruz
```

FILL IN: department, co-authors, and confirm the institution. Names inferred from the cluster
account and are a guess.

## 2. Abstract Summary Block (one paragraph, under 2,000 characters)

Ribosome profiling is the standard method for finding open reading frames that genome annotation
misses, including upstream ORFs and ORFs on long non-coding RNAs, but it is technically demanding
and exists for only a handful of cell types. RNA-seq is available almost everywhere. We built
RiboSignal, a dilated convolutional and transformer model of about 5 million parameters that
predicts the per-nucleotide ribosome P-site profile of a transcript from its sequence and RNA-seq
coverage alone, with no ribosome profiling experiment required. Feeding the predicted profile to a
standard ORF caller recovers annotated coding sequences at F1 0.90 on a human tissue held out of
training, and calls upstream and non-canonical ORFs at precision 0.77 and 0.59. A calibration step
that resamples the predicted profile as Poisson counts at a scaled sequencing depth is what makes
the non-canonical calls usable, roughly doubling their precision at fixed canonical recall. Applied
to mouse macrophages, a tissue and species never seen in training, the predicted ORFs form a
compact proteomic search database of 1,533 sequences that identifies 37 novel peptides by mass
spectrometry while costing only 36 canonical spectrum matches. Naively enumerating every possible
near-cognate ORF instead yields 2.5 million sequences and destroys 29,633 canonical spectrum
matches through false discovery rate inflation, showing that a principled, model selected ORF set
is not merely more convenient but necessary for the search to work at all. The same pattern holds
on a deep human proteome, where a 2,404 sequence predicted database found 11 novel peptides at zero
canonical cost against 11,278 lost by the 3.9 million sequence enumeration. The approach makes
non-canonical translation accessible in any tissue for which expression data already exist, and
runs on CPU without a ribosome profiling experiment.

## 3. Topic area

`AI and bioengineering`
(alternative: Synthetic, Systems, and Computational Biology)

## 4. PDF filename

`MALEKOS_Eric_AIandBioengineering_UCSC.pdf`

Pattern is LASTNAME_FirstName_BroadTopicArea_Institution.pdf

## Formatting checklist before export

- [ ] Exactly one page, 8.5 x 11, front only. Overflow is not processed.
- [ ] Margins 0.75 inches on all sides
- [ ] Body text Times New Roman 11 pt
- [ ] Title Times New Roman 12 pt bold
- [ ] Figure caption Times New Roman 10 pt
- [ ] At most ONE figure or table on the page
- [ ] No references (not permitted)
- [ ] No hyperlinks to supplemental content
- [ ] Export to PDF, then check the page count again

## Deadline

August 5, 2026.
