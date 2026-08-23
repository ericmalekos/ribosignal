# Caption veracity audit: `ribosignal_poster.pdf` (2026-08-16)

Checked every quantitative claim in the poster text against the artifacts on the cluster.
**Twelve claims verified exactly. Four problems, one of which I would fix before printing.**

The file that arrived was `ribosignal_poster.pdf`, not `poster_preview.svg` (md5
`2e057cd292e1856d6c607f515397a454`, 17.8 MB, 1 page). Text extracted cleanly, so this audit
covers the captions as written.

---

## 0. Check the orientation before anything else

The page is **2592 x 3456 pt = 36 in wide x 48 in tall, PORTRAIT**.

`docs/POSTER_LAYOUT.md` and the tracking task both specify **48 x 36 LANDSCAPE**. If the venue
specified landscape, the whole poster is rotated and everything below is moot until that is
settled. I do not know which is right, so flagging rather than assuming.

---

## 1. Verified correct

Every one of these reproduces exactly. No rounding disputes.

| Claim | Source of truth | Measured |
|---|---|---|
| F1 **0.831 to 0.911**, four human held-outs | `P5_values_attn.json`, `pred_preddepth` | THP-1 0.8307 (min) to hepatocyte 0.9106 (max) |
| Recall **0.996** annotated vs **0.098** dORF | `P5_by_class_attn.tsv` | pooled tp/n_ref = 0.9965 / 0.0984 |
| **5,071,106** parameters | summed from `best.pt` state_dict | 5,071,106 |
| **4,093 nt** receptive field | `1 + 2*(k-1)*sum(dilations)` = `1 + 4*1023` | 4,093 |
| CCL2 **100-codon** ORF, **99.0%** of P-sites, **85.2%** frame 0 | `human_ruizorera_attn_ranked.tsv` | 300 nt = 100 codons; 1347/1360 = 0.9904; f0 = 0.8515 |
| **69** Ribo-seq libraries | `dataset_census.tsv` | 6+6+32+6+5+3+11 = 69 |
| GSE182371 as the Ribo-seq accession | census `accession_rna` = GSE182372 | correct, and correctly distinguished |
| GTF2I ORF "runs to 3,323 nt" | per-nt TSV | called ORF 451..3324 |
| ... "window ends at 1,819, where the observed signal stops" | per-nt TSV | **1** P-site of 13,258 lies beyond 1,819 |
| Wang = **shallowest** of the three experiments | pack `target_counts.npy` | 74.5M < 106.1M (janich) < 143.2M (gse243134) |
| **18** mzML per population, **30 aa** ORF floor | `P12_values.json` | `mzml_per_population` 18, `min_aa` 30 |
| AUG null loses "**about a thousand** peptides per population" | `P11_values.json` | median dtotal **-1,129.5**, range -932 to -1,452 |

Three of these deserve a specific note:

**The recall aggregation is pooled, and that is the right call.** 0.996 / 0.098 matches
`sum(tp)/sum(n_ref)` across the four datasets. A mean-of-four would give 0.996 / 0.100. Pooled is
defensible and matches to the digit, so nothing to change -- just be ready to say which if asked.

**"Wang is the shallowest" does NOT contradict the RNA note I sent earlier.** Wang is the
shallowest arm on **Ribo-seq** (2 libraries, 74.5M P-sites) and the *deepest* on **RNA**
(6.88e9 coverage). Both captions are correct because they are on different axes. Worth knowing so
the two messages do not look like they disagree.

**The P-sites / RNA-coverage separate-axis treatment in panel 2 is exactly right.** Panel 2 already
says the two are not comparable. Hold that line: `n_rna_libraries` is also not a depth proxy --
on mouse liver, 19 RNA libraries buys roughly *half* the coverage of 2.

---

## 2. Four problems

### 2.1 Takeaway 3 breaks the bound that panel 5 carefully sets -- fix this one

> "In a comparison of ORF calling in mouse liver samples, the model had better recall than a real,
> high quality riboseq experiment (Figure 5)"

Panel 5's own caption is scrupulous about this:

> "The model does not beat deep Ribo-seq: it loses to both deeper experiments on recall, and to all
> three on precision and therefore F1. The win is bounded to a shallow experiment, at no
> experimental cost."

The takeaway drops "shallow" and substitutes "high quality". The experiment in question is **2
libraries at 75.4% frame-0** -- the shallowest of the three by a factor of two. Calling it high
quality inverts the very bound the panel spends a sentence establishing, and it is the claim a
reviewer standing at the poster will push on hardest.

Suggested replacement:

> "Against the shallowest of three mouse-liver experiments, the model recovers more non-canonical
> ORFs at no experimental cost; it still loses to the two deeper ones (Figure 5)"

### 2.2 "compared in panel 4a" -- that panel is numbered 5

Panel 2's caption cross-references **panel 4a** for the mouse-liver comparison. On this poster the
mouse-liver panel is **5. MOUSE LIVER: CROSS-STUDY COMPARISON**; panel 4 is human held-out
validation. Looks like a leftover from the old well naming. A reader who follows it lands on the
wrong panel.

### 2.3 Panel 2 loses Hepatocytes from its own arithmetic

> "Seven human cell types from one study (GSE182371), 69 Ribo-seq libraries; brain, the study's
> ninth, was dropped before training on low periodicity."

Seven, plus brain, is eight. The study has nine in the census, and the missing one is
**Hepatocytes (5 libraries)** -- the held-out set that the whole of panel 4 is built on. The 69
libraries are **training-only**, sitting under a heading that reads "TRAINING AND VALIDATION
DATASETS".

Full census for the nine:

| role | cell types | Ribo libraries |
|---|--:|--:|
| train | ES, Fat, Fibroblast, HA_EC, HCAEC, HUVEC, VSMC | 69 |
| holdout | Hepatocytes | 5 |
| dropped | Brain | 5 |

Suggested: "Seven human cell types used for training (69 Ribo-seq libraries) plus held-out
hepatocytes, from one nine-tissue study (GSE182371); brain was dropped before training on low
periodicity."

### 2.4 The reason given for recall-not-F1 is not the real reason

> "class values in (B) are RECALL, not F1, because this project defines precision over all ORFs
> only."

The practice is right; the justification is not, and the real one is a better answer if challenged.

`score_human_orf_calls.prf(pred, ref, subset=(c,))` filters **both** sides by class. An ORF that
both callers found but labelled differently -- uORF vs Overlap_uORF, internal vs uORF -- therefore
falls out of every class bucket. That makes the per-class `n_pred` wrong, and precision and F1
divide by it. **Recall is unaffected** because its denominator is the reference's own class
assignment.

Suggested: "class values in (B) are RECALL, not F1: the scorer filters both predictions and
reference by class, so ORFs the two callers class differently drop out of the precision
denominator. Recall is unaffected, since its denominator is the reference's own class assignment."

This is also documented at the top of `scripts/export_p5_by_class.py`.

### 2.5 Minor

References 4 and 5 are still `[author] et al. [year]`.

---

## 3. Not verifiable from here

Nothing material. Everything quantitative in the caption text traced to an artifact. The only
claims I could not check are rendering choices (which x-limit a panel draws, which rows appear in
a figure), and in the one case where a rendering choice carried a factual assertion -- GTF2I's
"where the observed signal stops" -- the assertion holds.
