# Message to the poster session: RNA columns + mouse-liver RNA provenance

Rsync when ready. Package md5 `663f0b0932a4f8c795c90065e9c8c791`.

```bash
rsync -av --delete --checksum \
  emalekos@mustard:/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/poster_package/ \
  ~/RNA_Futures/poster_package/
```

---

## 1. The three columns are in, and yes, it is a raw total

`dataset_census.{tsv,json}` now carries `n_rna_libraries`, `rna_coverage_total`,
`rna_source` and `accession_rna` for all 16 rows.

**You can label `rna_coverage_total` as depth.** I checked rather than assumed: it is
exactly `coverage.npy.sum()` over int32 per-nt counts. No normalisation, no
depth-correction. The pack's `global_mean_coverage` is exactly this value divided by
`sum_L`, i.e. the per-nt mean of the same quantity.

**One caveat that does affect the axis.** It is summed per-nt read depth over transcript
positions, whereas `psites` is a count of footprints. Different units. Put them on
separate axes and do not scale them against each other.

Both definitions are written into `dataset_census.json` under `quantities`, so the
labelling travels with the data instead of living in this message.

## 2. The trap you were right to worry about

`coverage_norm.json` has an `n_rna_samples` field that looks exactly like the column you
asked for. **It reports 0 for every `_canon` pack**, including Fibroblast, which holds
1.2 trillion coverage. Those packs reused their union twin's RNA verbatim, which is
correct because the final recipe does not touch RNA, but the rebuild never carried the
count forward. Reading that field directly would have printed "0 RNA libraries" across
the entire training set. Same failure mode as `psites_at_cds`.

Counts now resolve in order: the provenance file list, then `n_rna_samples` when it is
non-zero, then the union twin accepted only when `coverage_total` matches exactly.
Every row states which of those applied in `rna_source`.

**Unresolved rows are empty, never 0.** That is why `HBL-1` is blank: it has no pack in
this project, so there is nothing to read. Please keep it blank rather than rendering a
zero bar.

15 of 16 rows are filled.

Two rows behave in ways worth knowing before plotting:

- **GSE39561 shows GSE208041's exact RNA totals.** Correct, not a bug. It has no RNA of
  its own and borrows them. The `note` field says so.
- **B721's RNA is its entire model input.** Its pack is coverage-only by design. Its
  census row has `pack=None` because its *depth* comes from the external Ouspenskaia
  reference, so the RNA lookup had to be pointed at the pack explicitly, or the one
  dataset where RNA is everything would have shipped blank.

`accession_rna` is **GSE182372** for all nine Chothani tissues, which is a different
series from the Ribo-seq **GSE182371**.

## 3. Your question: the RNA behind the mouse-liver model call

**GSE243134.** 19 RNA-seq libraries, `SRR26055477` through `SRR26055*`, pooled into
`data/liver3x3/pool_gse243134_canon/coverage.npy`.

The scored cell is fixed, not an argmax: `attn` / `rna_input=gse243134` /
`arm=pred_preddepth`. Because `pred_preddepth` is the standalone arm, the model sees no
Ribo-seq for the query sample at all, so this RNA is the *only* experimental input
behind those recall numbers. The `ribo-*` half of the results directory name is the
scoring reference, not a model input.

**Do not use library count as a depth proxy on this panel.** The two run in opposite
directions here. All three arms sit on the same 22,974-tx universe, so these are
directly comparable:

| arm | RNA libraries | RNA coverage total | mean/nt | tx > 0 | top-100 share |
|---|--:|--:|--:|--:|--:|
| gse243134 | 19 | 3,242,315,370 | 57.54 | 22,503 | 31.1% |
| janich | 7 | 6,853,032,580 | 121.62 | 22,569 | 30.7% |
| wang | 2 | 6,884,416,652 | 122.17 | 22,605 | 39.6% |

19 libraries buys the **least** coverage, by roughly half.

There is a defensible point in this for the poster, if you want it: the model cell that
wins is running on the **shallowest** RNA of the three arms. The result is not bought
with deeper input.

## 4. A correction to E3's stated mechanism

E3's rationale used to say Wang RNA collapses uORF recall because "Wang is the
shallowest RNA arm (2 samples), and uORFs sit in 5'UTRs where coverage is thinnest".
The sample count is right and **the depth claim is backwards**. Wang RNA is the deepest
arm by total and by mean/nt, and covers the most transcripts.

What actually separates it is evenness: its top 100 transcripts hold 39.6% of all
coverage against about 31% for the other two, out of only 2 libraries. So the uORF
collapse is a concentration or replicate effect rather than a depth one.

I have corrected the comment in `make_shallow_experiment_win.py`. **No numbers change
and no figure needs regenerating.** But if that sentence made it into a caption or into
speaker notes, it needs replacing. The honest version is that the mechanism is not
established by these numbers, so state it as unexplained rather than swapping in a
second guess.
