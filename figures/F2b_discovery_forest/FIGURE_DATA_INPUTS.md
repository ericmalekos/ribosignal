# MAIN FIG 2b -- discovery-density forest plot

**What it shows:** the model's search-space advantage is a REPRODUCIBLE PATTERN, not a one-off.
FIGURES_PLAN names Fig 2's statistical power as its key weakness (the case originally rested on
HBL-1 alone) and specifies a forest plot of the model/null discovery-rate ratio as the fix. This is
that panel: **5 datasets x 2 released models x 2 threshold arms = 20 points, every one above 1.0.**

---

## `make_f2b_forest.py` -> `F2b_discovery_forest.{png,pdf}` + `F2b_values.json`

### Data inputs

    proteogenomics/data/<DS>_pilot/pgx_<MODEL>/report/table_all.json   (falls back to table.json)
    DS    in {A549, HBL1, SUDHL4, DoHH2, THP1}
    MODEL in {mamba4, attn}

Fields read per arm: `novel_peptides`, `db_novel_seqs`. All from the released models under frozen
search parameters with pgx RiboCode-called databases. Missing dataset/model/arm combinations are
**printed and recorded in the JSON**, never imputed.

### The metric is a RATE RATIO, and that is the whole argument

Discovery density = novel peptides per 1,000 database sequences, at 1% class-specific FDR. The
plotted value is `density(model arm) / density(null_atg)` within the same dataset and model.

Raw novel-peptide count is the wrong axis here and using it would flatter the null: a bigger database
absorbs more spectra, and `null_atg` does win on raw count in 2 of 5 datasets (HBL-1 20 vs 14,
SU-DHL-4 28 vs 22) while carrying **25x to 180x more sequences** and paying a canonical-ID cost the
model arms do not. Density asks the question that actually matters for search-space selection: per
sequence you commit to searching, how much do you find?

### Values (2026-08-07)

| dataset | substrate | mamba4 Poisson | attn Poisson | mamba4 theta=1 | attn theta=1 |
|---|---|--:|--:|--:|--:|
| A549 | tryptic | 168.8x | 112.5x | 32.2x | 27.6x |
| HBL-1 | HLA-I | 105.4x | 91.3x | 10.7x | 18.2x |
| SU-DHL-4 | HLA-I | 57.2x | 40.9x | 19.7x | 22.3x |
| DoHH2 | HLA-I | 133.6x | 166.2x | 72.8x | 95.9x |
| THP-1 | HLA-I | 529.3x | **2744.1x** | 179.4x | 259.2x |

**n = 20, all > 1.0. Median 93.6x, range 10.7x to 2,744x.** The Poisson arm leads the theta=1 arm in
every dataset, which is the calibration result restated on a fifth axis.

### Honest reading, to carry into the caption

- **THP-1's 2,744x is real but should not anchor the claim.** Its null found a single novel peptide
  from 301,855 sequences, so the denominator is one count away from being undefined. THP-1 also has
  the strictest FDR cutoff in the panel (23.5 vs 17.4-18.2) because a fast-scanning Orbitrap Fusion
  produced many marginal MS2. The median (93.6x) is the number to quote.
- **No error bars, deliberately.** These are single-search point estimates, not replicated
  measurements. The honest uncertainty statement is the peptide count printed beside each Poisson
  point (7 to 16 across datasets). A Poisson CI on those counts would imply a replication structure
  that does not exist.
- **The two released models are interleaved throughout**, consistent with every other head-to-head in
  the project. THP-1 is the one dataset where they diverge sharply (2,744x vs 529x); at 15 vs 3
  peptides that is noise, and the results.md entry says so.

### Figure-construction notes

- Log x-axis: the ratios span more than two orders of magnitude.
- `xlim` lower bound is pinned at **0.6** so the x=1 reference line stays on canvas. An earlier
  version scaled the lower bound off the data (3.6), which hid the line while keeping its
  "no advantage" label, implying the line was somewhere off-panel.
- Legend sits upper-left, where the panel is empty because every point exceeds 10x. A lower-right
  legend occluded the THP-1 row and truncated its annotation.

### Rebuild

    cd figures/F2b_discovery_forest && \
      /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_f2b_forest.py

Deterministic; re-reads every report JSON, so it picks up new datasets automatically once their
reports exist.
