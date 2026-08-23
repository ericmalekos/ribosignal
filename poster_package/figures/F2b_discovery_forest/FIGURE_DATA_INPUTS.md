# MAIN FIG 2b -- discovery-density forest plot

> **HISTORICAL, PRE-2026-08-15.** Every table and number from here until the "A549 REMOVED" section below was computed WITH A549 and is retained as the record of what it showed. It does NOT describe the current figure. The live numbers are in this panel's `*_values.json`; the removal and its cost are documented at the bottom of this file.

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


## Panel (b): TOTAL unique peptides vs gencode-only (added 2026-08-14)

Derived per arm from the same `table_all.json` / `table.json` the ratio panel reads:
`d_gencode_peptides + novel_peptides - gencode_arm.novel_peptides`. The baseline-novel term is
subtracted explicitly rather than assumed zero.

**Why it was added.** Panel (a) is a RATE RATIO. An arm can hold a 2,700x density advantage and
still hand back fewer peptides than not searching a novel database at all, and the ratio cannot say
so. Standing rule: all proteomics work reports absolute totals, including where they disfavour the
model.

| arm | above zero | values |
|---|--:|---|
| **model (Poisson)** | **9 / 10** | A549 +11 / **-35**, HBL-1 +14 / +12, SU-DHL-4 +10 / +13, DoHH2 +21 / +23, THP-1 +3 / +15 (mamba4 / attn) |
| model (theta=1) | 6 / 10 | A549 -182 / -185, HBL-1 -14 / -95, SU-DHL-4 +22 / +17, DoHH2 +95 / +35, THP-1 +5 / +8 |
| **naive AUG null** | **2 / 5** | A549 **-2,884**, HBL-1 -104, SU-DHL-4 +47, DoHH2 +17, THP-1 -5 |

The null does not depend on the model or the threshold arm, so it is drawn once per dataset as an
open X.

**The honest reading.** The Poisson-calibrated arm nets positive on 9 of 10 dataset x model
combinations; `theta=1` on 6 of 10; the naive null on 2 of 5. The single Poisson failure is
**A549 / attn = -35** (the tryptic proteome), and it must not be dropped from any caption -- the
matching mamba4 arm is +11, so the sign flips on model choice for that dataset.

**Do NOT merge with P11.** Across 12 mouse macrophage populations the model's database is above the
GENCODE-only baseline in only 5 of 12 (median -8). Five human datasets here, 12 mouse populations
there; the joint statement that survives both is that the model's database is roughly cost-neutral
to slightly positive while the naive nulls are reliably negative.

## ORF-length floor: this panel uses 7 AA, the ORF-call panels use 30 AA (declared 2026-08-14)

The ORF-call track filters at 90 nt = exactly 30 aa (`build_loader`, applied identically to model
and real Ribo-seq calls). Proteogenomics databases are built at `--min-aa 7` ("detectability"),
because HLA-I peptides are 8-11 aa. **So this panel spans ORFs the ORF-call panels exclude by rule**
-- 22% of the model database and 65% of `null_atg` are under 30 aa.

**IT MATTERS HERE.** Novel peptides at 1% class-specific FDR whose supporting ORFs are ALL under
30 aa -- the ones a 30-aa floor would delete -- across these 5 human datasets: **model Poisson 14 of
52 (26.9%)**, model theta=1 25 of 76 (32.9%), `null_atg` 22 of 72 (30.6%), `null_nc` 2 of 38 (5.3%),
but **CPAT 0 of 39 (0.0%) and CPC2 0 of 55 (0.0%)**. (The mouse macrophage sweep gives 0.0% for the
model -- the opposite answer, which is why it must not be quoted for these panels.)

**Against the null the effect is near-symmetric** (26.9% vs 30.6%), so the 7-aa floor does not
flatter the model there. **Against CPAT/CPC2 it is not symmetric**: their databases are 0-2% short
by construction, because composition-based coding-potential scoring selects long ORF-like sequences.
A 30-aa floor cuts the model ~30% and leaves them untouched.

**Consequence for D12's claim specifically: "the model wins all three HLA-I immunopeptidomes"
becomes "wins 2 of 3" under a 30-aa floor** -- HBL-1 flips (model 14 -> 8 vs CPC2 10 -> 10);
SU-DHL-4 (16 vs 14) and DoHH2 (15 vs 8) survive. Any caption comparing the model to a
coding-potential selector must state that the advantage depends on sub-30-aa ORFs.

Measured by `proteogenomics/scripts/orf_length_audit.py`; declared in `docs/PIPELINE_POLICY.md`.
Re-run it if the databases are rebuilt -- do not assume the 0% holds.

## A549 REMOVED 2026-08-15 -- this panel is now HLA-I ONLY

A549 was archived and removed from every figure
(`proteogenomics/data/_archive_a549_2026_08_15/README.md`). It was the project's **only tryptic
human dataset**, so what remains here is exclusively HLA-I immunopeptidome data.

**Never describe any number in this panel as a whole-proteome or tryptic result**, and do not
generalise it beyond HLA-I immunopeptidomics. The remaining datasets (HBL-1, SU-DHL-4, DoHH2, THP-1)
are all `nonspecific` / `num_enzyme_termini = 0`, where a 7-aa ORF floor is correct and permanent
because HLA-I peptides are 8-11 residues.

**What removal cost this figure specifically is recorded below.** Removing the only dataset on which
the model performed WORST is a real change to what these panels can claim, and it is documented
rather than absorbed silently.

**Cost of the A549 removal.** 20 points -> **16** (4 datasets x 2 models x 2 arms). The median
discovery-density ratio is unchanged at **93.6x**. The Poisson arm is now positive in **8/8**
rather than 9/10, because the two A549 points were the only failures: attn **-35** total peptides
and mamba4 **+11** on identical data. **The caveat that the sign flipped on model choice can no
longer be made from this figure** -- with A549 gone there is no dataset where the arms disagree.
A 8/8 result is cleaner than 9/10 and is also less informative; do not present it as a
strengthening of the evidence.
