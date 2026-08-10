# FIGURE D12 -- model-selected ORFs vs CPAT / CPC2 coding potential

**What it shows:** the comparison a reviewer will ask for. Every other figure scores the model
against naive enumeration, which only establishes "better than listing every ORF". CPAT and CPC2 are
the standard sequence-only coding-potential selectors and they shrink the search space by the same
order of magnitude as the model, from sequence composition alone. Closes locked decision **D3**.

- **(a) raw novel peptides**, 1% class-specific FDR. `*` marks the best of the four SELECTION RULES.
- **(b) discovery efficiency**, novel peptides per 1,000 database sequences -- normalises out the
  fact that the arms build databases of very different sizes.
- **(c) canonical cost**, change in GENCODE peptides vs the gencode-only baseline (symlog).

---

## `make_cpat_cpc2.py` -> `D12_cpat_cpc2.{png,pdf}` + `D12_values.json`

### Data inputs

    proteogenomics/data/<DS>_pilot/pgx_mamba4/report/table_all.json     DS in {A549, HBL1, SUDHL4, DoHH2}

Fields read per arm: `novel_peptides`, `db_novel_seqs`, `d_gencode_peptides`. Released **mamba4**,
frozen search parameters, pgx RiboCode-called databases. Datasets or arms whose `table_all.json` is
absent are **printed as MISSING** and recorded in the JSON, never treated as zero.

Upstream of those tables: `pgx/coding_potential.py` builds the `cpat` and `cpc2` arms by enumerating
from **exactly the same candidate pool** as `null_atg` and then applying each tool's own classifier.
That equality is enforced by `tests/test_invariants.py::test_coding_potential_pool_matches_null_arm`,
and it is what makes this a controlled comparison of selection rules rather than of pipelines.

### Values (2026-08-07)

| dataset | substrate | null_atg | CPAT | CPC2 | model theta=1 | model Poisson |
|---|---|--:|--:|--:|--:|--:|
| A549 | tryptic | 10 | **25** | 23 | 9 | 11 |
| HBL-1 | HLA-I | 20 | 3 | 10 | 8 | **14** |
| SU-DHL-4 | HLA-I | 28 | 7 | 14 | **22** | 10 |
| DoHH2 | HLA-I | 13 | 4 | 8 | **32** | 14 |

Discovery efficiency (novel peptides per 1,000 DB sequences), the panel (b) quantity:

| dataset | null_atg | CPAT | CPC2 | model theta=1 | model Poisson |
|---|--:|--:|--:|--:|--:|
| A549 | 0.03 | 1.17 | 1.11 | 0.87 | **4.58** |
| HBL-1 | 0.07 | 0.19 | 0.66 | 0.77 | **7.62** |
| SU-DHL-4 | 0.11 | 0.85 | 1.67 | 2.13 | **6.17** |
| DoHH2 | 0.05 | 0.53 | 1.04 | 3.90 | **7.16** |

---

## The reading, including the part that does not favour the model

**Between the four selection rules, the result splits by assay.** CPAT and CPC2 win the tryptic
whole proteome (25 / 23 vs the model's 11); the model wins all three HLA-I immunopeptidomes (14 vs
10, 22 vs 14, 32 vs 8 -- a 4x margin on DoHH2). This reproduces, on an independent pipeline and
independent data, the pattern recorded for the human MS work in the biotype-probe project (tryptic:
FM ~ CPAT; HLA: FM > CPAT).

The mechanism is visible in the database sizes. CPAT and CPC2 score a *transcript's* coding potential
from sequence composition, so they keep long, codon-biased, ORF-like sequences -- exactly what a
tryptic digest samples well and exactly the least novel population available. The model scores
*per-nucleotide translation* from the sample's own RNA-seq, so it keeps short, non-canonical,
cell-type-specific ORFs, which is what HLA-I presentation samples.

**On raw count the unselected `null_atg` arm beats every selection rule on HBL-1 (20) and SU-DHL-4
(28).** Keeping every candidate does find more peptides. It does so at 25x to 30x the database size,
and on A549 at a cost of 2,894 canonical peptides. That is why panel (b) exists and why the `*` in
panel (a) ranks only the selection rules: `null_atg` is the reference they are trying to beat per
unit of database, not a competitor on raw count. On efficiency the model's Poisson arm leads all
four datasets by 4x to 100x.

Panel (c) shows the cost side is not what separates the selective arms: `model_poisson` costs exactly
zero canonical peptides on all four datasets, CPAT/CPC2 cost -100/-145 on A549 and near zero
elsewhere, and `model_standard` costs -191 on A549 but GAINS +63 on DoHH2. `null_atg` is the
expensive one (-2,894 on A549).

**Caveat, stated plainly:** single-digit to low-double-digit peptide counts throughout. The direction
is consistent across three independent immunopeptidomes, but no individual dataset carries the claim.
That is the motivation for decision D2 (expand the panel to 6).

### Figure-construction notes

- Panel (a) sets `ylim` to 1.45x the maximum before drawing, or the DoHH2 asterisk and the legend
  both clip.
- Panel (c) uses `symlog(linthresh=10)`. On a linear axis A549's -2,894 is ~40x every other bar and
  flattens the panel to a single spike, making the selective arms indistinguishable from zero.

### Rebuild

    cd figures/D12_cpat_cpc2 && \
      /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_cpat_cpc2.py

`--model attn` switches to the other released model once its CPAT/CPC2 arms have been searched (only
`mamba4` has them as of 2026-08-07).
