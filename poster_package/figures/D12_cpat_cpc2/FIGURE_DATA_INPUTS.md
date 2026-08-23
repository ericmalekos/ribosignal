# FIGURE D12 -- model-selected ORFs vs CPAT / CPC2 coding potential

> **HISTORICAL, PRE-2026-08-15.** Every table and number from here until the "A549 REMOVED" section below was computed WITH A549 and is retained as the record of what it showed. It does NOT describe the current figure. The live numbers are in this panel's `*_values.json`; the removal and its cost are documented at the bottom of this file.

**What it shows:** the comparison a reviewer will ask for. Every other figure scores the model
against naive enumeration, which only establishes "better than listing every ORF". CPAT and CPC2 are
the standard sequence-only coding-potential selectors and they shrink the search space by the same
order of magnitude as the model, from sequence composition alone. Closes locked decision **D3**.

- **(a) raw novel peptides**, 1% class-specific FDR. `*` marks the best of the four SELECTION RULES.
- **(b) discovery efficiency**, novel peptides per 1,000 database sequences -- normalises out the
  fact that the arms build databases of very different sizes.
- **(c) canonical cost**, change in GENCODE peptides vs the gencode-only baseline (symlog).
- **(d) TOTAL unique peptides** vs the gencode-only baseline (symlog) = (c) + novel gained. The
  number an experimentalist takes home. Added 2026-08-14 under the standing rule that all
  proteomics work reports absolute totals, including where they disfavour the model.

---

## `make_cpat_cpc2.py` -> `D12_cpat_cpc2.{png,pdf}` + `D12_values.json`

### Data inputs

    proteogenomics/data/<DS>_pilot/pgx_mamba4/report/table_all.json     DS in {A549, HBL1, SUDHL4, DoHH2}

Fields read per arm: `novel_peptides`, `db_novel_seqs`, `d_gencode_peptides`. Panel (d) is derived:
`d_gencode_peptides + novel_peptides - gencode_arm.novel_peptides`. The subtraction is explicit
rather than assumed-zero; verified identical to a direct
`(canon+novel)_arm - (canon+novel)_baseline` recomputation on the frozen reports. Released **mamba4**,
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


## Panel (d): TOTAL unique peptides vs gencode-only (2026-08-14)

| dataset | substrate | null_atg | CPAT | CPC2 | model theta=1 | model Poisson |
|---|---|--:|--:|--:|--:|--:|
| A549 | tryptic | **-2,884** | -75 | -122 | -182 | **+11** |
| HBL-1 | HLA-I | -104 | -20 | +9 | -14 | **+14** |
| SU-DHL-4 | HLA-I | **+47** | +7 | +14 | +22 | +10 |
| DoHH2 | HLA-I | +17 | +11 | +15 | **+95** | +21 |

**The claim this supports is "never negative", not "always best".** The model's Poisson arm is the
only arm positive on all four datasets. `null_atg` wins SU-DHL-4 outright (+47) and `model_standard`
wins DoHH2 (+95) -- both are stated in the generator docstring rather than omitted.

**Do NOT merge this with P11.** Across 12 mouse macrophage populations the model's database is above
the GENCODE-only baseline in only 5 of 12 (median -8). Four cell lines here, 12 populations there;
the honest joint statement is that the model's database is roughly cost-neutral and the naive nulls
are not.

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

**Cost of the A549 removal -- this is the largest of any panel.** The figure previously split by
assay: **CPAT/CPC2 beat the model on the tryptic whole proteome, the model beat them on the HLA-I
immunopeptidomes**. A549 was the tryptic side. With it removed, only the half the model WINS
survives.

**This panel can no longer support 'the model beats CPAT/CPC2' as a general claim.** It shows the
model beating them on HLA-I immunopeptidomes, a substrate where short ORFs are the biology and
CPAT/CPC2's databases are 0-2% short by construction. The honest caption is substrate-scoped.
The tryptic comparison, in which the model lost, is preserved in
`proteogenomics/data/_archive_a549_2026_08_15/`.
