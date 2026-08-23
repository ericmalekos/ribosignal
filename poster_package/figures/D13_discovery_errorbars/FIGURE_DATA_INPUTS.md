# FIGURE D13 -- the database-design tradeoff

> **HISTORICAL, PRE-2026-08-15.** Every table and number from here until the "A549 REMOVED" section below was computed WITH A549 and is retained as the record of what it showed. It does NOT describe the current figure. The live numbers are in this panel's `*_values.json`; the removal and its cost are documented at the bottom of this file.

**What it shows:** three panels, one story. Predicting WHICH ORFs to include keeps the search
database compact, which mitigates the FDR-recalibration penalty a naive enumeration pays.

- **(a) GAIN** -- discovery efficiency: novel peptides per 1,000 novel DB sequences, 1% class-specific FDR.
- **(b) COST** -- canonical displacement: change in GENCODE PSMs vs the gencode-only baseline (dPSM).
- **(c) DRIVER** -- database size: novel sequences added, i.e. added decoy load. This causes both (a) and (b).

---

## CURRENT: `make_d13_pgx.py` -> `D13_discovery_tradeoff_pgx.{png,pdf}` + `D13_pgx_values.json`

Built from the pgx pipeline on the RELEASED models. Four datasets x both shipped models.

### Data inputs

- `proteogenomics/data/<DATASET>_pilot/pgx_{mamba4,attn}/report/table.json`
  for DATASET in {HBL1, DoHH2, SUDHL4, A549}, written by `pgx.report` via `pgx.run --report-only`.

The script reads those JSONs directly, so the figure tracks whatever the pipeline last produced.

### Regenerate

```bash
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd $NEW/figures/D13_discovery_errorbars && $PY make_d13_pgx.py
```

### Design decisions worth knowing

- **Both released models are drawn** (solid = mamba4, hatched = attn). The MS application is
  model-agnostic -- 49 vs 48 novel peptides across the four datasets -- so drawing only the GPU-only
  primary would imply a GPU is required for the discovery result. It is not.
- **(a) and (c) are log-scaled.** On a linear axis the null arms (0.01-0.11 per 1k) are visually
  indistinguishable from zero next to the model (3-9 per 1k), which reads as missing data rather
  than as a 100-470x difference.
- **(b) plots dPSM SIGNED, and marks exact zeros with a diamond plus a corner note.** Zero canonical
  cost is the model's headline result on 6 of the 8 model arms; an absent bar would read as "not
  measured". Same 0-means-missing-vs-0-means-zero trap the pipeline hit elsewhere. (The note sits in
  the axes corner because an annotation placed near the bars ran behind them.)
- **Absent arms are drawn as gaps, never as zeros.** `null_nc` exists only for A549: under
  nonspecific digestion the ~3M-sequence near-cognate null generates 2.86e9 peptides against
  MSFragger's ~2e9 cap, so the immunopeptidomes have no such arm (see proteogenomics/methods.md).
  The script prints a `[warn]` listing anything it could not load.

### Headline numbers (model_poisson)

| dataset | mamba4 density /1k | attn density /1k | mamba4 dPSM | attn dPSM |
|---|--:|--:|--:|--:|
| HBL-1 | 7.621 | 6.597 | +0 | +0 |
| DoHH2 | 7.161 | 8.909 | +18 | +18 |
| SU-DHL-4 | 6.169 | 4.410 | +0 | +0 |
| A549 | 4.576 | 3.049 | +0 | -60 |

A549 is the clearest case because it is ~16x deeper than the immunopeptidomes: the model finds 11
novel peptides at zero canonical cost from 2,404 sequences, while the near-cognate null finds 38 but
destroys 11,278 canonical PSMs (12.6% of baseline) from 3.9M sequences.

### Caveats

- n = 4 datasets, single-digit-to-low-double-digit novel peptide counts per arm. The density
  separation is large and consistent; the absolute counts are small.
- Novel columns are 1% CLASS-specific FDR (novel targets vs `REV_nuORF|` decoys only); dPSM is 1%
  GLOBAL FDR against the gencode baseline, so a DB-size cost stays separable from a DB-quality cost.

---

## SUPERSEDED: `make_discovery_errorbars.py` -> `D13_discovery_errorbars.{png,pdf}` + `D13_values.json`

The original, retained for provenance. **Do not use for the manuscript.** Wrong in three independent
ways, all fixed above:

1. **Old model** -- `orf_v2_attn_onehot_holdout_Hepatocytes` with the heuristic-Kozak ORF track, five
   changes behind the released primary (architecture, union universe, Brain dropped, `--kozak none`,
   mm1 coverage).
2. **f0-threshold ORF selection** -- every candidate ORF with `pred_frame0 >= 0.5`; no RiboCode, no
   Poisson calibration, no significance test. Databases ~44x larger than the pgx equivalents
   (HBL-1: 80,608 vs 1,837 novel sequences).
3. **Unfrozen `calibrate_mass = 2`** -- MSFragger re-derived search parameters per DATABASE, so the
   arms were not scored under the same rules. On HBL-1 the null arm ran `use_topN_peaks=100` /
   `intensity_transform=0` against the model arm's 150 / 1 (score cutoff 16.742 vs 18.483).

It also covered only 3 datasets and one model. Its inputs (`results/raw_search_fdr.json`,
`results/mokapot_stochastic.json`) come from the retired f0 search tree.

Its own historical rationale, kept because the reasoning is still sound and applies to the current
figure too: the deterministic target-decoy FDR was chosen over the mokapot class-FDR because the
mokapot per-run count is stochastic (SU-DHL-4 model bimodal) and its confident-canonical count is
non-monotonic across DBs (mokapot re-learns per DB), so it could not cleanly show the
FDR-recalibration mechanism. The sparse-class rescoring instability is itself a real finding -- see
the stochasticity study and `results/mokapot_stabilized.json`.


## Panel (d): TOTAL unique peptides vs gencode-only (added 2026-08-14)

Derived per arm as `d_gencode_peptides + novel_peptides - gencode_arm.novel_peptides` from the same
`table.json` the other panels read. The subtraction is explicit rather than assumed-zero.

**Why it was added.** Panels (a) and (b) can both favour an arm that nonetheless leaves the
experimentalist with fewer peptides than not searching a novel database at all, and (b)'s `dPSM` is
a PSM count, not a peptide count. Standing rule: all proteomics work reports absolute totals,
including where they disfavour the model.

| model | arm | HBL-1 | DoHH2 | SU-DHL-4 | A549 |
|---|---|--:|--:|--:|--:|
| mamba4 | model (Poisson) | **+14** | **+21** | **+10** | **+11** |
| mamba4 | model (theta=1) | -14 | +95 | +22 | -182 |
| mamba4 | null AUG | -104 | +17 | +47 | **-2,884** |
| mamba4 | null near-cognate | -- | -- | -- | **-7,780** |
| attn | model (Poisson) | +12 | +23 | +13 | **-35** |
| attn | model (theta=1) | -95 | +35 | +17 | -185 |
| attn | null AUG | -104 | +17 | +47 | -2,884 |
| attn | null near-cognate | -- | -- | -- | -7,780 |

**The "positive on all four datasets" claim holds for mamba4 Poisson ONLY.** attn's Poisson arm is
-35 on A549. D12 makes the all-four claim because it defaults to `--model mamba4`; its docstring
now carries this qualifier. Do not state it model-agnostically.

**null near-cognate exists only for A549** (tryptic): nonspecific digestion of a ~3M-sequence
near-cognate null exceeds MSFragger's ~2e9 peptide cap. Its -7,780 is the largest loss measured
anywhere in this project.

**Do NOT merge with P11.** Across 12 mouse macrophage populations the model's database is above the
GENCODE-only baseline in only 5 of 12 (median -8). Four cell lines here, 12 populations there.

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

**Cost of the A549 removal.** 4 datasets -> **3**, and an **entire ARM disappears**: `null_nc`
(near-cognate null) existed only for A549, because nonspecific digestion of a ~3M-sequence
near-cognate null exceeds MSFragger's ~2e9 peptide cap. The near-cognate null is therefore no
longer represented anywhere in the project's human proteogenomics, and the DRIVER panel now
compares the model against the AUG null only.
