# FIGURE D13 -- the database-design tradeoff

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
