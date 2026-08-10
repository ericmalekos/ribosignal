# MAIN FIGURES 1 and 2 -- composition and provenance

**Generator:** `make_main_figures.py` (cas12a env; needs `pypdf`, installed 2026-08-09).
**Outputs:** `Figure1_validity.pdf`, `Figure2_utility.pdf`, plus a `*_manifest.json` per figure
recording every panel's source path, its mtime, and where it was placed.

Panels are placed as their ORIGINAL PDF pages, scaled and translated -- **not re-rastered**. Text stays
selectable and lines stay vector. Panel letters are drawn into a matplotlib overlay of the same page
size and merged on top.

## Why this is a script and not a hand-assembled file

Every panel regenerates from result JSONs whenever the model changes, and in August 2026 five of the
six Fig 1 panels were found to have been built from a superseded checkpoint AND from the wrong model.
A figure assembled by hand goes stale exactly the same way and cannot be re-derived. This can be re-run
after any retrain, and the manifest makes a stale panel visible as a stale mtime.

The layout is a plain data structure (`FIGURES` in the generator): rows of `(panel_path,
relative_width)`. Re-specifying a figure means editing that list and nothing else.

## Figure 1 -- validity (7.20 x 10.19 in, 5 panels)

| letter | panel | source | model |
|---|---|---|---|
| a | observed vs predicted per-nt profiles (uORF / CDS / lncRNA exemplars) | `prediction_examples/` | mamba4 |
| b | drop-in ORF calling, F1 0.913 | `A2_ribocode_dropin/` | mamba4 |
| c | localization AUROC vs the observed Ribo-seq ceiling | `A1_localization_ceiling/` | mamba4 |
| d | depth crossover: predicting beats measuring below ~29M P-sites | `depth_crossover/` | (deployed) |
| e | 9-fold LOTO spread tracks held-out target quality, not a generalization gradient | `A3_loto_spread/` | **pre-union recipe** |

## Figure 2 -- utility (7.20 x 8.67 in, 3 panels)

| letter | panel | source | datasets |
|---|---|---|---|
| a | model/null discovery-density ratio, all points above 1.0 | `F2b_discovery_forest/` | 5 |
| b | model vs CPAT / CPC2, split by assay | `D12_cpat_cpc2/` | 4 |
| c | DB-design tradeoff: gain, cost, driver | `D13_discovery_errorbars/` | 3 |

## Deviations from FIGURES_PLAN, and why

The plan's Fig 1 list was (a) LOTO spread, (b) exemplars, (c) depth crossover, (d) one-hot ~= token-emb.

- **Panel (d) of the plan, "one-hot ~= token-emb", is NOT included.** No standalone figure was ever
  built for it, and its numbers in the plan (onehot 0.923 / rinalmo 0.915 / orthrus 0.916) are
  pre-union. `A2_ribocode_dropin` is used instead: it is the strongest validity statement available
  and is on the released model. If the encoding panel is wanted in the main figure, it needs building
  from the Task 41 union FM sweep first.
- **Panel order differs**: exemplars lead, LOTO spread comes last. The exemplars are the visual
  centerpiece and read better at the top; A3 is the one panel not on the released recipe.
- **Fig 2 has no panel (c) "highlighted novel PSMs" (D14)**, and no cross-species MS arm. Task 62 was
  superseded 2026-08-09 at 5 of 6 datasets -- see below.

## Known limitations of this assembly (state them, do not hide them)

1. **Panel 1e (A3) is on the pre-union recipe.** There is no union 9-fold LOTO and building one is nine
   retrains. Its CLAIM -- that the fold spread tracks held-out target periodicity at r=0.81 rather than
   being a generalization gradient -- is a property of the design and survives. Its ABSOLUTE Pearson
   values are checkpoint-dependent and must not be quoted beside released-model numbers. Label it.
2. **Panels 1d and 1e carry their own internal sub-panel labels**, which collide visually with the
   outer a-e letters. Fixing that means editing those two generators' internal labelling.
3. **1e is small** (3.43 x 1.44 in) and its axis text is at the edge of legibility at print size.
   Equal widths in that row were chosen so panel (d), which is nearly square, stayed readable; the
   two cannot both be comfortable in one row at 7.2 in.
4. **Fig 2's three panels cover 5, 4 and 3 datasets respectively.** That is a property of the
   underlying analyses, not the assembly, but a reader will notice the counts differ.
5. **Figure 1 is 10.2 in tall.** That exceeds a single printed page at full width and will need either
   a taller format or a panel moved to supplemental.

## Task 62 supersession (recorded here because it changes Fig 2)

Closed 2026-08-09 by user instruction at **5 of 6** immunopeptidome datasets: A549 (tryptic) plus
HBL-1, SU-DHL-4, DoHH2 and THP-1 (HLA-I). Five is enough for panel 2a to show a reproducible pattern
(20 points, all above 1.0, median 93.6x) rather than a one-off. The two not delivered are blocked on
data availability rather than effort:

- **B721.221 MS** -- MassIVE FTP port 21 is filtered from prism. B721 still contributes the drop-in
  ground-truth check (F1 0.666 = 75% of the 0.889 split-half ceiling).
- **A mouse immunopeptidome** -- PXD008733 has 39 raw files but no matched RNA-seq; Rospo CT26 MMRd has
  ENA RNA-seq but no raw MS accession. Matched RNA-seq is non-negotiable because it IS the model input.

**Consequence for wording:** Fig 2 has no cross-species MS arm. The cross-species claim rests on
Ribo-seq (mouse Wang / Janich / GSE243134 drop-in), not on mass spectrometry, and must say so.

## Regenerate

```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/main
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
$PY make_main_figures.py                       # both figures at 7.2 in wide
$PY make_main_figures.py --only Figure2_utility --width 6.9
```

Rebuild the underlying panels first if the model changed; each panel folder has its own
`FIGURE_DATA_INPUTS.md` with its regeneration command. `figures/README.md` "Checkpoint provenance"
lists which panels are on which model and which take `FIG_MODEL`.
