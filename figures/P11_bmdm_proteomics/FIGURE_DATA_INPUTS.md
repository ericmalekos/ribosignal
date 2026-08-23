# P11 -- mouse macrophage proteogenomic lift, 12 populations: data inputs

**CURRENT: the ATTENTION model at a 30-AA ORF floor** (regenerated 2026-08-15).

| | |
|---|---|
| reports | `proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30/reports/report_<POP>.json` (12) |
| model arm | `model_predicted_attn` |
| ORF floor | **30 aa** -- Rule 5, tryptic whole-cell lysate (trypsin/termini=2, verified from fragger.params) |
| search | 7 arms per population against identical spectra, frozen parameters |

Reproduce the superseded 7-aa mamba panel with
`--source .../pgx_xsubtype/frozen_reports --arm model_predicted`.

## TWO HEADLINE RATIOS FELL WHEN THE FLOOR WAS CORRECTED

|  | 7 aa, mamba (superseded) | **30 aa, attn (current)** |
|---|--:|--:|
| DB size vs null | 156x smaller | **78x smaller** |
| discovery density vs null | 111x higher | **43x higher** |
| model above baseline | 5/12 (median -8) | **6/12 (median +3)** |
| null above baseline | 0/12 (median -1,784) | **0/12 (median -1,130)** |

At 7 aa the null held 235,686 sequences, **~65% of them sub-30-aa ORFs the ORF-call track excludes by
rule**. The correct floor shrinks it to ~84,000 and the ratios fall with it. **The earlier 156x/111x
were inflated by counting sequences that should not have been in the comparison.** The model's own
net benefit moved the other way (+3 vs -8).

## Per-population (attn, 30 aa)

| population | novel (model) | novel (null) | dTotal (model) | dTotal (null) | DB ratio | density ratio |
|---|--:|--:|--:|--:|--:|--:|
| SmallIntestinal | 37 | 53 | +36 | -932 | 80x | 56x |
| Peritoneal | 30 | 38 | +28 | -1,097 | 91x | 72x |
| LargeIntestinal | 19 | 43 | +17 | -1,165 | 92x | 41x |
| BMDM | 34 | 54 | +11 | -1,042 | 81x | 51x |
| Kupffer | 19 | 38 | +9 | -1,329 | 87x | 43x |
| SpleenRecruited | 30 | 49 | +9 | -1,098 | 67x | 41x |
| SpleenResident | 31 | 43 | -3 | -1,092 | 61x | 44x |
| Microglia | 8 | 24 | -6 | -1,151 | 213x | 71x |
| LungRecruited | 20 | 39 | -7 | -1,108 | 66x | 34x |
| LiverRecruited | 24 | 37 | -16 | -1,452 | 66x | 42x |
| LungResident | 27 | 48 | -18 | -1,333 | 76x | 43x |
| RAW264 | 9 | 58 | -25 | -1,214 | 73x | 11x |

## Caveats a caption must carry

- **Cost-neutral, not a gain.** Median +3 of ~55,000 is inside the +-50 frozen-vs-unfrozen
  search-parameter swing measured on BMDM. The uniformly one-sided comparison is model vs NULL
  (0/12), not model vs baseline.
- **Absolute novel counts favour the nulls** and must be shown: they search ~78x more sequences.
- RAW264 is the weakest population (density 11x, -25 total); Microglia the most selective
  (213x smaller DB, only 8 novel peptides). Both retained.
- attn calls ~15-20% fewer ORFs than mamba at matched calibration, yielding a smaller, denser
  database with equal yield -- see results.md 2026-08-15.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/P11_bmdm_proteomics && $PY make_bmdm_proteomics.py
```
