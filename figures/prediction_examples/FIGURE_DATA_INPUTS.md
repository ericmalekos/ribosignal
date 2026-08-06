# FIGURE_DATA_INPUTS: prediction_examples*.png

**What it shows.** Observed vs predicted per-nucleotide Ribo-seq P-site profile for three representative
GOOD predictions -- one uORF, one CDS, one lncRNA ORF -- on the **held-out Hepatocytes** test set (the
tissue is never seen in training, so these are genuine generalization examples, not training fits). Three
rows, each = context view (left) + feature zoom (right).

- **Context views (left panels) + the uORF zoom:** blue fill = observed (real pooled P-sites, as a
  per-transcript fraction); red line = predicted (model softmax profile).
- **Periodicity zooms (right panels of the CDS and lncRNA rows):** one bar per nucleotide of the PREDICTED
  profile, colored by codon frame relative to the ORF start (frame 0 = in-frame, green; 1 orange; 2 purple);
  observed overlaid as a black step line. Frame-0 dominance is the visual signature of 3-nt periodicity.

## Variants on disk

| File | Model | Status |
|---|---|---|
| `prediction_examples.{png,pdf}` | `orf_v2_attn_onehot_holdout_Hepatocytes` | **superseded** -- pre no-Kozak / mm1 / union checkpoint |
| `prediction_examples_mamba4_union.{png,pdf}` | `orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` | current, **primary** shipping model (GPU-only) |
| `prediction_examples_attn_union.{png,pdf}` | `orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` | current, **CPU-compatible** shipping model |
| `prediction_examples_split_mamba4_union.{png,pdf}` | mamba4 union | **split-panel variant**, see below |
| `prediction_examples_split_attn_union.{png,pdf}` | attn union | split-panel variant |

## Two presentations of the same examples

`plot_prediction_examples.py` -- **context + zoom.** Three rows (uORF / CDS / lncRNA), each a
whole-transcript context view (observed fill, predicted line) beside a 60 nt zoom drawing the
PREDICTED profile as frame-coloured bars with observed overlaid as a black step. Compact, and shows
where in the transcript the feature sits.

`plot_prediction_examples_split.py` -- **observed above predicted, both frame-coloured.** Three
columns (uORF / CDS / lncRNA), two rows: observed on top, prediction directly beneath, sharing x
and y. The overlay above is asymmetric -- the observed trace is a single colour, so only the
PREDICTED profile's reading frame is legible and the reader must take the observed periodicity on
trust. Giving each series its own frame-coloured panel makes both independently readable and puts
the frame-0 percentage on each, so the comparison is quantitative:

| model | uORF obs / pred | CDS obs / pred | lncRNA obs / pred |
|---|---|---|---|
| mamba4 union | 88% / 88% | 83% / 86% | 68% / 80% |
| attn union | 90% / 96% | 89% / 86% | 68% / 82% |

Both scripts import selection, gates and zoom-window choice from `plot_prediction_examples`, so for
a given dump the two figures always showcase the same transcripts.

## Rows (re-selected per model with `--auto`)

| Model | uORF | CDS | lncRNA |
|---|---|---|---|
| superseded `attn` | POLR1D `ENST00000399697.7` r=0.993 | RPL32 `ENST00000429711.7` r=0.995 | LINC02693 `ENST00000647872.1` r=0.806 |
| **mamba4 union** | C1orf43 `ENST00000368521.10` r=0.999 | ZDHHC7 `ENST00000564466.5` r=0.993 | C19orf48P `ENST00000641834.3` r=0.728 |
| **attn union** | POLR1D `ENST00000399697.7` r=0.999 | GAS6 `ENST00000881732.1` r=0.992 | C19orf48P `ENST00000641834.3` r=0.838 |

`r` is whole-transcript profile Pearson, except the uORF row, which is scored on the 5'UTR alone (where the
uORF evidence lives; a whole-transcript correlation would be dominated by the far larger CDS). Both current
models independently select the SAME lncRNA exemplar (C19orf48P), a useful consistency check; `attn` fits
it better (0.838 vs 0.728).

For a side-by-side model comparison, pin the same uORF in both with
`--uorf_tx ENST00000399697.7` (POLR1D): it is in the top 3 for both models and has the highest uORF
prominence of the near-ties (40.1%).

## Selection rule, and why Pearson alone is not it

`--auto` re-picks the exemplars for whichever dump is passed, because the best-fitting transcripts are a
property of the model, not of the annotation -- carrying one checkpoint's picks into another model's figure
would show examples chosen for someone else.

**Ranking on Pearson alone fails systematically.** A transcript whose signal is one dominant spike scores
r ~ 1.0 trivially, since observed and predicted agree on the single position carrying all the variance.
Selecting on raw Pearson against the mamba4 dump returned exactly that: a CDS at **r = 1.000** that was a
lone spike over a flat body, and a lncRNA whose spike sat at the edge of its ORF. Both were useless as
illustrations, and both scored *better* than the distributed examples they displaced.

So `autoselect()` gates on what makes a good ILLUSTRATION and uses Pearson only to rank the survivors:

| gate | default | applies to | purpose |
|---|--:|---|---|
| `min_obs` | 200 | all | enough real signal to not be fitting noise |
| `min_nonzero` | 40 | CDS, lncRNA | breadth: the region must be covered, not one position |
| `max_spike` | 0.15 | CDS, lncRNA | no single position may carry more than this share of the region |
| `min_frame0` | 0.50 | CDS, lncRNA | the observed signal must actually be periodic, the thing being shown |
| `min_uorf_share` | 0.15 | uORF | the uORF must be a legible feature of the transcript |

uORFs are exempt from the SPIKE gate -- they are short, and a genuine uORF legitimately reads as a sharp
peak -- but they need `min_uorf_share` instead, because 5'UTR Pearson saturates even harder than the
whole-transcript version. Measured on both current models, ~1,961 uORF candidates yield only 8-13 above
r = 0.99, and the top few differ in the **fourth decimal place**, so the ranking there is noise. Those
near-ties are not interchangeable as illustrations: the uORF's share of transcript signal among them runs
from 46.3% (TBPL1) and 40.1% (POLR1D) down to 5.5% (PTK2) and 0.5% (EIF5A). Before this gate, a 0.0002
difference in r picked PTK2 for the `attn` figure -- a uORF carrying 5.5% of the signal, visually dwarfed
by an unrelated CDS spike -- over C1orf43 at 36.6%. That is why the two models' uORF rows originally looked
so unlike each other.

## Data sources

- **Predicted + observed profiles**: `results/loto/<run>/dropin/pred_profiles.npz` (`pred_flat` = model
  softmax profile, `obs_flat` = real pooled Hepatocytes P-sites, `tx_ids`/`lengths` for the ragged layout).
  Produced by `dump_pred_profiles.py` on the held-out-Hepatocytes LOTO model.
- **CDS bounds** (5'UTR / CDS / 3'UTR for the shading): `data/tx2cds.tsv`.
- **Gene names + biotype** (protein_coding vs lncRNA): `data/tx2biotype.tsv`.
- **ORF windows** (uORF / novel-ORF start-stop, for the shaded region + type): Hepatocytes RiboCode calls
  `expression_context_human/data/ribocode_per_tissue/Hepatocytes/Hepatocytes_collapsed.txt`
  (`ORF_type`, `ORF_tstart`, `ORF_tstop`).

## Regeneration

```
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

# current primary model (writes .png and .pdf side by side)
$PY $NEW/scripts/plot_prediction_examples.py --auto \
  --npz $NEW/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/pred_profiles.npz \
  --title "..." --out $NEW/figures/prediction_examples/prediction_examples_mamba4_union.png

# pin specific transcripts instead of auto-selecting
$PY $NEW/scripts/plot_prediction_examples.py --uorf_tx ENST... --cds_tx ENST... --lnc_tx ENST...
```

Without `--npz` the script reads the superseded dump, and without `--auto` it uses the superseded model's
hardcoded transcript defaults -- pass **both** when regenerating for a current model.

The right-hand periodicity zoom is the frame-aligned 60 nt window inside the ORF/CDS with the most
sustained signal (max count of nonzero observed positions), so it lands on a periodic body region rather
than a single spike (`best_zoom` in the script). The uORF zoom is the full 5'UTR [0, utr5_len).

## Caveats

- These are curated GOOD examples (top of the per-category gated ranking), chosen to illustrate what a
  high-quality prediction looks like -- not a random or representative sample. The distribution of scores
  across all transcripts is in the eval metrics (Task 13 localization, Task 15 profile Pearson).
- Hepatocytes is a deep tissue, so the observed profiles here are clean; on shallower held-out data the
  observed is noisier (see Task 18 depth crossover) even where the prediction is unchanged.
