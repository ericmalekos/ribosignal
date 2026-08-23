# PROFILE EXEMPLARS -- observed vs predicted per-nt profiles, spike-robust selection

**What these are:** replacements for `figures/prediction_examples/`, chosen so the observed/predicted
agreement is NOT carried by a single enormous peak, and extended to mouse.

## The problem with picking exemplars on r alone

A per-transcript Pearson r can be produced almost entirely by one nucleotide. Measured on the two
exemplars the old panel used:

| old exemplar | r_full | r after dropping top 3 | **top1_frac** |
|---|--:|--:|--:|
| RPL32 CDS `ENST00000429711.7` | 0.984 | 0.821 | **0.538** |
| POLR1D uORF `ENST00000399697.7` | 0.945 | 0.792 | **0.356** |
| LINC02693 lncRNA `ENST00000647872.1` | 0.694 | 0.637 | 0.052 |

**A single nucleotide carries 54% of RPL32's observed P-sites.** Two of the three old exemplars were
showcasing one spike; the lncRNA was fine.

## Selection

`scripts/select_profile_exemplars.py` scores every transcript with >= 200 P-sites and a called ORF
>= 90 nt, in five dumps (human Hepatocytes + Ruiz-Orera, mouse Wang / Janich-decon / GSE243134), all
on the released mamba4 model. Per transcript it computes:

| column | meaning | want |
|---|---|---|
| `r_full` | Pearson over the whole transcript | high |
| `r_drop1`, `r_drop3` | Pearson with the top 1 / top 3 highest-OBSERVED positions deleted | high |
| `r_orf`, `r_drop_orf3` | the same, restricted to the called ORF (the region a figure zooms into) | high |
| `top1_frac` | `obs.max() / obs.sum()` -- spike dominance | **low** |
| `spread90` | fraction of nonzero positions needed to reach 90% of the signal | high |
| `f0_obs_in_orf` | observed in-frame fraction inside the ORF -- real periodicity | high |
| `score` | `r_drop_orf3 x (1 - top1_frac/0.35) x (f0/0.5)`, all terms clipped to [0,1] | ranking only |

The composite is deliberately crude and **every term is in the TSV**, so it can be ignored and the
table re-ranked on any single column. ORF class and transcript biotype come from each dump's own
`real_collapsed.txt`, so the exemplars are the same calls every other analysis uses, and a lncRNA
exemplar is distinguishable from a novel ORF on a coding transcript.

Ranked output: `results/profile_exemplars/{<dataset>,all}_exemplars.tsv`.

## Plotting

`scripts/plot_profile_exemplars.py` takes rows from those TSVs -- so the figure is re-specifiable
without editing code:

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model

# best exemplar per ORF class for one dataset
$PY scripts/plot_profile_exemplars.py \
    --tsv results/profile_exemplars/mouse_gse243134_liver_exemplars.tsv --per-class 1 \
    --out figures/profile_exemplars/mouse_gse243134_liver

# hand-picked transcripts, mixed species, from the combined table
$PY scripts/plot_profile_exemplars.py --tsv results/profile_exemplars/all_exemplars.tsv \
    --tx ENST00000908311.1 ENSMUST00000022416.15 --out figures/profile_exemplars/mixed
```

### Which window the ORF zoom shows (`--zoom-mode`)

The right-hand panel shows `--zoom-width` nt (default 60 = 20 codons) of the called ORF. If the ORF is
shorter than that, all of it is shown. Otherwise `--zoom-mode` decides which window, and **the choice
measurably changes how good the model looks**. Mean Pearson inside the displayed window across the 12
current exemplars, against 0.574 for the whole ORF:

| mode | mean r in window | vs whole ORF | |
|---|--:|--:|---|
| `start` **(default)** | 0.491 | **-0.082** | first 60 nt from the initiator |
| `median` | 0.511 | -0.062 | window of median observed density |
| `center` | 0.583 | +0.009 | middle of the ORF |
| `dense` | 0.636 | **+0.063** | highest observed signal -- the original behaviour |

`dense` was the only behaviour until 2026-08-09 and it is a **best-case view**: it carries 2.4-3.1x the
signal of an average window on long ORFs and beat the whole-ORF correlation in 9 of 12 exemplars, by up
to +0.26. It selects on observed counts only, never on agreement, but agreement tracks depth, so the
effect is the same.

`start` is now the default because it involves **no selection on the data at all**: same window
definition for every transcript, always shows the initiator codon, nothing to cherry-pick. It is not
merely neutral -- it is conservative, running 0.082 BELOW the whole-ORF correlation, because the 5' end
of an ORF has lower and less regular coverage. `center` is the closest to representative (+0.009) if a
like-for-like window is wanted instead.

**The mode is written into the panel title and into every row of the TSV (`zoom_mode`)**, so a
best-case window can never appear unlabelled.

### Sequence track

`--seq-track` puts the sequence under the ORF zoom: **`codon` (nt triplets, the default)**, `aa`
(single-letter amino acids), `both`, or `none`. Text sits at the CENTRE of its codon so each entry
lines up with the three bars it summarises. The initiator codon is green and bold; stop codons are
red and bold. Thinning is width-aware -- triplets are three characters against a 3-nt slot and so
crowd about three times sooner than single letters -- so codons thin past 22 in view and amino acids
past 26 (every 2nd, then every 5th).

**Numeric position ticks are removed from the ORF zoom** when a sequence track is drawn: the codons
occupy that space, and the panel title already states the window (`ORF zoom (1329-1389)`). Pass
`--zoom-xticks` to keep them. The whole-transcript panel on the left always keeps its position axis --
without it there is no way to see where in the transcript the ORF sits.

Sequences come from **the exact universe FASTA each dump was produced with**, not a generic
species-level file:

| dataset | FASTA |
|---|---|
| human Hepatocytes, human Ruiz-Orera | `data/union_universe.fa` |
| mouse Wang | `data/heldout_refs/mouse_wang_expressed_universe.fa` |
| mouse Janich | `data/heldout_refs/mouse_janich_decon_universe.fa` |
| mouse GSE243134 | `data/heldout_refs/mouse_gse243134_universe.fa` |

Using a generic FASTA would silently mis-pair transcripts whose universe differs between arms, and the
resulting sequence would look perfectly plausible.

**Automatic correctness check.** Every exemplar's called ORF is translated and checked for INTERNAL
stop codons. If the ORF coordinates, the reading frame, or the transcript-to-FASTA pairing were wrong
for a transcript, its translation would be full of stops; a clean `M...*` with none in between
confirms all three at once. The script warns loudly per transcript. All current exemplars across all
five datasets pass with zero warnings -- e.g. the Gm8883 lncRNA ORF translates to 127 codons,
`MNFQWESQRAVRANRRNDRV...PRLSHGC*`, no internal stops. The check runs on the translation regardless of
which `--seq-track` mode is displayed, so switching to codons does not switch it off.

`--style mirror` (default) draws observed upward and predicted downward. In `overlay` the two agree so
closely on a good exemplar that the predicted line completely hides the observed fill and the panel
reads as if only one track were plotted; `overlay` is kept for when that is the point.

Each run writes `<out>.tsv` with **one row per plotted nucleotide** (`obs_psites`, `pred_raw`,
`obs_frac`, `pred_frac`, `in_called_orf`, `in_zoom_window`, `frame_rel_orf`, `nt`, `codon`, `aa`,
`codon_idx_in_orf`), so the panel can be rebuilt in any tool without re-running the model.

`nt` is populated for every position, but `codon` / `aa` / `codon_idx_in_orf` ONLY inside the called
ORF. Extending the ORF's reading frame across the whole transcript would emit codons at negative
indices that translate to plausible-looking residues while meaning nothing, and anyone filtering on
`aa != ""` would silently pick up 5'UTR "protein".

## ATTN files (2026-08-16) -- what an all-attn build should read

| File | Model | Contents |
|---|---|---|
| `human_hepatocytes_attn_allclasses.tsv` | attn | per-nt, **7 ORF classes x top 3**, `zoom_mode=start` |
| `human_hepatocytes_attn_allclasses_center.tsv` | attn | the same 21 transcripts at `zoom_mode=center` |
| `human_hepatocytes_attn_ranked.tsv` | attn | the **ranked selector output**, all 13,715 scoreable transcripts |
| `human_ruizorera_attn_ranked.tsv` | attn | ditto, 10,895 transcripts |
| `human_hepatocytes_attn.tsv` | attn | the earlier 4-class file; a strict SUBSET of `_allclasses` |

`human_hepatocytes_attn_allclasses.tsv` is a **strict superset** of `human_hepatocytes_attn.tsv`:
all 4 original transcripts are re-picked as their class winners, and their 19,799 shared per-nt
rows are value-identical (verified 2026-08-16 across `obs_psites`, `pred_raw`, `obs_frac`,
`pred_frac`, `in_called_orf`, `frame_rel_orf`, `nt`, `codon`, `aa`). Swapping one for the other
changes no existing number.

### The missing ncORF row was a `--top` truncation, not a data gap

The shipped exemplar TSVs were written with the selector's default `--top 200`, which keeps the
200 highest-scoring rows **per dataset across all classes**. Annotated CDS dominate that ranking,
so every scarce class was cut. That is why `human_hepatocytes_attn.tsv` had no `novel` / lncRNA row.

Re-running with `--top 0` shows the candidates were always there. Human Hepatocytes, attn, all
transcripts clearing the gates (>= 200 P-sites, ORF >= 90 nt):

| class | n clearing the gates |
|---|--:|
| annotated | 12,172 |
| uORF | 737 |
| **lncRNA** | **413** |
| Overlap_uORF | 296 |
| internal | 41 |
| dORF | 32 |
| Overlap_dORF | 24 |

**Always pass `--top 0` when the output feeds a per-class figure.** The default silently
restricts a class-resolved table to whatever survives a global ranking.

### The lncRNA class winner is not the one with the best `r_full`

The selector's `score` is `r_drop_orf3 * (1 - min(top1_frac/0.35, 1)) * min(f0/0.5, 1)` -- it
rewards distributed agreement and punishes spike dominance, deliberately ignoring `r_full`. On the
lncRNA class the two rankings disagree sharply:

| tx_id | gene | r_full | r_drop_orf3 | top1_frac | P-sites | score | rank by |
|---|---|--:|--:|--:|--:|--:|---|
| `ENST00000687951.3` | TRAF3IP2-AS1 | 0.175 | 0.688 | 0.054 | 202 | 0.581 | **score #1** |
| `ENST00000660020.2` | SNHG29 | 0.851 | 0.859 | -- | 5,020 | -- | best `r_full` of the strong set |
| `ENST00000641834.3` | C19orf48P | 0.838 | 0.639 | -- | 1,030 | -- | the `--auto` pick in `prediction_examples_attn_union` |

TRAF3IP2-AS1 wins on score with `r_full = 0.175`: agreement is good inside the ORF and poor across
the rest of the transcript, and the score only looks inside. Put beside an annotated exemplar at
`r_full = 0.905` it will read as a failure. **The class winner is a defensible choice but a bad
poster row.** Pick from `human_hepatocytes_attn_ranked.tsv` on a stated column and say which
column -- do not present the score winner as "the best" without qualification.

## Current picks (best per class, released mamba4)

| dataset | class | transcript | gene | r_full | r_drop_orf3 | top1_frac |
|---|---|---|---|--:|--:|--:|
| human Hepatocytes | annotated | ENST00000908311.1 | MORF4L2 | 0.936 | 0.930 | 0.030 |
| human Hepatocytes | dORF | ENST00000453292.7 | GNAS | 0.932 | 0.915 | 0.030 |
| human Hepatocytes | uORF | ENST00000555462.5 | PPP4R3A | 0.721 | 0.908 | 0.029 |
| human Ruiz-Orera | lncRNA | ENST00000778459.1 | SNHG29 | 0.674 | 0.762 | 0.095 |
| mouse GSE243134 | annotated | ENSMUST00000022416.15 | Anxa11 | 0.636 | 0.643 | 0.025 |
| mouse GSE243134 | lncRNA | ENSMUST00000188676.3 | Gm8883 | 0.495 | 0.545 | 0.047 |
| mouse Janich | annotated | ENSMUST00000040828.7 | H2-Ab1 | 0.588 | 0.661 | 0.047 |
| mouse Wang | annotated | ENSMUST00000082437.11 | Selenof | 0.596 | 0.611 | 0.078 |

Note the several rows where `r_drop_orf3` EXCEEDS `r_full` (PPP4R3A 0.721 -> 0.908). Whole-transcript
r is dragged down by regions outside the called ORF; agreement inside the ORF, with its top peaks
removed, is much stronger. Those are the most honest exemplars in the table.

`top1_frac` for every pick is 0.02-0.10, against 0.54 for the RPL32 exemplar being replaced.

## Regenerate

```bash
$PY scripts/select_profile_exemplars.py            # rescore all five dumps
$PY scripts/plot_profile_exemplars.py --tsv <...>  # replot
```
Re-run the selector after any retrain: it reads the dumps, so its rankings follow the model.
