# Reply to the poster session: well 3a (architecture + observed vs predicted)

All five items are answered. Four shipped new portable files; one is a "no, and here is why" that
turned out to be a "yes, and here is the file".

---

## 0. Confirmed. `prediction_examples.tsv` is mamba4. Doc fixed.

**You were right, and the mechanism is worse than a typo.** `plot_prediction_examples.py` sets:

```python
MODEL = os.environ.get("FIG_MODEL", "mamba4")
SUF   = "" if MODEL == "mamba4" else f"_{MODEL}"
```

So **the unsuffixed filename IS the mamba4 output, by construction.** The `_attn` suffix is what
marks attn. There is no scenario in which `prediction_examples.tsv` is attn.

Your Pearson re-derivation is reproduced exactly. I get 0.9987 / 0.9839 / 0.6941, and agree the
lncRNA row is decisive.

**One extra thing you could not have seen from the outside.** That file is a *hybrid*, which is why
its numbers match neither doc row. The Aug 9 run passed neither `--npz` nor `--auto`, so it took
the **current mamba4 union dump** with the **superseded checkpoint's hardcoded transcript picks**
(POLR1D / RPL32 / LINC02693, the argparse defaults). The doc's "superseded attn" row quotes
0.993 / 0.995 / 0.806, which is a third set of numbers again. Both rows were wrong.

`FIGURE_DATA_INPUTS.md` now has a corrected variant table with a `Transcript picks` column,
`prediction_examples_attn.{png,pdf}` added, the derivation table above, and a marker that the
unsuffixed files are a hybrid and should not be used.

**On the AGGF1 / POLR1D discrepancy: I cannot re-derive it, and this is worth knowing.** The
script's `CALLS` path,
`expression_context_human/data/ribocode_per_tissue/Hepatocytes/Hepatocytes_collapsed.txt`, has been
deleted. `--auto` raises `FileNotFoundError`, so `prediction_examples*` **cannot be regenerated at
all today**. I have marked what the PNG draws as authoritative over that table row. (18 scripts and
docs reference that deleted annotation tree. It is a bigger open problem than we had logged.)

Agreed on not regenerating. Fixing the doc is the right fix, and it is done.

---

## 1. ⭐ The ncORF exemplar exists. There are 413 of them.

**The gap was a `--top` truncation, not a data gap, and not a gate.**

The shipped exemplar TSVs were written with the selector's default `--top 200`, which keeps the 200
highest-scoring rows **per dataset across all classes**. Annotated CDS dominate that ranking, so
every scarce class got cut. Re-running with `--top 0` on the attn dumps, human Hepatocytes:

| class | n clearing the gates (>= 200 P-sites, ORF >= 90 nt) |
|---|--:|
| annotated | 12,172 |
| uORF | 737 |
| **lncRNA** | **413** |
| Overlap_uORF | 296 |
| internal | 41 |
| dORF | 32 |
| Overlap_dORF | 24 |

So you get all four extra classes you asked about, including `internal` and `Overlap_dORF`.

### Shipped

| File | Contents |
|---|---|
| `figures/profile_exemplars/human_hepatocytes_attn_allclasses.tsv` | per-nt, **7 classes x top 3 = 21 transcripts**, `zoom_mode=start`, your exact schema |
| `figures/profile_exemplars/human_hepatocytes_attn_allclasses_center.tsv` | same 21 at `zoom_mode=center` (this is also §5) |

**It is a strict superset of `human_hepatocytes_attn.tsv`.** All 4 of your current transcripts are
re-picked as their class winners, and their **19,799 shared per-nt rows are value-identical**
(checked across `obs_psites`, `pred_raw`, `obs_frac`, `pred_frac`, `in_called_orf`,
`frame_rel_orf`, `nt`, `codon`, `aa`). Swapping files changes no number you have already drawn. I
left the original in place rather than overwriting it, so nothing breaks under you.

### Read this before you pick the lncRNA row

**The class winner has `r_full = 0.175`.** The selector's `score` is
`r_drop_orf3 * (1 - min(top1_frac/0.35, 1)) * min(f0/0.5, 1)` -- it deliberately ignores `r_full`
and only looks inside the ORF. On lncRNAs the two rankings come apart badly:

| tx_id | gene | r_full | r_drop_orf3 | P-sites | why it appears |
|---|---|--:|--:|--:|---|
| `ENST00000687951.3` | TRAF3IP2-AS1 | **0.175** | 0.688 | 202 | score #1, so `--per-class` picks it |
| `ENST00000660020.2` | SNHG29 | **0.851** | 0.859 | 5,020 | strong on both, 25x the depth |
| `ENST00000641834.3` | C19orf48P | 0.838 | 0.639 | 1,030 | the `--auto` pick in `prediction_examples_attn_union` |

TRAF3IP2-AS1 agrees well inside the ORF and poorly across the rest of the transcript. Beside an
annotated exemplar at `r_full = 0.905` it will read as a failure and undercut exactly the point the
panel is making. **SNHG29 `ENST00000660020.2` is the row I would put on the poster** -- it is
strong on both metrics and is not a depth outlier. All three are in the shipped per-nt file, so it
costs you nothing to choose. State which column you ranked on in the caption.

Ruiz-Orera attn is also shipped and does have `novel`, `internal` and `Overlap_dORF`, if you want a
second dataset. It needs labelling as a different dataset, as you said.

---

## 2. Ranked selector output -- shipped, and bigger than you asked for

| File | Rows |
|---|--:|
| `figures/profile_exemplars/human_hepatocytes_attn_ranked.tsv` | 13,715 |
| `figures/profile_exemplars/human_ruizorera_attn_ranked.tsv` | 10,895 |

Every column you listed, plus `species`, `orf_tstart`, `orf_tstop`, `tx_len`, `psites_in_orf`,
`r_drop1`, `r_drop3`, `r_orf`. These are the **full scoreable set**, not a top-N, so "the exemplars
are the best in their class" is now auditable against everything they beat.

Note the class column is `orf_type`, and `tx_biotype` is separate: RiboCode labels an ORF on a
lncRNA as `novel`, so **lncRNA is `orf_type == "novel"` AND `tx_biotype == "lncRNA"`**. That is how
the plotting script derives its `lncRNA` class, and you will need the same rule to join.

Your reconciliation of the convention is correct: **r** = `r_full` (whole transcript), **r(top-3
removed)** = `r_drop_orf3` (inside the called ORF, top 3 highest-**observed** deleted).

---

## 3. `arch_spec_attn.json` -- shipped, and your derivation is right

At `figures/arch_attn/arch_spec_attn.json`. Everything is **read from artifacts, not transcribed**:
parameter count summed from `best.pt`'s state_dict, channel and kernel shapes from the tensor
shapes, hyperparameters from the run's `args.json`. So it drifts with the code instead of silently
going stale.

**Confirmations:**

- **`n_params = 5,071,106`** exactly. Breakdown now in the JSON: body blocks 3,947,520 / attention
  1,054,208 / count head 66,305 / in_proj 2,816 / profile head 257.
- **Receptive field 4,093 nt, derivation correct.** `1 + 2*(k-1)*sum(dilations)` = `1 + 4*1023`.
  Two k=3 convs per block each widening by `2d`; `in_proj` is k=1 and adds nothing.
- **`d_ff = 512` correct**, but it is `channels * 2`, not an independent hyperparameter.
- `n_heads` 8, `n_attn_layers` 2, `norm_first` true, GroupNorm(8), dropout 0.1, `count_weight` 0.1
  all confirmed.

### Two corrections that change what the panel should say

**1. "Receptive field ~4 kb" understates the attn model.** 4,093 nt is the **convolutional**
receptive field. The mixer is full self-attention masked for padding alone
(`self.attn(xt, src_key_padding_mask=~mask)`) -- no causal mask, no local window. **The attn
model's effective context is the entire transcript.** Suggested label: *"convolutional receptive
field 4,093 nt; attention is global over the transcript."*

**2. The coverage channel is not raw depth.** It is
`log1p(coverage / global_mean_coverage)` (`cov_norm=global_mean`) -- per-nt RNA depth divided by
the pack's own global mean *before* `log1p`. Channel order is one-hot(4) | orf_track(5) |
coverage(1), coverage always last.

Also in the JSON: the count head takes the masked mean of body features **concatenated with
`log1p(total RNA coverage)`**, which is not obvious from a block diagram, and
`peakiness_weight = 0.0` confirming the anti-smoothing term is not in this checkpoint.

On `arch_attn_values.json`'s `ENST00000958170.1` at Pearson 0.911: agreed, no portable source.
Leave it unquoted. I have not dropped it yet in case it is load-bearing somewhere.

---

## 4. All 10 start weights -- shipped. Your reimplementation is exactly right.

At `figures/arch_attn/orf_track_spec.json`.

```
ATG 1.0   CTG 0.5   GTG 0.35   ACG 0.35   TTG 0.3
ATA 0.25  ATT 0.2   ATC 0.2    AAG 0.15   AGG 0.15
```

Your guess at the missing four was correct on identity: UUG 0.3, AUA 0.25, AUU 0.2, AUC 0.2.

**Your 24-nt matrix reproduces exactly.** I ran the real `build_orf_track.orf_track()` on
`ATGCATGCTGAAGCCTAACTAGCC` and compared it against what you printed. All five channels match,
including `f2` being empty and the graded starts at 0 / 4 / 7 / 10. The JSON carries the computed
arrays so your assert can check against numbers instead of against a calibrated PNG.

**Your reimplementation is correct rather than approximate, and here is why**: the deployed model
is built with `--kozak none`, so the start channel is the raw codon weight with **no context
multiplier**. Omitting Kozak is the right behaviour, not a simplification.

**On the 24x10 matrix: only 24x9 is derivable, and I did not fake the tenth.** Channels 1-4 are
one-hot and 5-9 are the ORF track, both sequence-derived and both shipped. The 10th is RNA
coverage, which is *data* -- it has no defined value for a synthetic display window. Draw it as a
placeholder or omit it. The JSON says this explicitly under `rna_coverage_note` so nobody later
fills it with plausible numbers.

---

## 5. Center-zoom variant -- shipped

`figures/profile_exemplars/human_hepatocytes_attn_allclasses_center.tsv`, same 21 transcripts.

Verified it differs from the `start` file **only** in `in_zoom_window` (2,118 rows move); every
value column is identical. So the two are a clean like-for-like pair and the choice is purely
which window to frame, with no other variable moving.

Your reasoning for keeping `start` as the default is sound and I would keep it as the primary.

---

## Not done

- Nothing was regenerated in `prediction_examples/`, per your request.
- No new model run, retraining or eval.
- `make_arch_panel.py` left non-portable, as you asked.
