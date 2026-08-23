# Poster figure build — data corrections and requests

Context: I'm building the RiboSignal conference poster (36 × 48 in) from
`poster_package/`. All figures are regenerated at poster type from the
**portable** `*_values.json` / `*.tsv` only — the `make_*.py` generators can't
run off-cluster. It is an **all-`attn`** poster.

Four panels have drafts (`1a` background, `3b` held-out human, `4a` mouse liver,
`5a` proteogenomics). Three are blocked or degraded by the items below.

---

## A. Corrections — things that look wrong in the current package

**A1. `README.md` line 192 (THP-1) carries a mamba4 number inside the attn table.**
The table is headed "**attn, standalone arm**" and gives THP-1 F1 = **0.850**.
`PROJECT_results.md` (GSE208041, `pred_preddepth`) gives:

| model | F1 |
|---|---|
| attn | **0.851** |
| mamba4 | **0.850** |

So the README row appears to have picked up the mamba4 value. Same table's CAR-T
row says 0.904, where results.md gives attn 0.905 / mamba4 0.907 — matches
neither. Hepatocyte (0.909) and iPSC-CM (0.876) both cross-check fine.
Numerically trivial, but a mamba4 number inside a table labelled attn is exactly
the kind of thing this project's conventions exist to prevent. **Please confirm
which is authoritative and correct the source.**

**A2. `E2_class_ceiling` model values silently pool BOTH architectures.**
`E2_values.json` exposes only `model_pred_obsdepth_mean` /
`model_pred_preddepth_mean`. Per `E2/FIGURE_DATA_INPUTS.md` line 11 these are
means over **"all 18 cells per arm (2 models × 9 RNA/Ribo combinations)"**, and
`make_class_ceiling.py` confirms it (`np.mean(mod[c][arm])` over the pooled
list; its own legend says "2 models × 9 RNA/Ribo = 18 per arm").

There is **no `model` field on the file to catch this** — a model-provenance
check passes while the number is mixed. `PROJECT_POSTER_LAYOUT.md` then quotes
"Annotated CDS 0.975 against a ceiling of 0.991–0.995" as though it were a
single-model result. It is a two-architecture mean.

Consequence: I could only plot E2's **ceiling** band (`ribo_vs_ribo`, no model
involved). The strongest version of that panel — *annotated CDS sits at the
ceiling while dORF and iORF fail badly* — is **not currently plottable as attn**.

**A3. Two `FIGURE_DATA_INPUTS.md` files disagree with their own JSON.**
Small, no direction changes, unexplained — but they should agree:

| file | .md says | .json says |
|---|---|---|
| E2, uORF ceiling | 0.693–0.840 | 0.6817–0.8321 |
| E3, novel vs gse243134 recall | 0.581 → 0.709 | 0.5747 → 0.6953 |

**A4. `PROJECT_POSTER_LAYOUT.md` quotes the F2b discovery-density median as 96×.**
`F2b_values.json` and F2b's own `FIGURE_DATA_INPUTS.md` both give **93.6×**
(all 20 records, and the 10 attn records independently). I used 93.6.

---

## B. Data I need, and cannot get from the package

**B1. An attn `profile_exemplars` dump — this is the biggest blocker.**
Panel `3a` is "architecture + inputs + observed vs predicted signal". For the
observed-vs-predicted half there is currently **no source that is both current
and attn**:

- `prediction_examples.tsv` — attn, but the **superseded** pre-no-Kozak/mm1/union
  checkpoint (its own doc says so).
- `profile_exemplars/*.tsv` — current and spike-robust, but **released mamba4**.
- `prediction_examples_split_attn_union.pdf` — current and attn, but a PDF with
  **no underlying data**, so it can only be placed as-is at ~1.37×.

**Request:** run the `profile_exemplars` selector + dump against
`orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`, same schema
as the existing files (`dataset, tx_id, gene, class, zoom_mode, position_nt,
in_called_orf, in_zoom_window, frame_rel_orf, nt, codon, aa, codon_idx_in_orf,
obs_psites, pred_raw, obs_frac, pred_frac`). Human hepatocytes alone is enough.

**B2. Five P-site totals, and two QC rows.**
*(This supersedes the earlier metagene request — `S_metagene/metagene_pooled.tsv`
arrived, and panel `2a` has since been rebuilt as a dataset census: trained on /
held out / dropped. No metagene in it.)*

Panel `2a` plots depth and frame-0 share for every Ribo-seq dataset the project
touches. Here is where each number comes from, and where the holes are:

| rows | depth | frame-0 |
|---|---|---|
| 8 Chothani (incl. Hepatocytes) | `P3_values.json` → `tissues[].psites` | QC |
| 3 mouse livers | `E3_values.json` → `library.*.psites` | E3 |
| T cell, BMDM (NT) | `PROJECT_results.md` regeneration table — **prose** | QC |
| **iPSC-CM** | **nothing** | QC |
| **THP-1 (GSE208041)** | **nothing** | **nothing** |
| **CAR-T (GSE304796)** | **nothing** | **nothing** |
| **brain (dropped)** | **nothing** | `A3_values.json` `period_obs` — different statistic, printed as text |
| GSE39561 (dropped) | zero P-sites, coverage-only pack | `PROJECT_results.md` 40.7% — **prose** |

Four rows draw an em-dash for depth. **THP-1 and CAR-T are the only two datasets
in the whole project with no periodicity number at all**, and Eric asked for both
on the poster specifically, so they are what I most want filled.

**Request 1 — four `SOURCES` entries in `make_riboseq_qc.py`.** This is the cheap
one and it fills every frame-0 hole at once. `SOURCES` is an **explicit list**,
not a glob (I said "glob" in an earlier draft of this document — wrong, and the
distinction matters: nothing gets picked up automatically). It currently holds 8
Chothani + 5 held-out + 2 pgx. Four datasets need adding:

```python
("Brain",   "train",   CANON / "_canonical_Brain/canon_Brain_pre_config.txt"),
("THP-1 (GSE208041)", "heldout", ...),
("CAR-T (GSE304796)", "heldout", ...),
("THP-1 (GSE39561)",  "dropped", ...),
```

**Brain is the one-word fix.** The list is built by a comprehension over a tissue
tuple that already formats the exact path
`results/chothani_regeneration/_canonical_Brain/canon_Brain_pre_config.txt` —
adding `"Brain"` to that tuple is the whole change, and the file landed
2026-08-15.

⚠️ **Two docs disagree about whether brain belongs there,** and you should settle
it rather than take my word: the script's own comment says *"Brain is
intentionally absent (dropped for low periodicity before any pack was built), so
8 is the correct count, not 9"*, while `PROJECT_results.md` (2026-08-15) says the
panel *"becomes complete at 9 tissues on its next run."* The comment predates the
brain fetch. If 9 is now right, the comment needs updating too — it is the kind
of stale guard that will otherwise get brain removed again by whoever reads it.

**GSE39561 should parse fine** even though metaplots selected zero read lengths:
the `ROW` regex reads the commented per-read-length table, which exists
regardless of what was *selected* — that is presumably how the 40.7% pooled
number was derived in the first place. If it does come back with 0 rows, it lands
in `missing` and I will keep using the prose value.

**Request 2 — three P-site totals.** Per-nt over the pack's own universe, the
same quantity as `P3_values.json`'s `psites`: `human_ruizorera`, `gse208041`,
`cart`. `PROJECT_results.md:1119` gives iPSC-CM only as "global-mean depth 743
(~10× shallower than Fibroblast)", too loose to draw.

**Brain depth is the one that needs actual work.** The script's comment says brain
was dropped *"before any pack was built"*, so there is no pack to read a total
from. If a raw P-site count falls out of the canonical alignment cheaply, I will
take it; if it means building a pack for a tissue nobody trains on, say so and I
will leave that cell as an em-dash. **GSE39561 needs nothing here — its zero is
the measurement**, not a gap, and the panel prints it as *zero P-sites*.

Two notes on what I used in the meantime, so you can check them:

- **GSE39561 frame-0 = 40.7%** is the *pooled* frame concentration at annotated
  CDS starts, from "GSE39561: closed as a negative result". I read that as the
  same quantity as QC's `f0_frac`. If it is not, tell me — it is the bar that
  makes the dropped section work, sitting just above the 33% floor where every
  usable library sits at 74.6–87.7%.
- **GSE208041's 91.2%** from the same section is "at its dominant read length",
  not pooled, so I did **not** plot it — it is not comparable to the summed
  values in every other row. That is why THP-1's frame-0 cell is empty.

**B2c. Brain and GSE39561 now have their own section on the poster,** built from
"Brain and GSE39561 fail DIFFERENTLY, and low mapping does not predict which"
(2026-08-15). Please sanity-check the two claims I am printing: metaplots select
periodic read lengths in **5 of 5** brain libraries and **0 of 3** GSE39561, and
both sit at **~30% unique mapping**. Brain is shown as excluded on signal
strength (`period_obs` 0.044) rather than as a broken library, per that section's
correction.

**B2b. `psites_at_cds` in `riboseq_qc_values.json` is misnamed, and it cost me a
figure.** `make_riboseq_qc.py` sets it to the sum of the frame-0/1/2 counts in
RiboCode's `metaplots` table — P-sites in the window around **annotated start
codons**, not across the CDS and not the library. It runs ~0.5% of real depth
(Fibroblast 5.03M there vs 1,069.8M in the pack). I built a first draft of `2a`
plotting it as sequencing depth before catching this.

It is also not usable as a depth *proxy*: its ratio to real pack depth swings
0.33%–0.65% across the eight Chothani tissues, which is enough to **invert the
ES/Fat ordering** (QC says ES 2.10M > Fat 1.86M; real depth says ES 494M < Fat
565M).

**Request:** rename it `psites_at_start_window` (or similar) in the JSON, TSV and
the printed table. No recomputation needed — the number is fine, the name is what
misleads. Worth doing before anyone else reads it as depth.

**B3. An attn-only E2 (`model_vs_ribo_by_class`), if A2 is confirmed.**
Either a per-model breakdown in `E2_values.json`
(`model_pred_preddepth_mean_attn` / `_mamba4`), or the underlying
`model_vs_ribo_by_class.tsv` with its `model` column intact so I can filter.

---

## C. Two things I'd like confirmed rather than assumed

**C1. ORF-length floor mismatch between panels.** The proteogenomics data
(`F2b`) is built on a **7-aa** ORF floor; the ORF-call panels use **30 aa**, and
26.9% of model-Poisson novel peptides rest on ORFs a 30-aa floor would delete.
On the poster these panels sit two bands apart and a reader may take them as the
same ORF set. Is the 30-aa proteomics rebuild (Task #90) likely to land before
the conference? If so I'd rather wait for it than caption around it.

**C2. Annotated-CDS recall for iPSC-CM (0.997) and THP-1 (0.994)** appears only
in the `README.md` prose table — I found no machine-readable attn artifact for
those two arms. Hepatocyte's does cross-check against
`A2_ribocode_dropin/A2_values_attn.json`. Is there a values file I missed?

---

## D. Conventions that make a file usable to me

Not a complaint about existing files, just what saves a round trip:

1. **Put `"model"` on every values file**, even when it looks obvious. A2 was
   the case that made this concrete: an unmarked pooled mean reads as
   single-model.
2. **When a value is an aggregate, say what it aggregates over** in the JSON
   itself, not only in the sibling `.md`.
3. **Keep the `.md` and the `.json` numerically identical** (see A3) — I plot
   the JSON, so a divergent `.md` silently becomes wrong documentation.
4. Portable `*.tsv` beats a PDF every time: anything with a data file, I can
   regenerate at poster type and on the poster's palette. Anything PDF-only has
   to be placed as shipped, at whatever internal type size it was authored with.
