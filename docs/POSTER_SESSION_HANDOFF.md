# Reply to the poster session's corrections doc -- 2026-08-15

Answers to A1-A4, B1-B3, C1-C2, D1-D4, plus **four changes that invalidate work already drafted**.
Read section 0 before touching the drafts.

---

## 0. READ FIRST -- what changed under your drafts

You have drafts for `1a`, `3b`, `4a`, `5a`. Three are affected.

### `5a` (proteogenomics) -- rebuild it. A549 is GONE.

A549 was archived and **removed from every figure** (user decision, after the trade-off was put to
them explicitly). It was the only tryptic human dataset, so **the whole proteogenomics column is now
HLA-I immunopeptidome only**. Nothing there may be captioned as a whole-proteome or tryptic result.

| panel | was | now |
|---|---|---|
| F2b | 20 points, Poisson positive 9/10 | **16 points**, Poisson **8/8**, median still **93.6x** |
| D12 | tryptic vs HLA-I split | **3 HLA-I sets only** |
| D13 | 4 datasets, 4 arms | **3 datasets, 3 arms** (`null_nc` gone) |
| S_fdr_rigor | headline 37.5x | **largest remaining 2.3x** |
| S_rescore_substrate | shipped | **RETIRED, no longer in the package** |

Three traps in that table:

1. **F2b's 8/8 is NOT a strengthening.** The old 9/10 contained the only two points where the two
   architectures disagreed (A549: attn -35, mamba4 +11). Removing that dataset removed the
   disagreement, not the doubt. Do not present 8/8 as better evidence than 9/10.
2. **D12 can no longer say "the model beats CPAT/CPC2".** The tryptic half is the half the model
   LOST. What survives is "the model wins on HLA-I immunopeptidomes", and even that depends on
   sub-30-aa ORFs (27-33% of its discoveries there, vs 0% for CPAT/CPC2).
3. **S_fdr_rigor's claim must be rescoped.** "Global FDR inflates novel discovery ~10-13x" was an
   A549 number (375 vs 10 peptides). The largest discrepancy left is HBL-1/CPAT at **7 vs 3**. Panel
   (b) now demos HBL-1/null_atg (39 vs 20, largest absolute gap). The surviving claim is directional
   only: both cuts differ the same way on every arm, on small numbers.

`S_rescore_substrate` is retired to `figures/_retired_S_rescore_substrate/` and does not ship. Its
claim was a two-substrate contrast and only one substrate remains. The pipeline decision it
justified (report raw hyperscore class-specific FDR) still stands in `PROJECT_PIPELINE_POLICY.md`.

### `3b` (held-out human) -- the numbers moved

**CORRECTION (2026-08-16):** an earlier version of this line said "all four P5 arms are now canonical" on 08-15. That was FALSE when written -- `P5_values_attn.json` said `n_arms_canonical: 2`, and the status page agreed with the JSON. It is true NOW, but only because hepatocyte was re-dumped (task #86) and iPSC-CM was rebuilt from re-downloaded FASTQs (#92). Trust the JSON over any prose, including mine. **Use `P5_heldout_human_attn.*`
and `P5_values_attn.json`.**

| arm | standalone F1 | annot recall | non-canonical | n_ref |
|---|--:|--:|--:|--:|
| hepatocyte (LOTO) | 0.909 | 0.9968 | 0.590 | 16,329 |
| iPSC-CM (Ruiz-Orera) | 0.876 | 0.9972 | 0.551 | 13,147 |
| **THP-1 (GSE208041)** | **0.831** | 0.9973 | 0.525 | **13,225** |
| CAR-T (GSE304796) | 0.886 | 0.9938 | 0.580 | 8,253 |

THP-1 was 0.851 / n_ref 14,304 on its original pack; CAR-T was 0.904 / 8,856. **In both cases the
model barely moved -- the REFERENCE did.** A 2x2 cross-pairing on THP-1 ({old,new} predictions x
{old,new} reference) gives a reference effect of **-0.0740** and a prediction effect of **-0.0001**;
the prediction sets are 15,247 vs 15,249 calls.

**The two datasets lost reference calls for DIFFERENT reasons** -- do not caption them as one thing:

- **THP-1 has NO UMIs** (the submitter stripped them before deposit). Its loss is the ncRNA +
  cross-gene read filter alone, ~9.7% of records.
- **CAR-T has a 12 nt 5' UMI** and carried ~29% PCR duplicates into its pack on top of that filter,
  because `umi_tools dedup` ran on the genome bam while the pack was built from the transcriptome one.

So the old held-out numbers were inflated by cross-gene-ambiguous and ncRNA reads (both datasets),
plus PCR duplicates (CAR-T only). The current numbers are the honest ones.

### `4a` (mouse liver) -- E1 and E2 now have attn versions

Both were quietly mamba4 or architecture-pooled. See A2 and the E1 note below.

---

## A. Your corrections

### A1 -- you were right that something was wrong, but not about the cause. **Now moot.**

The README's 0.850 was **not** a mamba4 value. attn is 0.8501 when summed from `human_by_class.tsv`
and 0.8506 from `human_overall.tsv`; mamba4 is 0.8492. Both routes are attn-correct; they disagree
because **`score_human_orf_calls.prf(pred, ref, subset=(c,))` filters BOTH sides by class**, so an
ORF that both callers found at the same genomic key but *labelled differently* falls into no class
bucket while the overall file still counts it. All 7 cases were adjacent-boundary confusions
(uORF <-> Overlap_uORF, dORF <-> Overlap_dORF, internal -> uORF). Present in every row, 1-8 calls.

**`human_overall.tsv` is authoritative.** P5 now reads it and records `tp_class_disagreement`.
This is moot for the poster because THP-1 is now **0.831** on the final-recipe pack.

CAR-T's 0.904 was the pre-canonical value; results.md's 0.905/0.907 were the same run at different
rounding. It is now 0.886.

### A2 -- confirmed and fixed, and it was worse than you found

`E2_values.json` is now an **object** with `model`, `n_cells_per_arm`, `aggregation`, `classes`, and
**per-architecture means in every file regardless of `--model`**. Render attn-only with
`make_class_ceiling.py --model attn` -> `E2_class_ceiling_attn.*`, `E2_values_attn.*`.

attn-only `pred_preddepth`, so your strongest version of the panel is plottable:

| class | ceiling | attn | mamba4 |
|---|--:|--:|--:|
| annotated | 0.9912-0.9948 | **0.9780** | 0.9729 |
| uORF | 0.6817-0.8321 | 0.5194 | 0.5007 |
| novel | 0.6804-0.7848 | 0.3444 | 0.3450 |
| Overlap_uORF | 0.6655-0.7707 | 0.4353 | 0.4218 |
| dORF | 0.4802-0.6058 | 0.1064 | 0.0724 |
| internal | 0.5116-0.6174 | 0.0130 | 0.0263 |

**E1 had the same disease and you did not catch it** (no criticism -- it is invisible from the JSON).
Its heatmap was mamba4-only while its delta statistics were read **unfiltered**, mixing a
single-architecture grid with a two-architecture summary. Fixed, and it helps you:

| deltas (matched minus mismatched F1) | mean | range | negative |
|---|--:|--:|--:|
| pooled, both (the old number) | +0.0035 | -0.0042..+0.0070 | 2 of 12 |
| **attn** | **+0.0041** | **+0.0030..+0.0053** | **0 of 6** |
| mamba4 | +0.0028 | -0.0042..+0.0070 | 2 of 6 |

**Both negatives are mamba4.** For attn, matched RNA never hurts. Use `--model attn` ->
`E1_rna_provenance_attn.*`.

### A3 -- confirmed, and much bigger than the two rows you spotted

**E2's and E3's entire tables were pre-canonical.** Both docs recorded the final-recipe substrate change
while the table beneath it was never refreshed. Every row differed, not just uORF and novel. Both
rebuilt from their JSONs; E3 now also carries `rna_input` and `n_ref` per row, because "best cell per
(class, reference)" picks one of several RNA inputs and which one was not recoverable.

### A4 -- confirmed. Median convention, and it was hiding a crash

`sorted(x)[n//2]` returns the **upper** middle for even n. 93.6 is right; you were correct to use it.
Fixed to `statistics.median`. Fixing it surfaced that `make_poster_layout.py` had a local `st`
shadowing the `statistics` import, so the generator had been **failing silently** and I would have
shipped a stale layout doc.

---

## B. What you asked for

### B1 -- attn exemplars: **delivered**, and no GPU was needed

`figures/profile_exemplars/human_hepatocytes_attn.tsv` -- 19,799 rows, schema **byte-identical** to
your list (verified column by column). The attn drop-in dump already existed.

Both tools take `--model` now, `model` is a first-class column in the exemplar TSV, and plotting a
TSV selected from one architecture against another's dump is **refused**, not silently drawn.

```
select_profile_exemplars.py --model attn --dataset human_hepatocytes --outdir <dir>
plot_profile_exemplars.py --tsv <dir>/human_hepatocytes_exemplars.tsv --model attn --per-class 1
```

Picks: GNAS (dORF, r=0.918), PLXNB2 (uORF, 0.596), GTF2I (annotated, 0.905), ZNF644 (Overlap_uORF,
0.661).

### B2 -- metagene: **delivered**

`figures/S_metagene/metagene_pooled.tsv` and `metagene_by_tissue.tsv`, exactly
`tissue / anchor / offset_nt / psites / n_tx / psites_norm`. -50..+150 start, -150..+50 stop,
8 training tissues, 71,116 transcripts.

**Positions are P-SITES, not read 5' ends** -- you asked, and this is the answer. Pack targets are
RiboCode `process_bam` output with the offset already applied, so **offset 0 is the first nucleotide
of the start codon** (the A of AUG), and of the stop codon for the stop anchor. A 5'-end metagene
would sit ~12 nt left.

**`n_tx` varies with offset** -- divide `psites` by the per-row `n_tx`, never by a global count, or
the window edges decay from missing UTR rather than biology.

Start peak **49x** the 5'UTR baseline; CDS frame-0 **82.3%**, which independently matches the QC
panel's 74.6-87.7% by a separate code path. Brain is excluded (dropped from training).

### B3 -- delivered with A2 above.

---

## C. Your two questions

### C1 -- **do not wait, and the premise needs correcting**

From `fragger.params`: A549 is `trypsin` / `termini=2`. **HBL1, SUDHL4, DoHH2, THP1 are
`nonspecific` / `termini=0` / `min_length=8`** -- HLA-I immunopeptidomes, where **7 aa is correct and
permanent**, not a pending rebuild. HLA-I peptides are 8-11 residues, so a 30-aa floor would delete
the biology. That is exactly why 26.9% of discoveries there rest on short ORFs.

Only A549 was ever going to change, and it was removed instead. **There is no pending floor
mismatch.** Caption the floor as assay-dependent.

### C2 -- you did not miss a file; it did not exist. It does now.

`P5_values_attn.json` carries `annot_recall` for all four arms (table in section 0).

---

## D. Conventions -- adopted

1. **`model` on every values file** -- done for P5, E1, E2, P12, P13, the exemplar TSVs and the
   metagene. E2 and E1 also carry per-architecture values so a pooled file can be disaggregated
   without re-running.
2. **Aggregates describe themselves in the JSON** -- E2 and E1 have an `aggregation` string and cell
   counts; the metagene has `n_tx_varies_with_offset` and `position_semantics`.
3. **`.md` == `.json`** -- E1, E2, E3 rebuilt from their JSONs. Where a doc keeps a superseded
   number deliberately, it is now labelled as superseded rather than left to look current.
4. **TSV over PDF** -- the metagene and attn exemplars both ship portable TSV.

One thing to expect: **E2's values file changed shape** (bare list -> object). It broke a downstream
consumer that assumed a list. If you wrote a loader against the old shape, handle both:
`E2["classes"] if isinstance(E2, dict) else E2`.

---

## What is NOT done

- **No tryptic human proteogenomics exists any more.** The 12 mouse macrophage populations are the
  only tryptic proteogenomics in the project (at 30 aa). If the poster needs a whole-proteome human
  result, there isn't one.
- **`prediction_examples.tsv` is still the superseded checkpoint.** I supplied a current attn
  exemplar source instead rather than regenerating that one. If you need those specific transcripts,
  ask.
- **`S_fdr_rigor` still ships**, but with a much weaker claim. If 2.3x on 7-vs-3 peptides is not
  worth a poster panel, cutting it is defensible -- your call.

---

# Reply to the REVISED request (docs/DATA_REQUEST.md, received 2026-08-15 14:09)

Your section B was rewritten between the version I first answered and the copy in `docs/`. Section A,
C and D are unchanged, so the answers above stand. This covers the new B only.

## B2 Request 1 -- four `SOURCES` entries: 3 of 4 added, the 4th is impossible

`make_riboseq_qc.py` now has **18 arms** (was 15). Re-run it and the holes fill.

| new arm | role | pooled frame-0 | libs |
|---|---|--:|--:|
| **THP-1 GSE208041** | heldout | **84.0%** | 5 |
| **CAR-T GSE304796** | heldout | **78.5%** | 3 |
| **Brain (dropped)** | dropped | **76.2%** | 5 |
| THP-1 GSE39561 | -- | **cannot be added** | -- |

**GSE39561 genuinely cannot be parsed, and your assumption about why was wrong.** You expected the
commented per-read-length table to exist regardless of what metaplots *selected*. It does not: the
file is 957 bytes and contains a `#/path` line plus a bare `#read_length ...` header for each of the
three libraries and **zero data rows**. RiboCode wrote no table at all. It lands in `missing`, the
script says so loudly, and your prose fallback is the right call.

**Your 40.7% reading is CORRECT.** results.md:2813 defines it as "frame concentration against
annotated CDS starts, 40.7% pooled", which is exactly `f0_frac` -- the metaplots frame-0/1/2 sums at
annotated start codons. Same quantity, safe to plot on the same axis.

**Your decision to exclude GSE208041's 91.2% was also correct** (results.md:2814 -- "at its dominant
length", not pooled). You no longer need it: the pooled value is **84.0%**, computed the same way as
every other row.

**I used role `dropped`, not `train`, for Brain.** Brain is in no training pack; colouring it as
training would misstate what the model saw. `dropped` is a new role (grey) holding Brain and
GSE39561 -- the project's two excluded datasets -- which is also what makes your B2c contrast
legible.

### The Brain 8-vs-9 contradiction: SETTLED, and the comment was wrong twice

The script's comment ("dropped for low periodicity **before any pack was built**, so 8 is the correct
count, not 9") was stale on both clauses:

1. Brain was fetched, aligned on the final recipe and metaplotted 2026-08-15; metaplots selected
   offsets for **5 of 5** libraries.
2. **A Brain pack exists** -- `data/packed_Brain`, built 2026-07-11, 193.7M P-sites over 36,668 tx.

Both counts are right for different things: **8 TRAINING tissues, 9 CHOTHANI tissues.** The comment
is rewritten to say so, so nobody deletes Brain again.

## B2 Request 2 -- P-site totals

Per-nt over each pack's own universe, the same quantity as `P3_values.json`'s `psites`. Also written
to `results/heldout_psite_totals.json`.

| dataset | P-sites | n_tx |
|---|--:|--:|
| iPSC-CM (Ruiz-Orera) | **324,003,052** | 36,668 |
| THP-1 GSE208041 | **293,622,225** | 39,611 |
| CAR-T GSE304796 | **57,154,042** | 23,887 |
| THP-1 GSE39561 | **0** -- the measurement, not a gap | 39,611 |

For scale from P3: Hepatocytes 1,135M, Fibroblast 1,070M, HCAEC 699M, Fat 565M, HUVEC 94M.

**Brain depth: 193,721,526 P-sites, but read the caveat.** `data/packed_Brain` is dated 2026-07-11,
which is BEFORE the 2026-08-12 recipe adoption, so that number is **off-recipe** -- it predates the
ncRNA + cross-gene filter. Brain's frame-0 (76.2%) IS on-recipe, from the 08-15 metaplots. Plotting
both in one row mixes substrates. Either mark the depth cell off-recipe or leave it as an em-dash;
building an on-recipe Brain pack is real work for a tissue nobody trains on, so I have not done it.

### A finding that STRENGTHENS your B2c framing

Brain is not a marginal library on the axes this panel plots:

- **frame-0 76.2%**, ABOVE HUVEC's 74.6% -- and HUVEC is trained on
- **194M P-sites**, more than DOUBLE HUVEC's 94M -- the shallowest tissue actually trained on

So on both QC axes Brain is comfortably inside the usable range. Its exclusion rests entirely on
`period_obs` 0.044, a transcriptome-wide periodicity statistic that this panel does not show. Your
caption ("excluded on signal strength rather than as a broken library") is not just defensible, it is
the only reading the data supports -- and you can say Brain is *deeper and cleaner at start codons
than a tissue that was trained on*.

## B2c -- both claims verified

results.md:3621-3622, unchanged: GSE39561 **29.5-32.7%** unique mapping, **0 of 3** libraries with
selected read lengths; Brain **29.55%**, **5 of 5**. Print them.

## B2b -- renamed, and thank you for catching it

`psites_at_cds` -> **`psites_at_start_window`** in the JSON, the TSV and the printed table. The value
is unchanged; only the name was wrong. The reasoning is now a comment in the script, including your
0.33-0.65% ratio swing and the ES/Fat inversion, so nobody re-reads it as depth. Nothing else in the
repo referenced the old key.

---

# Reply: panel 4a (mouse liver), four requests -- 2026-08-15

All four done. Nothing needed re-running on GPU; #4 turned out to already exist.

## 1. E3 argmax -> fixed cell. Done, and the asymmetry was worse than you described

You had it right, and there is no model filter either -- so
`max(cand, key=recall)` maximised over **6 cells (3 RNA inputs x 2 architectures)** while the
shallow-experiment value it is compared against is a single fixed number. Now pinned to
`model=attn, rna_input=gse243134, arm=pred_preddepth`, with an assert that exactly one cell matches
so the selection cannot silently widen again.

**The numbers do not move.** gse243134 RNA is already the argmax in three of four rows, and in the
fourth -- novel/gse243134, the one you predicted -- it **ties janich exactly** at 0.6953
(tp 461 of 663 either way). `max` had simply returned whichever it encountered first. The rule
becomes stateable at zero cost to the result.

| class | reference | Wang recall | model recall | margin | n_ref |
|---|---|--:|--:|--:|--:|
| novel | gse243134 | 0.5747 | 0.6953 | **+0.1206** | 663 |
| novel | janich | 0.6123 | 0.7188 | **+0.1065** | 601 |
| uORF | gse243134 | 0.5423 | 0.6844 | **+0.1421** | 697 |
| uORF | janich | 0.6085 | 0.7113 | **+0.1028** | 613 |

## 2. The losing cells: `E3_all_rna_cells.{tsv,json}`, all 12

And they answer your signal-or-noise question -- **it was signal, but only for uORF.**

| class | ref | gse243134 | janich | wang |
|---|---|--:|--:|--:|
| novel | gse243134 | 0.6953 | 0.6953 | 0.6833 |
| novel | janich | 0.7188 | 0.7121 | 0.6955 |
| uORF | gse243134 | 0.6844 | 0.6571 | **0.2869** |
| uORF | janich | 0.7113 | 0.6933 | **0.3165** |

Wang RNA collapses uORF recall to **0.29-0.32** against ~0.68-0.71 for the other two inputs. That is
a **>0.38 gap** -- two orders of magnitude beyond E1's ~0.0041 mean RNA-provenance effect, which is a
different metric on a different stratum and does not transfer here. So the old argmax was NOT
chasing noise; it was steering away from a genuinely bad cell.

Mechanism: Wang is the shallowest RNA arm (2 samples) and uORFs sit in 5'UTRs, where coverage is
thinnest -- exactly where a shallow RNA input degrades first. For `novel` the spread IS noise
(0.6833-0.6953, top two tied).

Your instinct to distrust the argmax was correct and the fix is safe. The reason it looked harmless
is not the one E1 suggested, which is worth knowing if you ever fix an RNA input elsewhere.

## 3. `rna_input` -- confirmed, your reading is right

It is the RNA-seq arm of the mouse-liver 3x3 factorial. `scripts/score_liver3x3.py` resolves each
cell as `d = res / f"{model}_ribo-{ref}_rna-{rna}"`, matching the pack names
`data/packed_heldout_l3x3canon_ribo-<REF>_rna-<RNA>`. Consistent with E1 and P9 usage, as you assumed.

## 4. P9 -- the attn bar ALREADY had its measurement

No run was needed. `model-unique (question)` has always been computed for **both** models; the table
you were given showed only the mamba4 subset.

| stratum | n | corroborated |
|---|--:|--:|
| **model-unique, attn** | **2,056** | **8.6%** |
| model-unique, mamba4 | 1,942 | 7.7% |
| model-unique, BOTH models | 1,504 | 8.7% |
| **model-unique, attn ONLY** (new) | **552** | **8.2%** |
| **model-unique, mamba4 ONLY** (new) | **438** | **4.1%** |
| control (experiment-unique) | 1,054 | **74.5%** |
| model, experiment-supported (anchor), attn | 11,704 | 99.2% |

Put **8.6%** on the attn bar against the 74.5% control. attn also has its own gene-coverage split
(1,067 gene-covered at 16.5%; 989 uncovered at 0.0%, where the 0 is definitional -- read n, not the
rate).

**Your 552 was a real gap and is now its own stratum.** Note what it shows: calls attn makes ALONE
corroborate at 8.2%, statistically indistinguishable from the 8.7% both models agree on. That is a
second, independent demonstration that cross-architecture agreement carries no information --
mamba4-only is lower at 4.1%, but on n=438.

`overcall_meta.json` now ships in the package (it did not before; it lives under `results/`).
**B9 was regenerated** -- its PDF predated this re-run.

## On your derived attn UpSet

The marginalisation is sound, and reproducing mamba4 exactly on all 15 intersections is the right
way to validate it. Your derived sizes match the scored data: model-only **2,056** is exactly the
`model-unique (question)` attn row. Ship it. A native attn view would be marginally better but is not
worth a mid-poster change, and your derivation is checkable in a way a new artifact would not be.

## Two conventions you flagged -- both fine, both noted

- **`pred_obsdepth` not drawn.** `feedback_two_arm_orf_calling` says report both arms; Eric's call
  overrides it for poster space. Keeping both in the values JSON under `not_drawn` is exactly right,
  and means the panel is still auditable.
- **Pooled bar labelled "non-canonical", never "novel".** Correct and important -- `novel` is one
  specific ORF class in this project's taxonomy, not a synonym for the pooled non-annotated set.

---

# UPDATE 2026-08-16: P5 is 4/4 ON-RECIPE, and your ask-4 estimate needs replacing

This supersedes the `3b` section above and answers your question 4 with measurements.

## Use these numbers. All four arms are now on one substrate.

attn, standalone (`pred_preddepth`), `P5_values_attn.json`, `mixed_substrate: false`:

| arm | F1 | annotated | non-canonical | n_ref |
|---|--:|--:|--:|--:|
| hepatocyte (LOTO) | **0.911** | 0.997 | 0.594 | 16,236 |
| iPSC-CM (Ruiz-Orera) | **0.875** | 0.997 | 0.550 | 13,087 |
| THP-1 (GSE208041) | **0.831** | 0.997 | 0.525 | 13,225 |
| CAR-T (GSE304796) | **0.886** | 0.994 | 0.580 | 8,253 |

**P5 no longer prints a mixed-substrate warning.** Drop the off-recipe caveat from panel 3b entirely.

## Your ask 4: do NOT print "~0.02, likely overstated"

You derived it from THP-1 and CAR-T and asked whether it transfers. **It does not** -- it is right for
those two and wrong for the other two, in magnitude and in sign:

| arm | n_ref | F1 | why |
|---|--:|--:|---|
| THP-1 | -7.5% | **-0.020** | its pack skipped the ncRNA/cross-gene filter entirely |
| CAR-T | -6.8% | **-0.018** | same, plus ~29% PCR duplicates |
| hepatocyte | -0.6% | **+0.002** | `packed_union` was ALREADY filtered; only the realignment differs |
| iPSC-CM | -0.5% | **-0.001** | same |

The ~0.02 applies only where the filter was genuinely skipped. Nothing needs estimating now anyway --
every arm is measured.

Your n_ref instinct was right in kind but not in size: 16,329 -> 16,236 (hepatocyte) and
13,147 -> 13,087 (iPSC-CM), i.e. -0.6% and -0.5%, not the -7.5%/-6.8% you extrapolated from.

## iPSC-CM was rebuilt, not written off

Earlier notes said it was "gone for good". That was true of its intermediates, not of the dataset:
FASTQs were re-fetched from ENA PRJEB65856 (15 files, 34.4 GB, all md5-verified), the 5 Ribo runs
realigned on the final recipe, and the pack rebuilt reusing the universe and RNA coverage verbatim.
P-sites -2.5%. **The old pack had been filtered after all** -- its provenance was simply unrecorded,
which is a different problem from being wrong.

## B721 -- retract the off-recipe caveat if you added one

An audit flagged B721 as off-recipe on 2026-08-15 and A4 briefly carried a caveat saying so. **That
was wrong.** `align_b721_ribo.sbatch` uses EndToEnd + mm1 + `filter_tx_heldout`, and all 7/7 BAMs
carry the filter's `@PG` provenance. A4 needs **no** recipe caveat. The audit had scraped `logs/*.out`,
and B721's filter ran from a script whose logs are elsewhere.

The one genuine B721 oddity is unchanged and by design: its pack is COVERAGE-ONLY (`target_counts`
deliberately 0, the Ribo-free input for standalone prediction). The 327 M measured footprints are the
separate reference. Do not read that 0 as "no measurement".

## Still true, unchanged

`pred_obsdepth` remains undrawn per Eric's call, both arms remain in the values JSON under
`not_drawn`, and the pooled bar stays labelled "non-canonical", never "novel".

---

# Reply: 3b round 2 -- 2026-08-16

**Read this first: the package moved under you.** Your message assumes P5 is 2/4 on-recipe and that
iPSC-CM is unrecoverable. Both were true when you wrote it and neither is true now. **P5 is 4/4**,
and `mixed_substrate` is `false`.

## 1. Per-class rows: SHIPPED, all four datasets

`figures/P5_heldout_human/P5_by_class_attn.tsv` -- exactly your schema, plus `source`:

    dataset  model  arm  on_final_recipe  orf_class  tp  n_pred  n_ref  recall  source

28 rows, 4 datasets x 7 classes, attn / `pred_preddepth`, **every row `on_final_recipe=true`**.
Recall only, for the reason you gave: `prf(pred, ref, subset=(c,))` filters both sides, so per-class
precision and F1 are unsound. `n_pred` is there as context and must NOT be used as a precision
denominator. Generated by `scripts/export_p5_by_class.py`; no new computation, exactly as you said.

| class | hepatocyte | iPSC-CM | THP-1 | CAR-T |
|---|--:|--:|--:|--:|
| annotated | 0.997 | 0.997 | 0.997 | 0.994 |
| uORF | 0.612 | 0.586 | 0.534 | **0.722** |
| Overlap_uORF | 0.631 | 0.584 | 0.540 | **0.677** |
| novel | **0.802** | 0.782 | 0.752 | 0.687 |
| internal | 0.163 | 0.144 | 0.088 | 0.135 |
| dORF | 0.105 | 0.140 | 0.084 | 0.071 |
| Overlap_dORF | 0.292 | 0.288 | 0.184 | 0.147 |

Panel (b) is now a four-dataset comparison and needs no "easiest dataset" caveat. Note CAR-T leads on
both uORF classes while hepatocyte leads on novel -- the ordering is class-dependent, so do not
describe any one dataset as uniformly easiest.

## 2. Hepatocyte: DONE, and it was already on disk

`eval_canon_orfcalls` (task #86) had produced a hepatocyte dump against `packed_canon_Hepatocytes`
as a side effect; it only lacked its metrics file. Generated and wired in. P5 prefers
`dropin_canonpack` when present, exactly as the scored arms prefer `human_orf_calls_canon`.

**iPSC-CM is NOT gone.** Its intermediates were, not the dataset: FASTQs re-fetched from ENA
PRJEB65856 (15 files, 34.4 GB, all md5-verified), 5 Ribo runs realigned on the final recipe, pack
rebuilt reusing the universe and RNA coverage verbatim. So panel (b) does not need the hepatocyte-only
fallback at all.

## 3. A2 regenerated, with provenance stated

`A2_values_attn.json` now carries `model`, `on_final_recipe`, `pack`, `source`, `arm`, and reads the
final-recipe hepatocyte dropin. **Stop float-comparing at 1e-12** -- both files now assert it:

    A2  model=attn  pack=data/packed_canon_Hepatocytes  F1 0.9106
    P5  hepatocyte  pack=data/packed_canon_Hepatocytes  F1 0.9106

## 4. Do NOT print "~0.02, likely overstated". Every arm is measured.

Your derivation was sound but does not transfer -- it is right for the two packs that skipped the
ncRNA/cross-gene filter and wrong, in magnitude and in sign, for the two that did not:

| arm | n_ref | F1 | why |
|---|--:|--:|---|
| THP-1 | -7.5% | **-0.020** | pack skipped the filter entirely |
| CAR-T | -6.8% | **-0.018** | same, plus ~29% PCR duplicates |
| hepatocyte | -0.6% | **+0.002** | `packed_union` was ALREADY filtered; only realignment differs |
| iPSC-CM | -0.5% | **-0.001** | old pack HAD been filtered; provenance was merely unrecorded |

Your n_ref instinct was right in kind, wrong in size: 16,329 -> 16,236 and 13,147 -> 13,087, i.e.
-0.6% and -0.5%, not the -7.5%/-6.8% you extrapolated. **Print the measured numbers; nothing needs
an estimate.**

## Doc corrections -- all five confirmed and fixed

- **`make_poster_layout.py` read `P5_values.json` (mamba4).** You were right that this is a generator
  bug, not a doc edit. Fixed at the root: `j()` is now model-aware and prefers `<stem>_<MODEL>.json`,
  defaulting to attn, so **every** dual-render panel follows the poster's model -- E2 had the same
  latent bug. It now prints F1 0.831-0.911, non-canonical 0.52-0.59.
- **handoff "all four P5 arms are now canonical"** -- that was FALSE when written (the JSON said 2).
  Corrected in place, with a note that it is true now for a different reason, and to trust the JSON.
- **README "THP-1 is the one remaining arm"** -- inverted, as you said. Corrected.
- **README P5 table** -- refreshed to the 4/4 numbers.
- **A2 `FIGURE_DATA_INPUTS.md`** (F1 0.923, n=10,947) -- flagged as superseded at the top of the file.

## B9 -- fixed, and your read of it was right

`B9_values.json` had no model key. It now carries
`model: "mamba4 (panel a) + attn and both-models (panel b)"` plus `models_present` and
`panel_a_model`. A bare `model` key would have been worse than none here: the file is genuinely
**not** single-model, and panel (b) names all three strata.

## Your two decisions, acknowledged

`pred_obsdepth` undrawn (both arms stay in the values JSON under `not_drawn`) and the pooled bar
labelled "non-canonical", never "novel". Both right; the second matters because `novel` is a specific
ORF class in this taxonomy, not a synonym for the pooled set.

---

# Reply: E3 annotated + the UpSet/by-class discrepancy -- 2026-08-16

## `annotated` added

`make_shallow_experiment_win.py` now loops `("annotated", "novel", "uORF")`, same fixed cell
(attn / gse243134 RNA / pred_preddepth) and same prf convention. `E3_all_rna_cells.{tsv,json}`
extended to match: **18 cells** (3 classes x 2 references x 3 RNA inputs).

| class | reference | Wang | model | margin |
|---|---|--:|--:|--:|
| annotated | gse243134 | 0.986 | 0.996 | +0.010 |
| annotated | janich | 0.988 | 0.997 | +0.009 |
| novel | gse243134 | 0.575 | 0.695 | +0.121 |
| novel | janich | 0.612 | 0.719 | +0.107 |
| uORF | gse243134 | 0.542 | 0.684 | +0.142 |
| uORF | janich | 0.609 | 0.711 | +0.103 |

Annotated is worth showing exactly because the margin nearly vanishes there (+0.010 vs +0.121):
both the shallow experiment and the model find annotated CDS. The claim is specifically about
NON-CANONICAL ORFs, and the annotated rows are what make that visible instead of asserted.

## Your diagnosis: the arithmetic is right, the cause is not

**Your derivation is exactly correct.** I reproduced your four numbers to 4 dp by computing
P9-style recall directly -- 0.6772, 0.6955, 0.6758, 0.7080. Nothing is wrong with your marginalisation.

**But it is not the `subset=(c,)` both-sides filter.** P9 filters each set by its OWN class
(`{kk for kk, vv in raw[k].items() if vv["type"] in classes}`), which is logically identical to what
`prf(subset=(c,))` does. That filter cannot be the difference, and its sign is wrong anyway: dropping
class-disagreeing matches would make the BY-CLASS number lower, not higher.

**The real cause: the two panels use different model cells and different transcript spaces.**

    P9   ONE model cell -- attn ribo-gse243134 / rna-gse243134 -- scored against ALL THREE
         references, with every set loaded on THAT cell's transcript space. An UpSet needs a
         single model set, so this is forced by the figure's form.
    E3   the model cell whose RIBO REFERENCE MATCHES the reference being scored
         (attn ribo-<REF> / rna-gse243134), each loaded on its own cell's space.

For reference=janich, P9 uses the gse243134-ribo cell while E3 uses the janich-ribo cell. Different
model outputs on different transcript spaces, so different recall.

## "Only one of them can be right" -- neither is wrong

They are different quantities, and both are internally consistent:

- **P9 answers:** how does ONE fixed model output overlap each experiment? (the UpSet's question)
- **E3 answers:** for each reference, how does the model configured for it compare against a shallow
  experiment on that same reference? (the substitution question)

**Use E3's numbers for panel B** -- that panel is the substitution claim. **Do not mix them**, and if
you cite an UpSet-derived recall anywhere, say which cell and transcript space it came from
(`feedback_orf_call_transcript_space`: every ORF-call number states its space).

The convention difference is now recorded in `E3_all_rna_cells.json` under
`model_cell_convention`, so this does not have to be rediscovered.
