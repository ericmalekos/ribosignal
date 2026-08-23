# Mouse-liver 3x3 factorial: RNA universe x Ribo-seq dataset x model

**Question.** Within one fixed transcript universe, how do the three mouse-liver Ribo-seq experiments
agree with each other, and how close do the two released models get -- and does the answer change when
the universe itself is redefined by a different RNA set?

This replaces the earlier three-way comparison, which scored each model arm inside a *different*
universe and against a *different* reference, so RNA input, transcript space and reference all moved
at once. Nothing could be attributed.

## Design

Three RNA sets X in {Wang, Janich, GSE243134}. For each:

1. **Universe U_X** = X's salmon TPM >= 1 transcript list (13,041 / 14,887 / 17,225 tx).
2. **RNA coverage** inside U_X always comes from X. This is the axis under test.
3. **Three Ribo packs** on U_X, differing ONLY in `target_counts`: Wang, Janich, GSE243134 P-sites.
4. **RiboCode on each pack's real counts** -> 3 observed call sets inside U_X.
5. **Two model dumps** (attn, mamba4) on U_X with X's coverage. The three Ribo packs share tx_order,
   lengths, coverage and ORF track, so the prediction is identical across them -- 2 dumps per
   universe, not 6.
6. **Score inside U_X**: 3 observed-vs-observed pairs + 2 models x 3 references = 12 rows.

36 rows total, each with precision / recall / F1 for overall, annotated CDS, and non-canonical.

Totals: 9 packs, 3 ORF tracks, 6 dumps, 15 RiboCode runs.

## Uniform processing: everything lands on the UNION, everything slices from it

Mixing "slice Wang's pack" with "re-pool the others from BAMs" would put two code paths behind one
table. Instead every dataset is first brought onto one common reference space and every universe is a
slice of it.

**UNION = U_Wang + U_Janich + U_GSE243134 = 22,974 tx** (`data/heldout_refs/mouse_liver_union3_tx.txt`).
Verified: 100.00% covered by Wang's full-annotation pack, 0 chrM transcripts
(feedback_riboseq_exclude_mito_genes).

| dataset | how it reaches UNION | why |
|---|---|---|
| Wang | slice its existing full-annotation pack (221,835 tx) | no BAMs and no P-site hd5 survive; the pack is per-nt so the slice is exact |
| Janich | re-pool 7 archived BAMs | its pack covers only 63-78% of the other universes |
| GSE243134 | re-pool 40 archived BAMs | same |

RNA coverage reaches UNION the same way, then slices to U_X.

After that, **all nine Ribo packs and all three coverage tracks are produced by one operation: slice
UNION -> U_X.** No dataset gets a privileged code path.

## Validation gate before any of it is trusted

Re-pool Janich onto its OWN universe from BAMs and require the resulting `target_counts.npy` to match
its existing pack. If the BAM path and the historical path disagree, the factorial is measuring
processing differences rather than biology, and everything downstream is void. This runs first and
blocks the rest.

## Stated limitations

1. **Depth is not controlled.** The three Ribo sets stay at natural depth -- 60.6M / 94.5M / 133.4M
   P-sites over scored transcripts. Observed-vs-observed differences therefore confound lab and
   protocol with depth, and the depth ordering already tracks the ceiling ordering. Deferred by
   decision to task 75 (matched-depth arm); revisit if the depth ordering reappears here.
2. **Six of the nine cells are cross-study pairings** -- RNA from one study, Ribo from another, so
   different animals. Only the diagonal has matched biology. Accepted deliberately as "closer to real
   application", where a model is handed whatever RNA a user has. Every table marks diagonal vs
   off-diagonal; they are never pooled.
3. **Universe size is not held constant** across the three panels (13,041 / 14,887 / 17,225). Counts
   are comparable within a panel, not across panels. F1 across panels is comparable only in the sense
   that each is measured against its own within-universe reference.

## Execution order

```
0  validate: Janich BAMs -> own universe == existing pack        (BLOCKING)
1  UNION sources: Wang slice; Janich + GSE243134 from BAMs        3 jobs
2  ORF track on UNION (--kozak none), sliced per universe         1 job
3  9 packs = slice(UNION source, U_X) x RNA coverage X            1 array
4  15 RiboCode runs (9 real + 6 model)                            2 arrays
5  score -> results/mouse_liver_3x3/ + one TSV per panel
```

---

## 2026-08-10: step-0 gate FAILED, and the design changes because of it

`liver3x3_00_validate` (job 36624625) re-pooled the five Janich Ribo runs from re-aligned BAMs onto
the Janich universe and compared `target_counts.npy` against the surviving pack.

```
gold      94,483,232 P-sites over 37,658,319 nt
rebuilt   93,602,292 P-sites over 37,658,319 nt
DIFFER    442,926 nt differ, total abs diff 1,595,364 (1.6885% of signal)
```

Bar was <0.1%. This is 17x that, so: FAIL.

**Not a data problem.** The run accessions match the gold pack exactly (SRR1930188, 189, 193, 196,
197; `n_ribo_samples: 5` in both). The re-fetched FASTQs are md5-verified against ENA.

**Cause: RiboCode `metaplots` re-derives which read lengths carry 3-nt periodicity, per sample, from
the data.** From the rebuilt config:

```
SRR1930188 ... P-siteReadLength 29,30,31   P-siteLocations 12,12,12
SRR1930189 ... P-siteReadLength 29,30      P-siteLocations 12,12
```

Read length 31 is 13.19% of SRR1930188. A borderline length flipping in or out of the selected set
moves ~1% of total P-sites across five samples, which is the size of the observed divergence. The
gold pack's `metaplots` config lived in the deleted `expression_context_human` directory, so it
cannot be replayed. **Byte-continuity with the historical Janich pack is unrecoverable.**

### What changes

Chasing the old pack is the wrong target. The factorial does not need continuity with a historical
artifact; it needs all three Ribo-seq datasets on ONE pipeline, so that differences in the table are
differences between experiments rather than between code paths. The original design already conceded
an asymmetry here (Janich and GSE243134 re-pooled from BAMs, Wang sliced from a pre-existing
full-annotation pack) and the gate is precisely what surfaced how large that asymmetry is.

Wang turns out to be fixable: `data/external/wang2021_mouse/` still has both RNA FASTQs
(md5-verified). Both Wang RIBO runs had been deleted -- the same loss pattern as Janich, RNA kept and
Ribo gone. They are being re-fetched (`data/external/wang2021_ribo_refetch/`, ~4.7 GB).

Revised step 0/1: align Wang Ribo from FASTQ with the same recipe as Janich and GSE243134, then pool
ALL THREE from BAMs through `prepare_from_bams.py`. No dataset is sliced from a legacy pack.

Per-dataset `metaplots` offsets are kept (NOT pinned to a shared read-length set). That is the
standard RiboCode workflow and the correct choice: the three datasets have genuinely different
footprint-length distributions, and forcing one config across them would impose a different bias than
the one it removed.

### Stated limitation, carried into any writeup

Because `metaplots` selection is data-dependent, the rebuilt Janich pack is ~1.7% different from the
historical one. Factorial numbers are internally consistent across the nine cells but are NOT
directly comparable to previously reported Janich numbers computed on the old pack (e.g.
`results/liver_released/`). Those earlier results remain valid on their own pack; they are a
different measurement, not a contradicted one.

---

## RESOLVED 2026-08-13: the final-recipe rebuild is COMPLETE

> The OFF-RECIPE banner immediately below is retained as history. It is no longer the current state.
> Appended, not rewritten -- the original plan text is untouched, and the pre-append copy is at
> `mouse_liver_3x3_plan.md.bak.2026-08-14`.

The rebuild it describes as "in progress" finished: 40 Ribo alignments -> 3 pools -> 9 packs ->
18 dumps -> rescore, 0 failures. Current tables are `results/mouse_liver_3x3_canon/scored/`; the
off-recipe originals are archived at `results/_archive_offrecipe_2026_08_13/mouse_liver_3x3/`.

What changed, and what did not:

| | off-recipe | canonical |
|---|--:|--:|
| observed calls (janich / gse243134 / wang) | 12,804 / 13,145 / 12,104 | 12,289 / 12,584 / 11,646 |
| matched-minus-mismatched RNA, mean | +0.0030 | **+0.0035** |
| `pred_obsdepth` F1 range | 0.895-0.926 | 0.897-0.925 |
| `pred_preddepth` F1 range | 0.864-0.880 | 0.855-0.871 |

The headline conclusion is unchanged: matching the RNA input to the Ribo reference is worth ~0.003 F1
and goes negative in the same 2 of 12 cells.

**The one substantive finding the rebuild produced** is that final-recipe alignment strips 22.3% of
NOVEL and 12.5% of uORF reference calls while touching annotated CDS by only 1.4%. The model's own
standalone predictions are byte-identical across the switch, so the earlier non-canonical F1 was
INFLATED by the model being credited for matching alignment artifacts. See results.md, "The
off-recipe non-canonical F1 was INFLATED by alignment artifacts". Quote canonical non-canonical
numbers only, and never against previously reported ones.

Figures E1/E2/E3, `figures/S_riboseq_qc/`, and the tutorial's mouse tables all now read the canonical
scored tables.

## OFF-RECIPE, flagged 2026-08-13

> **OFF-RECIPE (flagged 2026-08-13).** The mouse-liver 3x3 numbers below were computed on packs built
> without `--alignEndsType EndToEnd` and without the ncRNA + cross-gene filter, both mandatory under
> `docs/PIPELINE_POLICY.md`. Re-pooling the same reads canonically removes 3.5-4.9% of P-sites
> (janich 110.0M -> 106.1M, wang 78.3M -> 74.5M, gse243134 150.5M -> 143.2M). Comparisons AMONG the
> nine cells stay valid (they share the deviation); absolute F1 is provisional until rescoring.

Canonical rebuild in progress: `data/liver3x3/pool_*_canon` and `pack_canon_*`, built by
`scripts/riboseq_align.sbatch` (EndToEnd + mm1 + filter). The originals are retained so the
effect on reported F1 can be measured rather than assumed.
