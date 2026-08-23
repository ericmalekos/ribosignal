# A4 -- B721.221: predicted vs MEASURED ORF calls, against the split-half ceiling

**What it shows:** the project's only NON-NULL benchmark, which had no figure until 2026-08-14
(results.md Task 66, landed 2026-08-07). Every other comparison scores the model against a null,
another selection rule, or a single-dataset reference. Here both arms are in one cell line.

## Panel -> data

| element | source |
|---|---|
| model arms (4) | `proteogenomics/data/B721_pilot/dropin_compare/b721_dropin_metrics.json` |
| split-half ceiling | `.../dropin_compare/ceiling_AvsB/b721_dropin_metrics.json` |
| prediction input | Sarkizova RNA-seq (the model never sees Ribo-seq) |
| measured reference | Ouspenskaia Ribo-seq, **327 M unique footprints** |

Genomic keying, restricted to the model's 11,527-gene space: 21,974 measured calls genome-wide ->
**16,331 in scope**. `min_len` 90 nt, pval 0.05.

## The ceiling, and why panel A alone would mislead

RiboCode was run independently on two **disjoint halves of the same Ribo-seq**, split by HLA allele
so each half holds complete libraries. Those halves reproduce each other at:

| | F1 | recall annotated | recall non-canonical |
|---|--:|--:|--:|
| split-half ceiling (A vs B, n_ref 15,703) | **0.889** | **0.973** | **0.498** |

**Two independent measurements of the same cells agree on barely half of each other's non-canonical
calls.** So the model's 0.24 non-canonical recall is not a quarter of perfect -- it is roughly half
of what a replicate EXPERIMENT achieves. Panel C states every model number against that denominator.

The ceiling is directional (A-vs-B and B-vs-A swap precision and recall; F1 is identical). The
generator uses **A-vs-B**, whose n_ref (15,703) is the closer match to the model's (16,331).

## Result

| arm | precision | recall | F1 | recall annot | recall non-canon | % ceiling F1 | % annot | % non-canon |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| mamba4 theta=1 | 0.780 | 0.580 | 0.665 | 0.754 | 0.241 | 75% | 78% | 48% |
| mamba4 Poisson | 0.887 | 0.506 | 0.644 | 0.714 | 0.099 | 72% | 73% | 20% |
| attn theta=1 | 0.761 | 0.592 | **0.666** | 0.759 | 0.267 | 75% | 78% | **54%** |
| attn Poisson | **0.918** | 0.491 | 0.640 | 0.705 | 0.074 | 72% | 72% | 15% |

Three readings, all in the generator docstring:

1. **The two-arm framework behaves as designed on real measured data.** Poisson buys precision
   (0.78 -> 0.89 mamba4, 0.76 -> 0.92 attn) and pays recall. This is the strongest independent
   support for the calibration dial, because the reference is measured translation, not a null.
2. **mamba4 and attn are indistinguishable again** (F1 0.665 vs 0.666) -- a third replication after
   the macrophage cross-subtype run and the 4-dataset MS panel.
3. **Canonical recall far exceeds non-canonical recall.** The class the proteogenomics application
   depends on is the one the model recovers worst. Panel B exists so this cannot be read off a
   pooled F1.

## Caveats a caption must carry

- **Never quote panel A's non-canonical recall without the ceiling.** 0.24 against a denominator of
  0.498 is a different claim from 0.24 against 1.0, and only the first is true.
- The fair headline is panel C's: from RNA-seq alone, against 327 M measured footprints, the model
  reaches **~75% of the assay's own self-agreement** on F1 and 77-78% on canonical ORFs.
- **The Poisson arms lose badly on non-canonical recall** (15-20% of ceiling). If the application is
  non-canonical discovery, theta=1 is the arm to use; the calibration dial trades exactly the class
  the application needs. Both arms are always shown (`feedback_two_arm_orf_calling`).
- n_ref is stated on the panel (16,331) per `feedback_orf_call_transcript_space`.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/A4_b721_ground_truth && $PY make_b721_ground_truth.py
```

## B721.221 IS ON-RECIPE -- a 2026-08-15 claim that this panel was off-recipe was WRONG

**Retracted the same day it was written.** An automated recipe audit reported B721 as off-recipe with
"zero filter_tx_heldout runs". That was an artifact of the audit METHOD, not a property of the data.

**The evidence B721 is on-recipe:**

- `proteogenomics/scripts/align_b721_ribo.sbatch` aligns with `--alignEndsType EndToEnd`,
  `--outFilterMultimapNmax 1`, `--quantMode TranscriptomeSAM`, and then calls `filter_tx_heldout.py`.
- All **7 of 7** BAMs in `B721_pilot/riboseq/` carry `@PG` provenance showing they were rewritten by
  that filter (`/data/tmp/emalekos/heldout_filter/<run>_<pid>/clean.bam` -> `sorted.bam`), which is
  exactly its temp layout.

**Why the audit was wrong, because the lesson generalises.** It scraped `logs/*.out` for the filter's
"DONE ->" line. B721's filter ran from a script under `proteogenomics/scripts/`, whose logs are not in
that directory, so it scored zero. **A dataset's absence from a log scrape is not evidence it was
skipped.** The BAM header is the authoritative record where the BAM survives.

The converse also holds and the audit got that wrong too: `mouse_janich_liver_ribo` and
`mouse_gse243134_liver` show NO filter evidence in their `@PG`, yet the logs prove they were filtered
-- their filtered copies live in `data/liver3x3/pool_*/`, while `heldout_bam/` keeps the raw
alignments. Neither logs nor headers alone settle it; the question is always which artifact fed which
pack.

**Nothing about this panel's numbers changes.** No caveat is needed on A4 for recipe provenance.

**The one thing that IS unusual here, and is by design:** `data/packed_heldout_human_b721` is
COVERAGE-ONLY -- `target_counts` is deliberately 0, because the pack is the Ribo-free input for
standalone prediction. The 327 M measured footprints are the separate REFERENCE this panel scores
against. Do not read that 0 as "no measurement".
