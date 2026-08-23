# E1 — the RNA-seq input need not come from the same experiment

Generator: `make_rna_provenance.py` (cas12a env). Outputs `E1_rna_provenance.{pdf,png}`,
`E1_values.{json,tsv}`.

## Sources

| element | file |
|---|---|
| heatmap cells (F1) | `results/mouse_liver_3x3_canon/scored/model_vs_ribo.tsv`, `model == mamba4` |
| caption delta statistics | `results/mouse_liver_3x3_canon/scored/matched_vs_mismatched_rna.tsv` (all 12 rows, both models) |
| RNA sample counts in the row labels | `data/liver3x3/pool_*/provenance` inputs; janich 7, gse243134 19, wang 2 |

`mamba4` is shown as the primary model (locked decision D1b). `attn` is in the source TSV and behaves
the same way; the caption statistics pool BOTH models so the headline number is not model-specific.

## Numbers

**Corrected 2026-08-15, twice over.** The numbers here were the pre-canonical ones (the substrate
change below was recorded while this line was not updated), AND the delta statistics were read
WITHOUT a model filter while the heatmap above them was mamba4-only -- one figure mixing a
single-architecture grid with a two-architecture summary. Both are fixed; the deltas are now
filtered to the rendered model and the counts appear in the values JSON.

| deltas (matched minus mismatched F1) | mean | range | negative |
|---|--:|--:|--:|
| pooled, BOTH architectures (the old, unfiltered read) | +0.0035 | -0.0042 to +0.0070 | 2 of 12 |
| **attn** (`--model attn`) | **+0.0041** | **+0.0030 to +0.0053** | **0 of 6** |
| mamba4 (default) | +0.0028 | -0.0042 to +0.0070 | 2 of 6 |

**Both negatives are mamba4.** For attn, matched RNA never hurts: the effect is positive in all 6
(3 Ribo references x 2 arms). An all-transformer figure can state that directly, and the pooled
"negative in 2 of 12" understated the result for the architecture the poster actually shows.

The effect is small either way -- a few thousandths of an F1 point -- and that smallness is the
panel's point: RNA provenance is not what limits this model.
model x reference x arm combinations.

## Reading it

The result is the ABSENCE of a pattern. If RNA provenance mattered, the red-outlined diagonal would be
the brightest cell in every column. It is not: in the `wang` column the matched cell (0.919) is the
LOWEST of the three, and in the `janich` column it ties.

What does vary is the ROW. Wang RNA (2 samples) is the weakest input against every reference in both
panels. Library quality matters; provenance does not.

## Caveats

- **Per-panel colour scale.** The within-panel spread is ~0.03 F1; a shared 0-1 scale would render
  both panels uniformly flat and unreadable. The scale therefore EXAGGERATES small differences by
  design, which is the correct choice for showing that the diagonal does not stand out, but means
  colour must not be read as an absolute magnitude. Cell values are printed for that reason.
- The two panels have different colour ranges and are not comparable to each other by colour.
- Depth is unmatched across the three RNA arms (7 / 19 / 2 samples), deliberately: unmatched depth is
  closer to real application. A matched-depth arm is tracked separately as task #75.

## Regenerating

```
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
  figures/E1_rna_provenance/make_rna_provenance.py
```

## Substrate change 2026-08-13: rebuilt on final-recipe alignments

This figure previously read `results/mouse_liver_3x3/scored/`, built from Ribo-seq alignments that
lacked `--alignEndsType EndToEnd` and the ncRNA + cross-gene filter (`docs/PIPELINE_POLICY.md`). It
now reads `results/mouse_liver_3x3_canon/scored/`, rebuilt from re-aligned reads: 40 alignments ->
3 pools -> 9 packs -> 18 dumps -> rescore.

Final-recipe pooling removed 3.5-4.9% of P-sites and 3.8-4.3% of observed ORF calls. Every qualitative
conclusion in this panel held. The off-recipe outputs are archived under
`results/_archive_offrecipe_2026_08_13/` rather than deleted, so the old and new numbers can be
compared instead of one silently replacing the other.
