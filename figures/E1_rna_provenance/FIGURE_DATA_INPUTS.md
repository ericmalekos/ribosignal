# E1 — the RNA-seq input need not come from the same experiment

Generator: `make_rna_provenance.py` (cas12a env). Outputs `E1_rna_provenance.{pdf,png}`,
`E1_values.{json,tsv}`.

## Sources

| element | file |
|---|---|
| heatmap cells (F1) | `results/mouse_liver_3x3/scored/model_vs_ribo.tsv`, `model == mamba4` |
| caption delta statistics | `results/mouse_liver_3x3/scored/matched_vs_mismatched_rna.tsv` (all 12 rows, both models) |
| RNA sample counts in the row labels | `data/liver3x3/pool_*/provenance` inputs; janich 7, gse243134 19, wang 2 |

`mamba4` is shown as the primary model (locked decision D1b). `attn` is in the source TSV and behaves
the same way; the caption statistics pool BOTH models so the headline number is not model-specific.

## Numbers

Matched minus mismatched: **+0.0030 F1 mean**, range **-0.0057 to +0.0079**, negative in **2 of 12**
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
