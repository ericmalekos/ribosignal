# E2 — model vs the between-experiment ceiling, by ORF class

Generator: `make_class_ceiling.py` (cas12a env). Outputs `E2_class_ceiling.{pdf,png}`,
`E2_values.{json,tsv}`.

## Sources

| series | file | selection |
|---|---|---|
| ceiling band (shaded) + median line | `results/mouse_liver_3x3/scored/ribo_vs_ribo_by_class.tsv` | all 6 ordered experiment pairs per class; band = min to max F1 |
| model points | `results/mouse_liver_3x3/scored/model_vs_ribo_by_class.tsv` | all 18 cells per arm (2 models x 9 RNA/Ribo combinations); large dot = mean |

Both produced by `scripts/score_liver3x3.py` via `compare_dropin_calls.build_loader`
(pval<=0.05, ORF>=90 nt, predicted enrichment>=0.5x uniform), genomic key `(gene_id, ORF_gstop)`.

## Numbers

| class | ceiling F1 | model obsdepth (mean) | model preddepth (mean) |
|---|--:|--:|--:|
| annotated | 0.992-0.995 | 0.988 | 0.982 |
| uORF | 0.693-0.840 | 0.579 | 0.546 |
| novel | 0.690-0.796 | 0.622 | 0.426 |
| Overlap_uORF | 0.668-0.793 | 0.456 | 0.453 |
| dORF | 0.516-0.635 | 0.134 | 0.115 |
| internal | 0.539-0.620 | 0.022 | 0.020 |

## Reading it

- **Annotated CDS is at the experimental ceiling.** The model points sit inside or against the band.
- **Every other class is below it**, and the deficit is class-specific rather than uniform.
- **`internal` and `dORF` are the real failures**: the band sits at 0.52-0.64 (two labs DO reproduce
  these) while the model is at 0.02 and 0.13. An overall-F1 plot hides this completely, which is why
  this figure exists in this form.

## Caveats

- `Overlap_dORF` is omitted from the panel: n_ref = 47 in the smallest reference, too few for a
  stable F1. It remains in the source TSV.
- Model points pool matched and mismatched RNA inputs. That is defensible here only because RNA
  provenance moves F1 by +0.0030 on average (figure E1); it would not be defensible otherwise.
- Mouse liver panel only. Human internal ORFs score 0.154-0.241, well above the near-zero here, so
  the internal-ORF result must NOT be generalised beyond this panel.

## Regenerating

```
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
  figures/E2_class_ceiling/make_class_ceiling.py
```
