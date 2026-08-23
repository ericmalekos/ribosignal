# E2 — model vs the between-experiment ceiling, by ORF class

Generator: `make_class_ceiling.py` (cas12a env). Outputs `E2_class_ceiling.{pdf,png}`,
`E2_values.{json,tsv}`.

## Sources

| series | file | selection |
|---|---|---|
| ceiling band (shaded) + median line | `results/mouse_liver_3x3_canon/scored/ribo_vs_ribo_by_class.tsv` | all 6 ordered experiment pairs per class; band = min to max F1 |
| model points | `results/mouse_liver_3x3_canon/scored/model_vs_ribo_by_class.tsv` | all 18 cells per arm (2 models x 9 RNA/Ribo combinations); large dot = mean |

Both produced by `scripts/score_liver3x3.py` via `compare_dropin_calls.build_loader`
(pval<=0.05, ORF>=90 nt, predicted enrichment>=0.5x uniform), genomic key `(gene_id, ORF_gstop)`.

## Numbers

**These are the CANONICAL numbers, and they are read from `E2_values.json`.** Until 2026-08-15 this
table still held the OFF-RECIPE values from before the final-recipe rebuild described at the bottom of
this file: the substrate change was documented while the table beneath it was not updated. Every row
differed (e.g. uORF ceiling 0.693-0.840 -> 0.682-0.832, annotated preddepth 0.982 -> 0.975). If you
edit this table by hand again, re-derive it from the JSON.

Pooled over both architectures (the historical default, 18 cells per arm):

| class | ceiling F1 | model obsdepth (mean) | model preddepth (mean) |
|---|--:|--:|--:|
| annotated | 0.9912-0.9948 | 0.9858 | 0.9754 |
| uORF | 0.6817-0.8321 | 0.5468 | 0.5100 |
| novel | 0.6804-0.7848 | 0.5756 | 0.3447 |
| Overlap_uORF | 0.6655-0.7707 | 0.4304 | 0.4286 |
| dORF | 0.4802-0.6058 | 0.1028 | 0.0894 |
| internal | 0.5116-0.6174 | 0.0222 | 0.0197 |

**Per architecture (9 cells per arm each), `pred_preddepth`.** The pooled column above is a mean over
BOTH architectures; nothing in the old JSON said so, so it read as a single-model result. Use these
for a single-architecture figure:

| class | ceiling F1 | attn | mamba4 |
|---|--:|--:|--:|
| annotated | 0.9912-0.9948 | **0.9780** | 0.9729 |
| uORF | 0.6817-0.8321 | **0.5194** | 0.5007 |
| novel | 0.6804-0.7848 | 0.3444 | **0.3450** |
| Overlap_uORF | 0.6655-0.7707 | **0.4353** | 0.4218 |
| dORF | 0.4802-0.6058 | **0.1064** | 0.0724 |
| internal | 0.5116-0.6174 | 0.0130 | **0.0263** |

Render a single-architecture panel with `make_class_ceiling.py --model attn`, which writes
`E2_class_ceiling_attn.*` and `E2_values_attn.*` without touching the pooled outputs. Per-model means
are present in **every** values file (`model_<arm>_mean_<attn|mamba4>`) regardless of `--model`, and
the JSON now carries `model`, `n_cells_per_arm` and an `aggregation` string.

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
  provenance moves F1 by +0.0041 on average for attn (+0.0028 mamba4, +0.0035 pooled -- figure E1,
  corrected 2026-08-15 from a stale +0.0030); it would not be defensible otherwise.
- Mouse liver panel only. Human internal ORFs score 0.154-0.241, well above the near-zero here, so
  the internal-ORF result must NOT be generalised beyond this panel.

## Regenerating

```
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
  figures/E2_class_ceiling/make_class_ceiling.py
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

## About a THIRD of this ceiling's WIDTH is sequencing depth, not protocol (measured 2026-08-16)

The band here is min-max F1 over the 6 ordered experiment pairs at their NATURAL depths, and those
depths differ a lot: wang 78.3M < janich 110.0M < gse243134 150.5M P-sites. Task #75 re-ran the whole
observed-vs-observed matrix with all three thinned to a common 75M (3 seeds; seed sd 0.002-0.017):

| class | natural band | matched-depth band | width |
|---|--:|--:|--:|
| annotated | 0.991-0.995 | 0.990-0.993 | 0.004 -> 0.003 |
| uORF | 0.682-0.832 | 0.720-0.812 | **0.150 -> 0.092** |
| novel | 0.680-0.785 | 0.702-0.788 | 0.104 -> 0.087 |
| Overlap_uORF | 0.665-0.771 | 0.660-0.757 | 0.105 -> 0.097 |
| dORF | 0.480-0.606 | 0.478-0.553 | **0.126 -> 0.075** |
| internal | 0.512-0.617 | 0.500-0.573 | 0.106 -> 0.073 |

**No conclusion in this panel changes.** The band's POSITION barely moves, so annotated still sits at
its ceiling, dORF and internal still fail badly, and the model's deficit is unchanged. What changes is
how the band should be DESCRIBED: it is narrower than the natural-depth version implies, and roughly a
third of its width was sequencing depth rather than between-lab variability.

Do not restate this as "the experiments only disagree because of depth" -- depth explains about half
the disagreement. At matched depth wang still calls ~5% fewer ORFs than the other two, and that
residual is the genuine protocol difference.

Data: `results/liver3x3_matched_depth/e2_ceiling_matched_depth.json` (per-seed bands included).
