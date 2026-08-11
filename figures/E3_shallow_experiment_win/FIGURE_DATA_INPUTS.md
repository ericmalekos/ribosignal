# E3 — standalone model out-recalls a shallow Ribo-seq experiment

Generator: `make_shallow_experiment_win.py` (cas12a env). Outputs
`E3_shallow_experiment_win.{pdf,png}`, `E3_values.json`, `E3_values.tsv`.

## What each panel reads

### Panel A — the three mouse-liver libraries

| quantity | source |
|---|---|
| pooled P-sites | `data/liver3x3/pool_{janich,gse243134,wang}/target_counts.npy`, summed |
| 3-nt periodicity (mean f0) | `data/liver3x3/pool_*/_work/psites/*_pre_config.txt`, the `f0_percent` column of every read length RiboCode `metaplots` selected, averaged |
| Ribo sample count | `data/heldout_bam/{mouse_janich_liver_ribo,mouse_gse243134_liver,mouse_wang_liver_ribo}/*Log.final.out` |

GSE243134's directory holds both arms, so the Ribo runs are filtered against
`data/liver3x3/gse243134_ribo_runs.txt` (21 of 40). Counting all 40 would inflate its sample count
with RNA runs.

### Panel B — shallow experiment vs standalone model

| series | source |
|---|---|
| Wang's recall (blue) | `results/mouse_liver_3x3/scored/ribo_vs_ribo_by_class.tsv`, rows where `compared == wang` |
| model recall (red) | `results/mouse_liver_3x3/scored/model_vs_ribo_by_class.tsv`, `arm == pred_preddepth`, best cell per (class, reference) |

Both TSVs come from `scripts/score_liver3x3.py`, which loads calls through
`compare_dropin_calls.build_loader` — the same filtering behind every other drop-in number in this
project (pval<=0.05, ORF>=90 nt, predicted enrichment>=0.5x uniform) — and keys ORFs genomically on
`(gene_id, ORF_gstop)`.

## Numbers in the figure

| class | reference | Wang recall | model recall | margin |
|---|---|--:|--:|--:|
| novel | gse243134 | 0.581 | 0.709 | +0.128 |
| novel | janich | 0.644 | 0.724 | +0.081 |
| uORF | gse243134 | 0.560 | 0.676 | +0.117 |
| uORF | janich | 0.626 | 0.688 | +0.061 |

The winning arm is `pred_preddepth` in all four, and the winning model is `attn` in all four.

## Caveats that must travel with this figure

- **The model does not beat deep Ribo-seq.** It loses to both deeper experiments on recall, and to
  all three on precision and therefore F1 (novel F1 0.595-0.656 and uORF 0.436-0.687 against
  experimental 0.690-0.796 and 0.693-0.840). The claim is bounded: better than a *shallow*
  experiment, at no experimental cost.
- **Every one of the 28 recall wins across novel and uORF is against Wang.** None is against Janich
  or GSE243134. That is stated in the caption rather than left for a reader to discover.
- **Wang is shown to be typical, not chosen to be weak.** Panel A exists for that purpose; without it
  the comparison would be unfalsifiable.
- **Internal ORFs are excluded from this figure and must not be implied by it.** The model emits
  11-38 internal calls against 147-223 in the reference and hits 0-4, i.e. F1 0.000-0.037 against a
  between-experiment ceiling of 0.539-0.620. See the `Model vs observed, by ORF class` table in the
  tutorial.
- **Not comparable to pre-2026-08-10 Janich numbers.** RiboCode `metaplots` re-derives periodic read
  lengths per sample, so the rebuilt Janich pack differs ~1.7% from its historical one; the
  superseded results are archived under `results/_archive_pre_2026_08_10_janich_pack/`.

## Regenerating

```
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
  figures/E3_shallow_experiment_win/make_shallow_experiment_win.py
```

Everything is read at build time from `results/mouse_liver_3x3/scored/` and `data/liver3x3/`, so the
figure follows the factorial if it is re-run. No numbers are transcribed into the script.
