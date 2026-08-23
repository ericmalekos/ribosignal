# P9 -- mouse-liver ORF-call UpSet: data inputs

Poster panel. Where does the model's disagreement with observed Ribo-seq actually live? The F1 tables
cannot separate "missed an ORF all three experiments agree on" from "called an ORF no experiment
saw"; this does.

## Eight outputs (2 model families x 4 ORF-class views)

| file | models drawn | ORF classes |
|---|---|---|
| `P9_upset_mouse3x3.{pdf,png}` | mamba4 | all |
| `P9_upset_mouse3x3_cds.*` | mamba4 | `annotated` |
| `P9_upset_mouse3x3_lncrna.*` | mamba4 | `novel` |
| `P9_upset_mouse3x3_uorf.*` | mamba4 | `uORF` |
| `P9_upset_mouse3x3_2models.*` | attn + mamba4 | all |
| `P9_upset_mouse3x3_2models_cds.*` | attn + mamba4 | `annotated` |
| `P9_upset_mouse3x3_2models_lncrna.*` | attn + mamba4 | `novel` |
| `P9_upset_mouse3x3_2models_uorf.*` | attn + mamba4 | `uORF` |

Each writes a sibling `P9_values*.json` with every intersection, including those not drawn.

## Set -> data

| set | source | note |
|---|---|---|
| Janich | `results/mouse_liver_3x3_canon/mamba4_ribo-janich_rna-janich/real_collapsed.txt` | 5 libraries |
| GSE243134 | `.../mamba4_ribo-gse243134_rna-gse243134/real_collapsed.txt` | 21 libraries |
| Wang | `.../mamba4_ribo-wang_rna-wang/real_collapsed.txt` | 2 libraries |
| mamba4 | `.../mamba4_ribo-gse243134_rna-gse243134/pred_preddepth_collapsed.txt` | standalone arm |
| attn | `.../attn_ribo-gse243134_rna-gse243134/pred_preddepth_collapsed.txt` | standalone arm |
| loader / universe | `.../mamba4_ribo-gse243134_rna-gse243134/pred_profiles.npz` | 22,974 tx |
| tx -> gene | `data/tx2biotype_mouse.tsv` | mouse vM38 |
| **FP rates on the plot** | `results/merged_liver_ribocode/overcall/overcall_meta.json` + `overcall_recovery_by_class.tsv` | `scripts/score_merged_liver_overcall.py` |

Filtering is `compare_dropin_calls.build_loader` (pval <= 0.05, ORF >= 90 nt, predicted enrichment
>= 0.5x uniform). Keying is genomic `(gene_id, ORF_gstop)`, so a different representative isoform of
one ORF is not scored as a disagreement (`feedback_orf_call_transcript_space`).

## The "model only" bar is now a MEASUREMENT, not a label

It used to read "candidate FPs". A joint RiboCode call over all 28 aligned with the final recipe mouse-liver
libraries (`results/merged_liver_ribocode`, 26,388 raw calls, +77% to +102% over any single dataset)
measured how many are corroborated, against a control of observed calls that are equally unsupported
by the rest of the 3x3:

| stratum | n | corroborated |
|---|--:|--:|
| model calls that ARE experiment-supported (anchor) | 11,441 | **99.3%** |
| experiment-unique observed calls (**control**) | 1,054 | **74.5%** |
| model-unique, mamba4 | **1,942** | **7.7%** |
| model-unique, both models | 1,504 | 8.7% |

**The annotation must always carry the control alongside the rate.** "8% corroborated" is
meaningless without "vs 74% for a genuinely real, singly-observed ORF"; the generator builds the
label so the two cannot be separated.

**Two architectures agreeing is NOT evidence.** The 5-set family was built expecting the shared bar
to identify real-but-shallow ORFs. It does not: 8.7% is indistinguishable from either model alone.
Read that bar as a shared systematic error, never as a confidence filter. The code comment beside
`FAMILIES` records the falsified hypothesis so it is not re-derived.

## Per-class FP rates drawn on each view

| view | mamba4 | both models |
|---|--:|--:|
| all classes | 92% | 91% |
| annotated CDS | 92% | 91% |
| lncRNA (`novel`) | 92% | 92% |
| uORF | 88% | 86% |

`internal` is the one class where the CONTROL is itself weak (33-54%), so a low model number there is
genuinely ambiguous. It has no view of its own; its counts are in the values JSON.

## Caveats a caption must carry

- `novel` and `lncRNA` are the SAME set here, not nested: 100% of `novel` calls sit on lncRNA
  transcripts in every set checked (954/954 observed GSE243134, 2,774/2,774 model, 649/649 Wang).
  RiboCode calls an ORF `novel` when its transcript has no annotated CDS, and in a pc + lncRNA
  universe those transcripts are lncRNAs.
- `Overlap_uORF` is deliberately EXCLUDED from the uORF view (it overlaps the annotated CDS and
  behaves differently). Its counts are printed in the run summary, not plotted.
- `pred_preddepth` differs per RNA input, so there is no single "the model" set. These use the
  GSE243134 RNA (deepest, 19 libraries); the alternatives are in the values JSON.
- At most 20 intersections are drawn; the omitted tail is annotated on the plot and listed in full
  in the values JSON.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
# FP rates come from the scorer, so run it first if the merged reference changed:
$PY scripts/score_merged_liver_overcall.py
cd figures/P9_upset_mouse3x3 && $PY make_upset_mouse3x3.py
```

If the scorer has not run, the panel falls back to the old "(candidate FPs)" label rather than
plotting a stale number.
