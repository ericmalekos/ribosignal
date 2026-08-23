# C12 -- codon occupancy vs the empirical ceiling: data inputs

Left: the model-free ceiling (observed vs observed). Right: the model against whichever ceiling
applies to it.

## Panel -> data

| element | source | produced by |
|---|---|---|
| both panels, every point | `results/codon_occupancy_canon/_compare_asite.json` | `scripts/codon_compare.py --site asite` |
| the per-codon vectors it compares | `results/codon_occupancy_canon/<label>.tsv` (13 dumps) | `scripts/codon_occupancy.py` via `scripts/run_codon_occupancy_all.sh` |
| CDS coordinates used to assign codons | `data/{human_ribocode_annot_primary,mouse_ribocode_annot}/transcripts_cds.txt` | `scripts/build_*_ribocode_annot.sbatch` |

## What the two ceilings are

Observed-vs-observed correlations between INDEPENDENT datasets, no model involved.

| ceiling | n pairs | median r | meaning |
|---|--:|--:|---|
| within species | 31 | 0.944 | codon occupancy is highly reproducible between experiments |
| across species | 36 | 0.269 | and is almost entirely species-specific (tRNA pools) |

Same-dataset pairs are EXCLUDED. They share one observed vector and correlate at 1.000 by
construction, which would inflate the within-species ceiling.

## Which ceiling applies to which model number

| substrate | model r | ceiling used | fraction |
|---|--:|--:|--:|
| mouse (transcripts never seen) | 0.263 | 0.269 across-species | **98%** |
| human (~99.7% seen in training) | 0.548 | 0.944 within-species | **58%** |

A human-trained model predicting MOUSE codon dwell cannot beat how much codon structure transfers
between species at all, so the across-species ceiling is its ceiling. The two percentages answer
different questions and must not be read as a ranking.

## Caveats a caption must carry

- **The human number is memorisation-contaminated.** Those transcripts were ~99.7% present in
  training (see C13). The mouse number is the clean one and is what the argument rests on.
- **The CAR-T half of the human number sits on a pack with ~29% PCR duplicates** -- that pipeline
  deduplicated the genome BAM but built the pack from the undeduplicated transcriptome BAM
  (`docs/PIPELINE_POLICY.md`). Human 0.548 is provisional pending a rebuild; mouse is unaffected.
- The within-species ceiling has 3 low outliers near r=0.345 (visible in the left panel), so its
  full range is 0.344-0.956 while the median is 0.944. The median is plotted; the range belongs in
  the caption.
- n is small: 13 dumps total (4 human, 9 mouse), and the across-species ceiling rests on 36 pairs
  drawn from those.

## Recipe robustness

The mouse side was recomputed after the Ribo-seq alignments were rebuilt canonically, which changed
total P-sites by 3.5-4.9%. Every median came back identical to three decimals; the two sets of codon
vectors differ by up to 0.0385 per codon yet correlate at r=0.99952. Per-transcript normalization to
the CDS mean divides out a roughly uniform depth change, so this result does not depend on which
alignment recipe a reader would have used.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/C12_codon_occupancy && $PY make_codon_occupancy.py
```

Upstream, if the compare JSON is missing:

```bash
bash scripts/run_codon_occupancy_all.sh
$PY scripts/codon_compare.py --dir results/codon_occupancy_canon \
  --out results/codon_occupancy_canon/_compare_asite.json --site asite
```
