# C11 -- ORF-track channel ablation: data inputs

One panel per model; cells are the change in ORF-calling F1 when one ORF-track channel is zeroed.
Red = the model calls that class BETTER without the channel.

## Panel -> data

| element | source | produced by |
|---|---|---|
| every cell (`delta`) | `results/orf_channel_ablation_canon/ablation_scores.tsv` | `scripts/score_channel_ablation.py` |
| the 12 ablated predictions behind that TSV | `results/orf_channel_ablation_canon/<model>_<arm>/pred_{obsdepth,preddepth}_collapsed.txt` | `scripts/orf_channel_ablation.sbatch` (array 0-11) |
| the unablated baseline (`base_f1`) | `results/mouse_liver_3x3_canon/<model>_ribo-gse243134_rna-gse243134/pred_obsdepth_collapsed.txt` | `scripts/liver3x3_04_dump.sbatch` |
| the observed reference the F1 is scored against | `results/mouse_liver_3x3_canon/<model>_ribo-gse243134_rna-gse243134/real_collapsed.txt` | same |
| pack the predictions were made on | `data/packed_heldout_l3x3canon_ribo-gse243134_rna-gse243134` | `scripts/liver3x3_03_packs.sbatch` |

## Which arm is plotted

`pred_obsdepth` (predicted shape scaled by the reference dataset's observed depth). `pred_preddepth`
is scored and stored in `C11_values.json` but not drawn -- both arms are always computed
(`feedback_two_arm_orf_calling`); the heatmap shows one to stay legible, and the signs agree.

## Filtering

Inherited from `compare_dropin_calls.build_loader`, the same code path as every other drop-in number
in this project: genomic keying `(gene_id, ORF_gstop)`, pval <= 0.05, ORF >= 90 nt, restricted to the
model's test transcripts, and on the predicted side a mean-density enrichment >= 0.5x uniform.
`n_ref` is recorded per class in `C11_values.json` (`feedback_orf_call_transcript_space`).

## Substrate

CANONICAL alignments (`scripts/riboseq_align.sbatch`: EndToEnd + mm1 + ncRNA/cross-gene filter), per
`docs/PIPELINE_POLICY.md`. The earlier off-recipe version of this ablation is archived under
`results/orf_channel_ablation/`; sign agreement between the two is 57/60 across 5 classes x 2 models
x 6 arms, with all three disagreements within 0.005 of zero.

## Caveats a caption must carry

- The `annotated` column is ~0.000 everywhere. That is a real result (the ORF track is not what
  drives canonical CDS calling, which is already F1 0.986-0.990) and not missing data.
- `internal` F1 is LOW in absolute terms (base 0.008 attn / 0.025 mamba4). A +0.05 delta is a large
  RELATIVE change on a small number; the figure shows the sign and the dissociation, not a claim that
  internal ORFs become well-called.
- `f012` is not the sum of ch0+ch1+ch2. The frame channels are partly redundant, so removing one
  understates the group; that is why the joint arm exists.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/C11_channel_ablation && $PY make_channel_ablation.py
```

Upstream, if the scores TSV is missing:

```bash
$PY scripts/score_channel_ablation.py \
  --ablation-dir results/orf_channel_ablation_canon \
  --base-dir results/mouse_liver_3x3_canon \
  --out results/orf_channel_ablation_canon/ablation_scores.tsv
```

## "So should the frame channels be removed?" -- no, and the figure alone could mislead

The red `internal` column invites that conclusion. Three things argue against it, and a caption
should carry at least the first.

**1. The internal gain is noise, not recovered ORFs.** In the `f012` arm, internal PRECISION is
**0.035** (attn) / **0.041** (mamba4): roughly 96% of the internal calls the model starts emitting
are wrong. F1 rises only because recall climbs from ~0 to ~0.06 on top of that. The frame channels
are not concealing good internal calls; they are suppressing a class the model cannot call well with
or without them.

**2. The classes that lose are larger and are the ones the application needs.** Removing all three
frame channels costs uORF -0.048/-0.037 and novel -0.038/-0.043 (attn/mamba4), over 697 and 663
reference ORFs, against 207 internal ORFs gaining. About 6.6x more reference ORFs sit on the losing
side, and uORF/novel are exactly the non-canonical classes the proteogenomics application depends
on. Aggregate F1 falls 0.024 on both models.

**3. This is an inference-time knockout, not a retrain.** Every arm zeroes a channel on a model that
was TRAINED with it, so the figure measures what the trained model RELIES on -- not what an
architecture without the channel would achieve. A model trained without frame channels could
redistribute onto sequence and recover the uORF/novel loss, and would never acquire the frame prior
that suppresses internal. The architecture question is genuinely open and this experiment does not
close it.

The clean test is a training arm with a frame-zeroed ORF track: `scripts/train_canon.sbatch` with
`RIBO_ORF_TRACK` pointed at such a track, ~14 h on one GPU. Not run as of 2026-08-13.
