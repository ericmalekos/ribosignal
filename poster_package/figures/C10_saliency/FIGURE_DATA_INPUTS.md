# FIGURE C10 -- input saliency: what sequence the model reads: data inputs

**What it shows:** gradient of the predicted in-frame profile mass over an ORF, with respect to the one-hot
SEQUENCE input, per nucleotide. This is "what sequence the model reads to decide to place periodic
ribosome signal here." The start codon is among the most salient positions in the entire transcript and
saliency is enriched inside the ORF -- the model reads the start context + reading frame, it does not just
copy RNA-seq coverage. Panel per exemplar: full-transcript saliency (ORF shaded) + a start-codon zoom
colored by frame.

**Generator:** `make_saliency.py` (cas12a, CPU torch 2.x). Loads the deployed model and backprops on 3
exemplar transcripts. Reusable; regenerates if the model is retrained.

## Data source
- Model: `results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/best.pt` (deployed one-hot orf_v2_attn,
  held-out Hepatocytes -- these transcripts were NOT in training).
- Inputs assembled by the project's `PackedStore` + `RiboDataset` from `data/packed_Hepatocytes/`
  (coverage, target) + the shared 5-channel non-AUG ORF track `data/packed/orf_track_v2.npy` (set via
  `RIBO_ORF_TRACK`; the LOTO tissue packs share one 36,668-tx universe + offsets, verified identical, so
  the Fibroblast-pack ORF track aligns to Hepatocytes).
- Exemplars (same as prediction_examples) + ORF windows from Hepatocytes RiboCode calls:
  POLR1D uORF ENST00000399697.7 (148-183), RPL32 CDS ENST00000429711.7 (78-485),
  LINC02693 novel lncRNA ORF ENST00000647872.1 (259-924).

## Method
- Target = sum over the ORF's frame-0 positions (tstart, tstart+3, ...) of log-softmax(profile logits).
- saliency_seq[i] = sum_c |d target / d feats[i, c]| over the 4 one-hot channels (input-gradient magnitude).
- Coverage-channel saliency is also captured in the JSON (feats[:, -1]) but not plotted.

## Numbers (start-codon prominence)
| exemplar | start-codon saliency rank (1 = highest in tx) | mean-in-ORF / mean saliency |
|---|---|---|
| POLR1D uORF | 8 / 2036 (top 0.4%) | 13.6x |
| RPL32 CDS | 137 / 2094 (top 6.5%) | 3.6x |
| LINC02693 lncRNA | 15 / 7914 (top 0.2%) | 5.1x |

## Caveats
- Curated GOOD exemplars (the same three as prediction_examples), chosen to illustrate; not a random sample.
- Saliency is a first-order (input-gradient) attribution of ONE target (in-frame ORF mass); attention-weight
  maps or integrated gradients would give complementary views. The start-codon prominence is robust across
  all three ORF classes here (uORF, CDS, lncRNA), which is the point.
- The RPL32 CDS start ranks lower (137th) because a long highly-periodic CDS spreads saliency across many
  in-frame body positions; the uORF and lncRNA (shorter / sparser context) concentrate more sharply on the
  start.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/C10_saliency
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_saliency.py
```
