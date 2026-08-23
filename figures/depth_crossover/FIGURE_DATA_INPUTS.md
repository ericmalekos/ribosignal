# FIGURE_DATA_INPUTS: depth_crossover.png

**Question.** At what Ribo-seq sequencing depth does PREDICTING the profile from RNA-seq (no experiment)
recover more of the deep ORF-call truth than MEASURING at that depth? Makes the LOTO "profile-Pearson
dispersion is a target-quality ceiling" claim concrete: below a crossover depth, the model beats the
experiment.

**Model / setting.** `orf_v2_attn_onehot`, Hepatocytes held out (LOTO; the model never saw Hepatocytes).
Ground truth for every comparison: the OFFICIAL deep Hepatocytes RiboCode calls, restricted to the
held-out test-tx genes, matched on the genomic ORF locus `(gene_id, ORF_gstop)` (same key + filters as the
main drop-in: raw pval_combined <= 0.05 after RiboCode's own BH gate, ORF length >= 90 nt, predicted calls
additionally enrichment >= 0.5x uniform).

## Panels

- **Top -- drop-in F1 vs depth.** y = ORF-call F1 of each call set vs the deep official calls.
- **Bottom -- non-canonical recall vs depth.** y = recall of uORF / Overlap_uORF / dORF / Overlap_dORF /
  novel / internal ORFs (the strata that depend on clean periodicity).
- x (both) = absolute test-set P-site count (log scale) = depth fraction f x full-depth. Full-depth test
  set = ~6.30e8 P-sites (33,918 test transcripts, mean ~18,573 P-sites/tx).

## Series and their data sources

| Series | What it is | Source file(s) |
|---|---|---|
| Measurement curve (blue, markers + line, std error bars over seeds) | RiboCode `detectORF` on the REAL Hepatocytes profile binomially thinned to fraction f of depth, vs deep official | `.../dropin/real_collapsed.txt` (f=1) + `.../dropin/subsample/real_depth<f>_seed<s>_collapsed.txt` (f<1), produced by `subsample_depth.sbatch` -> `ribocode_dropin.py --variant real --subsample <f> --seed <s>` |
| pred_preddepth line (red dashed) | Model's fully standalone prediction (predicted shape x count-head predicted depth; uses RNA-seq + sequence, ZERO Ribo-seq), vs deep official. Depth-independent. | `.../dropin/pred_preddepth_collapsed.txt` |
| pred_obsdepth line (purple dotted) | Predicted shape at full REAL depth, vs deep official (upper reference for the prediction approach) | `.../dropin/pred_obsdepth_collapsed.txt` |
| Crossover marker (vertical) | Interpolated depth where the measurement F1 curve crosses the pred_preddepth line | computed in `depth_curve_summary.json` -> `crossover` |

`.../dropin/` = `results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/dropin/`. The predicted collapsed
files come from the already-run baseline drop-in (`dump_pred_profiles.py` -> `ribocode_dropin.py`);
`pred_profiles.npz` in that dir supplies the real per-nt counts that are thinned, the predicted profile
(for the enrichment filter), and the full-depth P-site total for the x-axis.

## Regeneration

```
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
FIG=$NEW/figures/depth_crossover

# 1) thinned-depth RiboCode calls (SLURM array; expand DEPTHS/SEEDS + --array for the full grid)
sbatch $NEW/scripts/subsample_depth.sbatch                         # pilot: 6 depths x seed 0
# full: sbatch --array=0-23 --export=ALL,DEPTHS="0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002",SEEDS="0 1 2" subsample_depth.sbatch

# 2) aggregate -> curve + crossover
$PY $NEW/scripts/subsample_depth_curve.py \
  --dropin_dir $NEW/$RUN/dropin \
  --official /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/ribocode_per_tissue/Hepatocytes/Hepatocytes_collapsed.txt \
  --profiles $NEW/$RUN/dropin/pred_profiles.npz --out $FIG

# 3) figure
$PY $NEW/scripts/plot_depth_crossover.py --curve $FIG/depth_curve.tsv \
  --summary $FIG/depth_curve_summary.json --out $FIG/depth_crossover.png
```

## Caveats

- Thinning is binomial (each footprint kept i.i.d. with prob f) applied to the pooled per-nt P-site
  counts -- it models shallower SEQUENCING of the same library, not fewer biological replicates or a
  different protocol. Protocol/periodicity-quality differences across real datasets are a separate axis.
- The crossover is specific to Hepatocytes depth + this test-tx universe + the genomic-locus F1 operating
  point (90 nt, 0.5x enrichment). Absolute P-site counts are for these 33,918 test transcripts, not a
  whole library; multiply by (whole-transcriptome / test-set) coverage ratio to map to a real library size.
- pred_preddepth is depth-independent by construction (the count head reads RNA-seq coverage, never
  Ribo-seq), so its line is flat; that is the whole point of the comparison.
