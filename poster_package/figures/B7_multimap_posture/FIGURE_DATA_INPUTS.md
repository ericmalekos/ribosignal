# FIGURE B7 -- RNA-seq multimap posture does not change the model input: data inputs

**What it shows:** the deployed RNA-seq alignment uses STAR `--outFilterMultimapNmax 20` (mm20). A skeptic
asks whether single-mapper (mm1) coverage would change the model. Fast proxy: align one 30M-read DoHH2
RNA-seq subsample both ways, build per-nt transcript coverage each way, compare per-transcript on the DoHH2
expressed universe. Result: identical on ~95% of transcripts (85.7% at Pearson >= 0.999, median 1.0000);
the divergent 5.3% are paralog/repeat/multimap-heavy loci (mm20 inflates depth 7-145x). The model is
per-transcript-normalized and sees identical input on the vast majority of transcripts -> a full mm1
retrain is unwarranted. Panel a = per-tx Pearson distribution; panel b = divergent tail vs multimap
inflation.

**Generator:** `make_multimap_posture.py` (cas12a). Reads the per-tx TSV produced by the proxy compare.

## Data source / how to regenerate the inputs
- Proxy aligns (mm1 + mm20) -> per-nt coverage hd5:
  `proteogenomics/scripts/align_proxy_mm.sbatch` (export `MULTIMAP=1` and `MULTIMAP=20`; same 30M DoHH2
  reads) -> `proteogenomics/data/mm_proxy/DoHH2_mm{1,20}_coverage.hd5`.
- Per-transcript comparison + the TSV this figure reads:
  `proteogenomics/scripts/compare_coverage_mm.py` -> `proteogenomics/data/mm_proxy/coverage_mm_pertx.tsv`
  (tx_id, pearson, depth_mm1, depth_mm20, ratio_mm1_mm20) + `coverage_mm_comparison.txt` (summary line).

## Numbers (DoHH2 expressed universe, 30,950 tx with mm20 depth >= 50)
- median per-tx Pearson(mm1, mm20) = **1.0000**; 85.7% of tx >= 0.999; only 5.3% < 0.95.
- total depth ratio mm1/mm20 = 0.9236 (mm20 has ~8% more reads, all from multimappers).
- 983 tx (3.2%) are multimap-inflated (mm20 > 1.5x mm1) -- the divergent tail (histone clusters, paralogs,
  repeat-embedded transcripts).

## Why this closes the mm1 question (context)
- Deliberate posture: `--outFilterMultimapNmax 20` is for RNA-seq COVERAGE input (matches training). The
  single-mapper rule (`--outFilterMultimapNmax 1`) is a Ribo-seq-only rule (P-site periodicity), not an
  RNA-seq-coverage rule.
- Because coverage is identical on ~95% of transcripts and the model per-transcript-normalizes, retraining
  on mm1 coverage would not materially move outputs. The full Chothani mm1 re-download + retrain was
  therefore stopped (ENA was throttling it to ~2 days anyway).
- Policy consequence (user directive 2026-07-21): since the two postures are equivalent, mm1 (unique
  genomic mappers) is now the STANDARD for RNA-seq coverage alignment too -- smaller intermediates, faster
  STAR, one posture shared with Ribo-seq. All 8 coverage scripts switched to `--outFilterMultimapNmax 1`
  (methods.md sec 2). The deployed model (mm20-trained) accepts mm1 input without a retrain per this figure.
  Already-completed datasets keep their mm20 coverage (proxy-equivalent; re-running downstream would be
  churn); new datasets (B721/THP-1/mouse) are mm1.
- Consistent with Task 21's posture-B check (results.md): the ~18x isoform-multimap inflation is absorbed
  by per-transcript normalization; conclusions are posture-invariant.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/B7_multimap_posture
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_multimap_posture.py
```
