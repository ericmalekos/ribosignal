# Results: Ribo-seq signal prediction model

> **Which numbers are live?** -> **`docs/STATUS_CURRENT_VS_ARCHIVED.md`** (reviewed 2026-08-14).
> One page listing CURRENT vs SUPERSEDED vs IN-PROGRESS for models, ORF-call results, proteogenomics
> and figures, plus the claims that were measured and found FALSE. Read it before quoting anything
> from an older section of this file -- several entries below have been superseded in place and say
> so, but the status page is the index.


Status (updated 2026-07-16): model TRAINED + FULLY EVALUATED across Tasks 1-20. Deployment
recommendation = one-hot `orf_v2_attn` (~5M params, CPU-runnable; FM embeddings give no lift over
one-hot per Task 15). Landed: fold-0 + multi-fold held-out-chromosome baselines (Tasks 6, 9);
RiNALMo-vs-Orthrus and input/architecture ablations (Tasks 7, 8, 15, 17); uORF/dORF-aware eval
(Task 9c); 9-fold leave-one-tissue-out cross-tissue transfer (Task 12); localization + RiboCode
drop-in ORF-calling (Tasks 13, 14); cross-study human Ruiz-Orera + cross-species mouse Wang transfer
(Task 16); depth crossover, predict-vs-measure (Task 18); non-AUG ATG+CTG drop-in (Task 19); and the
Kozak start-context ablation (Task 20), whose finding -- the hand-picked Kozak factor is redundant
with what the sequence backbone learns -- was acted on by flipping the ORF-track default to
`--kozak none`. In progress: posture-B multimap sensitivity check (Task 21). See methods.md for the
per-task methodology and the auto-memory `project_riboseq_signal_model.md` for the running state.

## Task 1: Ribo-seq per-sample BAM cleaning posture (verified)

Verdict: the correct base for the per-nt P-site target is the 32 RiboCode
`*_psites.hd5`, which were built from the CLEANED `Aligned.toTranscriptome.out.bam`.
The three per-sample BAM flavors are in different states:

| per-sample file | posture | use for target? |
|---|---|---|
| `{SRR}.Aligned.sortedByCoord.out.bam` (genomic) | CLEANED (NH==1 + rRNA/tRNA/miRNA loci dropped) | genomic, not the transcript axis |
| `{SRR}.Aligned.toTranscriptome.out.bam` (transcriptome) | CLEANED (ncRNA-tx + cross-gene multimappers dropped; within-gene isoform multimappers retained) | YES (source of the hd5) |
| `{SRR}.toTranscriptome.sorted.bam` (transcriptome) | RAW / dirty (120M+ records, NH up to 1000+) | NO, do not use |

Provenance: Fibroblast hd5 timestamped Jun 28 00:13+, after the in-place clean of
`Aligned.toTranscriptome.out.bam` at Jun 27 17:42. STAR aligned with
`--outFilterMultimapNmax 20`; NH==1 filtering was a separate Phase-30 step on the
genomic BAMs. Target multimapper posture = A (reuse hd5 as-is, user-confirmed).

### Empirical verification (`scripts/verify_target_inputs.py`, logs/verify_target_inputs.txt)

```
input hd5 files            : 32
axis transcripts           : 509,650
per-sample P-site depth    : min 2,144,083  max 7,676,156  total 121,618,960

-- ncRNA drop-list posture --
drop-list transcripts      : 7,587
  present on axis          : 7,587   (present-but-zero, as expected)
  absent from axis         : 0
present ncRNA by gene_type : misc_RNA=2207, snRNA=1901, miRNA=1879, snoRNA=942,
                             rRNA_pseudogene=497, scaRNA=49, rRNA=47, unknown=24,
                             Mt_tRNA=22, ribozyme=8, sRNA=5, vault_RNA=4, Mt_rRNA=2
pooled P-sites over present ncRNA : sum=0  max_single_tx=0  nonzero_tx=0/7587

-- positive controls (housekeeping) --
ACTB  : isoforms=43  pooled_total=38,724,541  top=ENST00000462494.5 (1,234,486)  len_match=43/43
GAPDH : isoforms=40  pooled_total= 4,947,235  top=ENST00000229239.10 (162,424)    len_match=40/40
```

Interpretation: all 7,587 rRNA/tRNA/miRNA/sno/sn drop-list transcripts present on the
axis pool to exactly 0 P-sites (reads over them were removed upstream), while ACTB and
GAPDH pool to millions (the reader sees real signal, so the 0 is real, not an all-zero
bug), and every control transcript's hd5 per-nt array length equals its GENCODE mature
length (the per-nt axis is the transcript-nt axis the embeddings and RNAseq align to).

## Task 2: Fibroblast per-nt P-site target (transcriptome coords) -- BUILT

Build: `scripts/build_fibroblast_psite_target.py`, SLURM job 35119831 on phoenix-21,
18:10:16 to 18:17:08 (about 7 min), pooling the 32 Fibroblast `*_psites.hd5`
element-wise (posture A).

Outputs (`data/target/`):

| file | size | contents |
|---|---|---|
| `Fibroblast_psites_pooled.hd5` | 3.75 GB | `transcript_ids` (vlen-str) + `p_sites` (vlen-int32, per-nt pooled, gzip-4); attrs tissue, n_samples=32, n_transcripts=509650, psites_number_total=121,618,960, psites_number_per_sample (JSON), multimapper_posture=A |
| `Fibroblast_psites_summary.tsv` | 39.7 MB | one row/tx: tx_id, length, total_psites, n_nonzero_nt, max_psite, gene_id, gene_name, chrom, transcript_type |

Headline numbers:

- grand total pooled P-sites = 2,205,261,175 ; max single pooled nucleotide = 191,722
  (int32-safe).
- transcripts with nonzero signal = 343,360 / 509,650 (67.4%).
- signal-bearing transcripts by transcript_type: protein_coding 193,953 ; lncRNA
  71,362 ; retained_intron 30,429 ; protein_coding_CDS_not_defined 21,944 ;
  nonsense_mediated_decay 19,676 ; processed_pseudogene 4,360 ; (others smaller).
- top transcripts by pooled P-sites: THBS1 (2.80M) and FN1 (2.25M), both
  protein_coding, which are the canonical fibroblast-secreted ECM genes
  (thrombospondin-1, fibronectin-1) -- a biological sanity check that the deepest Ribo-seq
  signal is the genes a fibroblast translates most.

### Multimapping multiplicity (a property of posture A, documented)

The pooled per-nt total (2,205,261,175) is 18.1x the pooled distinct footprint count
(121,618,960). That factor is the mean within-gene isoform multiplicity: each footprint
is counted on every compatible isoform (posture A). It inflates absolute magnitude and,
because a read over a shared exon maps to more isoforms than a read over an
isoform-unique exon, can also distort the relative per-nt shape between shared and unique
regions. This is expected and acceptable for the first target: the model is scored with
scale-invariant per-transcript Pearson/Spearman + 3-nt periodicity, and per-transcript
normalization removes the magnitude inflation. It motivates (a) per-transcript
normalization at train/eval time and (b) the pre-registered posture-B (primary-only)
sensitivity check.

### Validations

- Build correctness: `pooled[tx] == sum over the 32 sources[tx]` element-wise
  (`np.array_equal` True) for FN1 (ENST00000354785.11), ACTB (ENST00000462494.5),
  GAPDH (ENST00000229239.10), and a drop-list ncRNA (ENST00000290239.7, 0==0). No
  double-counting.
- Axis cross-check: for all 443 RiNALMo per-token pilot transcripts, the pooled hd5
  per-nt array length equals the embedding row count L (443/443 match, 0 mismatch,
  0 off-axis). The target aligns 1:1 with the FM embeddings, position-for-position.

## Task 3: Fibroblast RNAseq per-nt coverage (input feature) -- DONE

Second model input: per-nt RNAseq read coverage on the same transcriptome axis as the
target. Strandedness ISR (reverse), posture A to match the target, sense filter
`is_read2 != is_reverse`.

Smoke (SRR15513233): STAR 55.4M read pairs, 91.79% uniquely mapped + 3.08% multi
(~95% mapped, adapter-free as expected). The per-nt counter walked 1,406,801,328
alignment records and counted 1,387,087,336 sense records (mean within-gene isoform
multiplicity ~12.7x per mate under posture A; the ~1.4% dropped are antisense mates,
confirming the strand filter), covering 341,305/509,650 transcripts (matches the
target's 343,360). Output hd5 3.75 GB.

Counter decision: pysam single-threaded (`scripts/rnaseq_coverage.py`), ~81 min/sample.
deepTools bamCoverage was benchmarked and rejected (>3.5 h/sample: 509,650-contig
`--binSize 1` coverage over a 1.4B-record BAM is pathological; see methods.md section 4).

Fleet: SLURM array 35123909 (`align_rnaseq_array.sbatch`, samples 234-264; 233 reused
from the smoke), raised to 32-wide. All 32 per-sample coverage hd5 produced (uniform
3.5 GB each; validated one: SRR15513236 counted 719.7M sense of 722.9M records, ACTB
total 17.9M over its 2,554-nt length matching the target length). SLURM packed 26 of the
32 tasks onto one node (phoenix-20), so the tail ran longer under I/O contention -- a
lesson to spread with a lower throttle or `--nodes`/`--spread-job` next time. The last
task (SRR15513256) ran 3 h 22 m: it is the deepest sample (1.92B sense records vs the
~0.47-1.67B range across the 32) counted single-threaded on the contended node, i.e.
genuine workload, not a hang. All 32 completed cleanly (0 failed); per-sample sense counts
range 468M-1,925M, every file covering all 509,650 transcripts.

Pooled (`pool_rnaseq_coverage.py`, job 35124893, 6 m 32 s): element-wise sum of the 32
per-sample coverage hd5 -> `Fibroblast_rnaseq_coverage_pooled.hd5` (3.5 GB) +
`Fibroblast_rnaseq_coverage_summary.tsv`. Grand total pooled coverage = 757,573,677,538
per-nt read-depth counts (posture A).

Axis note: the per-sample coverage `transcript_ids` follow the STAR `@SQ` order, which is
the SAME SET as the target's RiboCode order but a DIFFERENT permutation (160,133/509,650
positions differ; 0 set difference). All 32 coverage files share one axis (so positional
pooling is valid among them), but coverage-to-target and coverage-to-embedding joins must
be BY transcript id, never by position. All 36,668 universe transcripts are present in
both axes.

Validation (`scripts/validate_coverage.py`, `logs/validate_pooled_coverage.txt`): pooled
coverage and target share the same 509,650-transcript set; ACTB and GAPDH per-nt lengths
match the target (2,554 and 1,285) with large pooled coverage totals (714M, 488M).
Coupling total-coverage vs total-P-sites (log1p Pearson over the 36,668 universe) = 0.703,
in the expected 0.7-0.8 band -- coverage is a good transcript-level magnitude prior.
Mean-per-nt coverage vs salmon mean TPM = 0.490: the gap from a higher value is the
posture difference (raw multimap-counted coverage under posture A vs salmon's EM-resolved,
decoy-aware TPM), consistent and expected, not a defect.

## Task 4: expressed universe, per-token RiNALMo embeddings, held-out-chrom split -- DONE

Expressed universe (`scripts/define_universe_and_fasta.py`): protein_coding + lncRNA with
Fibroblast salmon meanTPM >= 1 and mature length <= 10,000 nt (RiNALMo cap) = 36,668
transcripts (34,248 pc + 2,420 lncRNA) -> `data/fibroblast_universe.tsv` + `.fa`.

Per-token RiNALMo embeddings (`scripts/extract_rinalmo_universe_array.sbatch`, GPU array
of 8 length-balanced chunks on A5500, ~1 h 05 m per chunk): 36,668/36,668 transcripts
embedded, one `{tx}_tokens.npy` per transcript of shape (L, 1280) float16 with L ==
transcript length, ~234 GB total under `data/rinalmo_token_emb/chunk_{1..8}/`. Indexed by
`scripts/build_embedding_index.py` -> `data/rinalmo_token_emb/tx_index.tsv` (36,668 rows,
0 missing). Axis validation: for 6/6 sampled transcripts the embedding row count L equals
the GENCODE length exactly (1:1 with the target and coverage per-nt axes).

Held-out-chromosome split (`scripts/build_chrom_split.py`): transcripts grouped by
chromosome (gene-disjoint by construction) into 5 count-balanced folds ->
`data/splits/fibroblast_chrom_kfold.json`. Fold 0 (chr1/7/8/22, ~7,473 test tx) is the
canonical first held-out-chromosome test.

## Task 5 (design): data-driven model decisions

From `scripts/inspect_distributions.py` (`logs/inspect_distributions.txt`; 1,500
expressed+translated tx sampled, whole-universe per-tx summaries):

| measurement | value | design consequence |
|---|---|---|
| per-nt zero fraction (target) | 74.2% | sparse count profile |
| per-nt var/mean (target) | ~2,355 | overdispersed -> NB / multinomial, NOT Poisson |
| 3-nt periodicity (autocorr lag 3/6/9/12 vs off-frame) | 0.32/0.24/0.23/0.21 vs ~0.02-0.05 | strong; reproduce it + use as eval metric |
| per-tx P-sites range | median 3,463, max 2.8M | huge -> per-transcript normalization essential |
| coverage per-nt zero fraction | 5.0% (mean 147, max 20k) | dense conditioning input; feed as log1p |
| within-tx coverage-vs-P-site Spearman | median 0.17 | coverage does NOT localize; embedding must |
| across-tx total coverage-vs-P-site (log1p Pearson) | 0.78 | coverage IS a good magnitude prior |

Decision: a BPNet-style dual head over `[emb ; log1p(coverage)]` -- residual dilated 1-D
convolutions (O(L), the FM embedding already carries long-range context) feeding a
**profile head** (log-softmax over positions + multinomial NLL; scale-free, transferable,
captures periodicity + CDS localization) and a **count head** (log-count / NB regression;
magnitude, dataset-depth-dependent). Evaluated by per-transcript Pearson/Spearman +
periodicity recovery on the held-out chromosome, split by biotype. Rationale in
methods.md section 5B.

### Task 5 (implementation): code complete, pipeline chained

Implemented (methods.md section 5C; all `ruff`-clean, py310):

- `scripts/model.py` -- `RiboSignalModel` dual-head dilated CNN + `profile_multinomial_nll`
  + `count_mse`.
- `scripts/dataset.py` -- numpy-only loader over ragged packed target/coverage +
  per-tx embedding `.npy`, token-budget batching, chrom-fold splits.
- `scripts/pack_target_coverage.py` -- ragged 1-D pack of target + coverage for the
  universe (numpy-only reads at train time).
- `scripts/train.py` -- training loop, cosine schedule, early stop, scale-invariant
  per-transcript eval split by biotype.
- `scripts/validate_coverage.py`, `scripts/plot_predictions.py` -- pooled-coverage
  validation and post-training figures.

End-to-end validation on a 1,500-transcript dry-run pack (per-sample coverage stand-in,
tiny model, CPU): loss fell 12.46 -> 7.75, val profile Pearson rose across epochs, and the
model already produced 3-nt periodicity in its predicted profile (period metric ~0.77 vs
observed ~0.21) -- the profile head captures frame structure. Checkpointing, early stopping,
and metric export all verified.

Pipeline chained with SLURM `afterok` dependencies (auto-runs unattended): pool
(35124893) -> post-pool validate+pack (35125105) -> fold-0 GPU training (35125106,
one A5500 in the CUDA SIF). Results land in `results/fold0_rinalmo/`
(`history.json`, `best.pt`, `test_metrics.json`).

## Task 6: fold-0 held-out-chromosome result (RiNALMo) -- FIRST BASELINE

Training (job 35125106, one A5500, 2 h 31 m): early-stopped at epoch 11, best val at
epoch 5 (val profile Pearson 0.595). The whole afterok chain ran unattended after the
coverage fleet finished. Test = fold 0 (chr1/7/8/22), gene-disjoint from train (folds 2-4)
and val (fold 1); 6,962 test transcripts scored (>= 50 total P-sites).

| biotype | n | profile Pearson (median) | Spearman (median) | periodicity pred / obs |
|---|---|---|---|---|
| all | 6,962 | 0.599 | 0.225 | 0.521 / 0.219 |
| protein_coding | 6,698 | **0.604** | 0.241 | 0.526 / 0.223 |
| lncRNA | 264 | 0.343 | -0.171 | 0.144 / 0.041 |

Reading of the numbers:

- The approach works: median per-transcript profile Pearson **0.60 for protein_coding on
  held-out chromosomes**. From (RiNALMo per-token embedding + per-nt RNAseq coverage) the
  model predicts where ribosomes sit along a transcript it has never seen, and the
  prediction generalizes across chromosomes.
- The model learned reading frame: the predicted profiles carry stronger 3-nt periodicity
  than the noisy observed P-site profiles (pc 0.526 predicted vs 0.223 observed) -- it
  effectively denoises toward the ideal periodic footprint pattern rather than copying the
  experimental noise, which is exactly the behaviour a signal/denoising model should show.
- Spearman (0.24) is far below Pearson (0.60) by construction: 74% of nt are zero, so the
  rank correlation is dominated by ranking the huge tie-block of near-zero positions
  (noise), while Pearson is dominated by placing the mass on the real peaks (which the
  model gets right). Pearson is the meaningful profile-shape metric here.
- lncRNA is appropriately weak (Pearson 0.34, observed periodicity only 0.04): most
  lncRNAs are not genuinely translated, so their measured "P-site" signal is largely
  non-periodic noise that no model can predict a coherent profile for. The gap between pc
  and lncRNA is itself a signal that the model captures real translation structure, not a
  generic coverage-shape artefact.

Figure `results/fold0_rinalmo/figures/fold0_predictions.png` (`scripts/plot_predictions.py`)
shows observed vs predicted per-nt profiles with periodicity zooms and the per-transcript
Pearson distribution by biotype.

### 5-fold cross-validation (RiNALMo, both inputs)

All 5 gene-disjoint chromosome folds (`train_allfolds.sbatch` folds 1-4 as array 35126435 +
fold 0 above), same config. The result is tight across folds.

| fold | pc profile Pearson | pc periodicity pred / obs | lncRNA Pearson |
|---|---|---|---|
| 0 | 0.604 | 0.526 / 0.223 | 0.343 |
| 1 | 0.597 | 0.502 / 0.218 | 0.382 |
| 2 | 0.589 | 0.505 / 0.219 | 0.370 |
| 3 | 0.596 | 0.458 / 0.225 | 0.378 |
| 4 | 0.597 | 0.487 / 0.224 | 0.402 |
| **mean +/- sd** | **0.597 +/- 0.005** | 0.496 / 0.222 | 0.375 +/- 0.021 |

The protein_coding profile Pearson is 0.597 +/- 0.005 across five held-out-chromosome folds:
the fold-0 number was not a lucky split. Predicted periodicity (~0.50) consistently exceeds
observed (~0.22), so the denoising-toward-frame behaviour holds fold to fold.

### Input ablations (fold 0)

`--input_mode` zeros one input (array 35126496). Tests whether the FM embedding adds value
beyond RNAseq coverage.

| inputs | pc profile Pearson | reading |
|---|---|---|
| both | 0.604 | baseline |
| embedding-only | 0.593 | embedding carries almost all profile shape (coverage adds +0.011) |
| coverage-only | 0.171 | coverage alone barely localizes; periodicity ~0 (no frame) |

The embedding-only result (0.593, matching the full 0.604) confirms the design: profile
shape -- where ribosomes sit -- is sequence-determined, so the embedding does the
localization. Coverage-only collapses to 0.171 with predicted periodicity -0.013 (no
reading frame at all), so RNAseq coverage carries none of the frame/localization signal.
Coverage's real contribution is magnitude (the count head), which the scale-free profile
Pearson does not measure -- hence coverage barely moves this metric (+0.011). The magnitude
(count-head) accuracy that credits coverage is reported below.

### Magnitude (count-head) accuracy -- a clean double dissociation

`scripts/eval_extra.py` (jobs 35128965 / 35129022) adds the metric the profile Pearson
misses: the per-transcript Pearson between the count head's predicted log-count and the
observed `log1p(total P-sites)` (protein_coding test transcripts). Coverage feeds the count
head, so this is where coverage should earn its keep.

| inputs | profile Pearson (shape) | count Pearson (magnitude) |
|---|---|---|
| both | 0.604 | 0.783 |
| embedding-only | 0.593 | 0.399 |
| coverage-only | 0.171 | 0.789 |

All three rows are fold 0 (n=6698 pc test transcripts) so the input configurations compare on
the same split; the count metric is stable across the 5 folds (0.785 +/- 0.025).

This is a clean double dissociation: each input dominates exactly one axis. **Embedding owns
shape** (embedding-only 0.593 ~= both 0.604, coverage-only collapses to 0.171), **coverage
owns magnitude** (coverage-only 0.789 ~= both 0.783, embedding-only drops to 0.399). The
two-input design is justified on both axes with numbers, not assertion: sequence sets where
and in what frame ribosomes sit, RNAseq abundance sets how much.

### Posture-B sensitivity (single-isoform-gene proxy)

Rather than rebuild the target primary-only, `eval_extra.py` reports profile Pearson split by
whether a test transcript's gene has a single annotated isoform (posture A == posture B by
definition, no multimap ambiguity) vs multiple isoforms (which bear the ~18x multimapping).

| subset | pc profile Pearson (5-fold) | n (per fold) |
|---|---|---|
| single-isoform gene | ~0.55 (0.51-0.58) | 35-64 |
| multi-isoform gene | ~0.597 (0.59-0.60) | 6.2k-6.8k |

The multi-isoform transcripts -- the ones actually affected by posture-A multimapping -- are
not worse than the clean single-isoform transcripts (if anything slightly better, being more
expressed). So the multimap inflation is not distorting the shape prediction, and a full
primary-only (posture-B) target rebuild is not warranted by this evidence (it remains
available as the pre-registered heavy check). Caveat: the single-isoform subset is small
(n=35-64/fold), so this is suggestive, not definitive.

A full code-first writeup of the whole experiment is at `tutorial/index.html`
(self-contained; motivation, all decisions, inputs/outputs, results, reproduce commands).

## Task 7: Orthrus vs RiNALMo embedding comparison (fold 0) -- DONE

Question: does the choice of foundation model matter for the per-nt profile? The same
Fibroblast universe was embedded per-token with Orthrus 4-track (a Mamba mRNA model, 512-d,
sequence-only; see methods.md section 4C) alongside RiNALMo (1280-d), and the identical
dual-head model was trained on RiNALMo alone, Orthrus alone, and the two concatenated
(1792-d) on the same gene-disjoint fold-0 split, seed, and hyperparameters
(`scripts/train_backends.sbatch`, aggregated by `scripts/aggregate_backends.py`;
results/backend_cmp/). `usable()` intersects the per-backend indices so all three train and
test on the identical 6698 pc + 264 lncRNA test transcripts.

| backend | d_emb | pc profile P | pc period_pred | count P | lncRNA profile P |
|---|---|---|---|---|---|
| rinalmo | 1280 | 0.612 | 0.485 | 0.780 | 0.346 |
| orthrus | 512  | 0.617 | 0.498 | 0.810 | 0.369 |
| concat  | 1792 | 0.614 | 0.449 | 0.763 | 0.341 |

(period_obs 0.223 for all; pc profile Spearman ~0.24 for all, the 74%-zeros effect.)

- Orthrus (mRNA-specialized Mamba, 512-d) matches and marginally beats RiNALMo (general RNA
  LM, 1280-d) on every axis -- pc profile, periodicity recovery, count magnitude, and lncRNA
  profile -- at 40% the embedding dimension. The pc gap (+0.005) is within run-to-run noise,
  but the edge is consistent across all four metrics and largest on the depth-relevant count
  head (0.810 vs 0.780) and on lncRNA (0.369 vs 0.346).
- Concatenating the two does NOT help: concat lands between the singles on pc (0.614) and is
  worst on lncRNA (0.341) and count (0.763). The two embeddings are largely redundant for
  this task, so stacking them adds parameters without information and mildly hurts the small
  lncRNA set.
- Practical read: a single Orthrus 4-track embedding is the efficient default (best signal at
  the lowest dimension); RiNALMo is an equivalent alternative; the concatenation is not worth
  its cost. One fold only, so a multi-fold confirmation is the natural follow-up.

## Task 8: input and architecture improvements (fold 0) -- DONE

Two goal-driven additions were tested against the RiNALMo fold-0 baseline (0.612 pc / 0.346
lncRNA), same split/seed/hyperparameters, only the input track and/or architecture changing
(`scripts/train_improve*.sbatch`, aggregated by `scripts/aggregate_improve.py`;
results/improve/):

- ORF-candidate track (`scripts/build_orf_track.py`): a sequence-derived, annotation-free
  per-nt track marking ATG..in-frame-stop occupancy in all 3 frames + start + stop. Two
  start-channel variants -- ATG-only (v1) and non-AUG-aware (v2: graded start_ext over AUG +
  the 9 near-cognate starts 1 nt from AUG, Kozak-weighted). Built annotation-free so it works
  for uORFs, dORFs, and lncRNA ORFs, not just canonical CDS.
- Attention: 2 self-attention layers after the dilated-conv body, giving full-transcript
  context past the CNN's ~4 kb receptive field (`--n_attn_layers`).

| config | ORF track | attn | pc profile P | lncRNA profile P |
|---|---|---|---|---|
| baseline    | none   | 0 | 0.612 | 0.346 |
| orf         | ATG    | 0 | 0.619 | 0.384 |
| orf_v2      | nonAUG | 0 | 0.618 | 0.334 |
| attn        | none   | 2 | 0.617 | 0.321 |
| orf_attn    | ATG    | 2 | 0.608 | 0.405 |
| orf_v2_attn | nonAUG | 2 | 0.614 | 0.362 |

- The ATG ORF track is the best single change: pc 0.619 (+0.007) and lncRNA 0.384 (+0.038). It
  supplies candidate ORF structure the FM embedding only approximates, and helps the small
  lncRNA set most.
- Non-AUG start information does NOT help and actively hurts lncRNA: orf_v2 lncRNA 0.334 (below
  even baseline 0.346) vs orf-ATG 0.384; with attention 0.362 vs 0.405. The graded near-cognate
  start channel adds noise rather than signal on this metric, worst on the class it was meant to
  help. Caveat: profile Pearson is a whole-transcript shape metric dominated by the main ORF, so
  it is not the right probe for non-AUG ORF detection -- the proper test is localization to a
  held-out non-canonical ORF truth set (not yet built), and the near-cognate weights are a
  first-pass ordinal.
- Attention is a pc-vs-lncRNA trade-off: attention-only nudges pc (0.617) but hurts lncRNA
  (0.321); attention + the ATG ORF track gives the best lncRNA of all (0.405, +0.059 over
  baseline) at a small pc cost (0.608). The ORF track and attention synergize for lncRNA --
  attention has useful transcript-wide ORF structure to attend to.
- Recommendation for the non-canonical-ORF goal (lncRNA-weighted): ATG ORF track + attention
  (best lncRNA 0.405), and drop the non-AUG channel. For pc-weighted profile, the ATG ORF track
  alone (0.619) is best. One fold only; multi-fold confirmation is the follow-up.

Count-head magnitude (P) is stable across all six configs (baseline 0.780, orf 0.791, orf_v2
0.787, attn 0.793, orf_attn 0.772, orf_v2_attn 0.761): it is driven by the RNA-seq coverage
input, not the ORF track or attention, so those changes move profile shape and lncRNA but not
magnitude. This matches the Task 6 double dissociation (coverage owns magnitude). The attention
rows' count was backfilled by `scripts/eval_extra_attn.sbatch` after the in-job eval_extra OOMed
at the old default budget 24000; eval_extra now caps its batch budget at the run's training
budget (16000 for attention runs).

## Task 9: multi-fold confirmation of the improvement matrix -- DONE

The Task 8 matrix was fold 0 only (lncRNA n=264). All 6 configs were extended to the 5
gene-disjoint chromosome folds so the lncRNA deltas rest on the pooled ~1,239 transcripts (each
scorable transcript held out once). `scripts/train_multifold.sbatch` (SLURM array 35137714, 25
new cells: baseline folds 0-4 fresh in the improve lineage + the other 5 configs folds 1-4, fold
0 reused; per-tx `pertx.tsv` backfill for the reused fold-0 runs was job 35137713). All 30 cells
COMPLETED, none failed. Settings match `backend_cmp/rinalmo_f0` exactly; `baseline_f0` reproduced
it (pc 0.611 vs 0.612), confirming the harness. Aggregated by `scripts/aggregate_multifold.py`
(`results/improve/summary_multifold.{tsv,json}`).

Pooled = median over the concatenated per-transcript profile Pearson from all 5 folds (each
scorable transcript held out once). The three ORF classes are tracked TOGETHER (deltas vs baseline
in parentheses): pc whole-transcript n=33,257; lncRNA whole-transcript n=1,239; uORF (pc 5'UTR,
translated-uORF tx) n=22,495; dORF (pc 3'UTR) n=12,286. The uORF/dORF metric and its rationale are
in Task 9c; all 8 configs including the CTG track (Task 9b) are shown.

| config | pc | lncRNA | uORF | dORF |
|---|---|---|---|---|
| baseline | 0.600 | 0.349 | 0.559 | 0.291 |
| orf (ATG) | 0.603 (+0.003) | 0.405 (+0.056) | 0.584 (+0.024) | 0.329 (+0.038) |
| orf_v2 (nonAUG) | 0.607 (+0.008) | 0.382 (+0.033) | 0.593 (+0.033) | 0.329 (+0.038) |
| orf_ctg (ATG+CTG) | 0.603 (+0.003) | 0.384 (+0.035) | 0.577 (+0.018) | 0.320 (+0.030) |
| attn | 0.607 (+0.008) | 0.378 (+0.029) | 0.581 (+0.022) | 0.307 (+0.016) |
| orf_attn (ATG+attn) | 0.612 (+0.012) | **0.409** (+0.060) | 0.604 (+0.045) | 0.340 (+0.049) |
| **orf_v2_attn (nonAUG+attn)** | **0.618** (+0.018) | 0.404 (+0.055) | **0.618** (+0.059) | 0.348 (+0.057) |
| orf_ctg_attn (ATG+CTG+attn) | 0.615 (+0.015) | 0.392 (+0.043) | 0.610 (+0.051) | **0.350** (+0.059) |

lncRNA periodicity (frame), whole-tx: baseline 0.072, orf 0.083, orf_attn 0.112 (best), orf_v2_attn
0.081 -- attention on the ATG track most improves lncRNA frame prediction.

**Best all-around config: `orf_v2_attn` (non-AUG start channel + attention).** Tracking the three
regions together resolves the earlier region-by-region picks into one. `orf_v2_attn` is best on pc
(0.618) and uORF (0.618), essentially best on dORF (0.348, tied with orf_ctg_attn 0.350), and
statistically tied for best lncRNA (0.404 vs orf_attn 0.409 -- a +0.005 gap on the noisy n=1,239
lncRNA set against a +0.014 uORF advantage on the tight n=22,495 set). Head to head it beats
`orf_attn` on pc (+0.006), uORF (+0.014), and dORF (+0.008), losing only within-noise lncRNA
(+0.006).

Findings, and how they revise the fold-0 (Task 8) conclusions:

- **The ATG ORF track's lncRNA benefit is confirmed and larger than fold 0 said: +0.056 pooled**
  (fold 0 reported +0.038). This is the single most robust effect. pc is essentially flat across
  all configs (0.600 to 0.618, within the +/-0.005 to 0.009 fold std).
- **`orf_attn` is the best lncRNA config, 0.409 (+0.060), unchanged recommendation** -- but the
  decomposition is corrected. The lift is ORF-track-dominated (+0.056); attention adds only a
  marginal +0.004 on top of the ORF track for lncRNA *shape* (fold 0 had inflated this to +0.021).
  Attention's real value is elsewhere: it lifts pc (orf 0.603 -> orf_attn 0.612) and, most
  clearly, lncRNA **periodicity/frame** (0.083 -> 0.112, the best of any config).
- **RETRACTED: "non-AUG hurts lncRNA."** Fold 0 put `orf_v2` at 0.334, *below* baseline (0.346),
  reading as harm. Pooled, `orf_v2` is 0.382, **+0.033 above baseline** -- non-AUG *helps*, just
  less than ATG-only (0.405). The near-cognate channel dilutes rather than reverses the benefit.
- **RETRACTED: "attention alone hurts lncRNA."** Fold 0 had `attn` at 0.321 (below baseline);
  pooled it is 0.378, **+0.029**. Another fold-0 artefact.
- **Attention compensates for the non-AUG dilution:** attention lifts `orf_v2` more (+0.022,
  0.382 -> 0.404) than it lifts `orf` (+0.004), so `orf_v2_attn` (0.404) nearly catches
  `orf_attn` (0.409), and has the best pc of all (0.618).

Root cause of the fold-0 errors: baseline lncRNA swings 0.306 / 0.379 / 0.307 / 0.376 / 0.422
across folds (fold 0 was the joint-lowest), so a single-fold delta on n=264 was an unreliable
reference. **Standing lesson: never quote a single-fold lncRNA delta as a finding here.**

Recommendation (per-region; the consolidated all-around pick is above and in Task 9c): for lncRNA
shape alone, `orf_attn` is marginally best (0.409, best frame 0.112); but once uORFs are scored
(Task 9c), the non-AUG channel earns its place and `orf_v2_attn` becomes the best all-around config.
The remaining check for alternative-start localization is the non-canonical-ORF truth-set eval
(not yet built).

### Task 9b: ATG+CTG ORF track (v3) -- DONE (see Task 9c uORF result)

Follow-up motivated by the multifold correction: the non-AUG track (`orf_v2`) does not hurt
lncRNA as fold 0 suggested (+0.033 pooled), but it is weaker than ATG-only (+0.056), consistent
with its 10-codon near-cognate channel diluting the signal. CTG is by far the dominant
alternative start biologically, so this tests a *selective* ATG+CTG track. Unlike `orf_v2` (which
kept ATG-only occupancy and only added a graded start channel), the v3 track (`build_orf_track.py
--mode atgctg` -> `orf_track_v3.npy`) EXTENDS the ORF occupancy to CTG-initiated ORFs, with a
2-level start channel (ATG=1.0, CTG=0.5). Occupancy density is 0.555/frame (vs ATG-only 0.377),
informative and far from the ~0.9 that opening all 10 near-cognates would give. Two configs
(`orf_ctg`, `orf_ctg_attn`) trained at all 5 folds (`scripts/train_ctg_multifold.sbatch`, array
35174489 gated `afterok` on the v3 build 35174488), directly comparable to `orf` and `orf_v2` on
the pooled basis. `aggregate_multifold.py` includes both. Question: does ATG+CTG beat ATG-only
(0.405) by capturing real CTG-initiated translation, or does CTG mostly add noise like the full
near-cognate set.

## Task 9c: uORF-aware evaluation -- DONE

uORFs are the most common alternative ORF, but the whole-transcript pc metric is CDS-dominated
and barely reflects them. Added a 5'UTR (uORF) and 3'UTR (dORF) profile-Pearson metric
(`build_tx2cds.py` -> `data/tx2cds.tsv`; `eval_extra.py` scores the window `[0, cds_start)`;
`aggregate_multifold.py` pools it). Training stays annotation-free -- CDS coords define the eval
window only. Scored set: pc tx with 5'-complete CDS, 5'UTR >= 30 nt, >= 20 5'UTR P-sites
(translated uORF); ~22,000 pooled (~4,400/fold, ~18x the lncRNA n, so tightly powered). Applied to
the already-trained models by re-eval (`reeval_uorf.sbatch` job 35176144 for the 6 core configs;
`reeval_ctg_uorf.sbatch` job 35176193, afterok the CTG array, for the 2 CTG configs). The pooled
uORF matrix (does the ORF track / CTG / attention improve uORF prediction, not just whole-tx
shape) lands here on completion. Note: `orf_ctg` (ATG+CTG) completed its 5 folds -- on
whole-transcript lncRNA it is 0.384 (+0.035), below ATG-only (0.405) and level with the diluted
non-AUG track (0.382); the uORF metric is the more relevant test for CTG (re-eval DONE in the uORF
result below -- all 8 configs).

### uORF result (all 8 configs, DONE)

Pooled 5'UTR (uORF) profile Pearson, n=22,495, baseline 0.559 (per-fold std 0.005 to 0.018 --
tightly powered):

| config | uORF pooled | Δ uORF | dORF pooled | (whole-tx lncRNA Δ) |
|---|---|---|---|---|
| baseline | 0.559 | -- | 0.291 | -- |
| orf (ATG) | 0.584 | +0.025 | 0.329 | (+0.056) |
| orf_v2 (nonAUG) | 0.593 | +0.034 | 0.329 | (+0.033) |
| orf_ctg (ATG+CTG) | 0.577 | +0.018 | 0.320 | (+0.035) |
| attn | 0.581 | +0.022 | 0.307 | (+0.029) |
| orf_attn (ATG+attn) | 0.604 | +0.045 | 0.340 | (+0.060) |
| **orf_v2_attn (nonAUG+attn)** | **0.618** | **+0.059** | 0.348 | (+0.055) |
| orf_ctg_attn (ATG+CTG+attn) | 0.610 | +0.051 | 0.350 | (+0.043) |

**The uORF metric reverses the non-AUG verdict.** On whole-transcript lncRNA the non-AUG track was
weaker than ATG-only (0.382 < 0.405). On the uORF metric it is **stronger**: orf_v2 0.593 > orf
0.584, and orf_v2_attn 0.618 > orf_attn 0.604. This is biologically expected -- uORFs are enriched
for non-AUG (especially CUG) initiation, so the non-AUG start channel helps precisely where uORFs
live (the 5'UTR), even though it diluted the whole-transcript signal. The whole-transcript metric
was simply the wrong probe for uORFs, which is exactly why this metric was needed.

- **Best uORF config: `orf_v2_attn` (non-AUG track + attention), 0.618 (+0.059)** -- also the best
  pc (0.618). For the uORF/non-canonical-start goal this beats the ATG-only `orf_attn` (0.604).
- **Attention adds more to uORF than to lncRNA:** on top of the ORF track it gives +0.020 (ATG) to
  +0.025 (non-AUG) for uORF, vs only +0.004 for whole-tx lncRNA. Attention helps uORF prediction.
- **dORF (3'UTR)** follows the same ordering at lower absolute values (baseline 0.291 to
  orf_v2_attn 0.348), consistent with rarer/weaker 3'UTR translation.

**CTG (ATG+CTG track) does not beat the full non-AUG channel on uORF, and the occupancy extension
slightly hurts.** `orf_ctg` (0.577) is *below* ATG-only (0.584) and well below the non-AUG channel
`orf_v2` (0.593); with attention `orf_ctg_attn` (0.610) sits between `orf_attn` (0.604) and
`orf_v2_attn` (0.618). The lesson is about *encoding*: the v3 track both extended the hard ORF
occupancy to CTG and shrank the start channel to 2 codons (ATG=1.0, CTG=0.5), whereas `orf_v2` keeps
ATG-only occupancy and adds a graded start-propensity channel over all 10 near-cognates. The non-AUG
benefit for uORFs comes through the soft graded start channel, not through extending the hard
occupancy to CTG (which adds density and slightly hurts). Caveat: this conflates the two changes; a
clean isolation (ATG occupancy + graded start channel restricted to {ATG,CTG}) was not run, so "CTG
alone vs all near-cognates" is not fully separated -- but extending occupancy to CTG is clearly not
the win.

Implication: the recommendation is use-case-dependent. For lncRNA-ORF shape, ATG track + attention
(`orf_attn`) is best; for uORFs (the most common alternative ORF), **non-AUG start channel +
attention (`orf_v2_attn`) is best (uORF 0.618)**, delivered as a graded near-cognate start channel
over ATG-only occupancy, not by extending the ORF occupancy to CTG.

### Task 9d: concentration (spikiness) robustness -- DONE

Profile Pearson is inflated by spiky transcripts: if one nucleotide holds most of the P-sites (a
strong initiation pause), matching that single spike gives r ~ 0.99 for free (e.g. CKS1B, 71% of
P-sites on one nt, scores 0.99). `scripts/robustness_concentration.py` measures each transcript's
max single-nt fraction on the SAME window as each metric (whole tx for pc/lncRNA, 5'UTR for uORF,
3'UTR for dORF, from the packed target counts) and re-computes the per-config pooled medians on the
distributed subset (top nt < 10%).

Composition (what fraction of scored transcripts are spiky): pc is 86.9% distributed / 1.8% spiky
(median top-nt 3.6%) -- pc scores are honest, not spike-driven. **uORF is the concentrated one:
only 10.4% distributed, 33.1% spiky, median top-nt 22.7%** -- the ~0.56 to 0.62 headline uORF
numbers are inflated by these. lncRNA 41% / 9%; dORF 24% / 21%.

Verdict: **the config rankings survive the distributed-only restriction.** On distributed uORFs
(top nt < 10%, n=2,331): baseline 0.493, orf 0.511, orf_v2 0.519, orf_ctg 0.504, attn 0.510,
orf_attn 0.525, orf_v2_attn 0.543 (best), orf_ctg_attn 0.530. `orf_v2_attn` is still best, non-AUG
still beats ATG (0.519 > 0.511; 0.543 > 0.525, delta if anything larger than on the full set), and
the CTG occupancy extension still slightly hurts (0.504 < 0.511). So the non-AUG-for-uORF result is
NOT a spike artefact. Distributed pc rankings hold (pc is 87% distributed already); distributed
lncRNA keeps ATG (orf 0.433 / orf_attn 0.431) ahead of non-AUG (orf_v2_attn 0.418). The one honest
correction is the ABSOLUTE level: distributed-only uORF performance is ~0.49 to 0.54, meaningfully
below the spike-inflated ~0.56 to 0.62 headline, so real 5'UTR profile prediction on non-trivial
transcripts is ~0.5, not ~0.6. `results/improve/robustness_concentration.{tsv,json}`.

## Task 10: RNA-seq input depth normalization (transferability) -- DONE

The RNA-seq coverage input was raw pooled per-nt read depth with only a `log1p` transform (no depth
normalization), so it carried the dataset's absolute sequencing depth: fine within Fibroblast
(train/test same depth) but not transferable across datasets of different depth (the stated
cross-dataset objective). Confirmed by inspection: `rnaseq_coverage.py` adds 1 per mate per position
(raw depth), pooling is an element-wise sum, packing stores raw int32, and `dataset.py` did
`feats[:,-1]=log1p(cov)`. ACTB per-nt coverage mean ~8,733 (raw); a 10x-shallower library shifts the
whole `log1p` channel by ~log(10)=2.3.

Fix: `train.py --cov_norm global_mean` (recorded in `args.json`; default for new runs). Coverage is
divided by the dataset global mean per-nt depth (7,733 for Fibroblast; `data/packed/
coverage_norm.json`) before `log1p` -- a CPM-like depth normalization; a new dataset computes its OWN
global mean over its packed coverage. `cov_norm=raw` is kept as the legacy path so the existing runs
(Tasks 6 to 9d, all raw) still evaluate correctly (`eval_extra` defaults to raw when `args.json` lacks
the key). Backward-compatible; no model change (the count head sums the now-normalized channel).

Depth-invariance proven by construction (CPU check): feeding a 10x-deeper library and recomputing the
mean, the `global_mean` input is IDENTICAL (max |diff| 1.2e-7), while the raw input shifts by +2.285
(~log 10). So the normalized model is invariant to sequencing depth; the raw model is not.

Within-dataset preservation check DONE (`scripts/train_covnorm.sbatch`, job 35178857): baseline +
orf_v2_attn retrained at all 5 folds with `cov_norm=global_mean`, pooled and compared to the raw runs
(`scripts/compare_covnorm.py`, median profile Pearson per region):

| config      | region | raw   | normalized | diff   |
|-------------|--------|-------|------------|--------|
| baseline    | pc     | 0.600 | 0.593      | -0.007 |
| baseline    | lncRNA | 0.349 | 0.363      | +0.014 |
| baseline    | uORF   | 0.559 | 0.544      | -0.015 |
| baseline    | dORF   | 0.291 | 0.285      | -0.006 |
| orf_v2_attn | pc     | 0.618 | 0.609      | -0.009 |
| orf_v2_attn | lncRNA | 0.404 | 0.387      | -0.017 |
| orf_v2_attn | uORF   | 0.618 | 0.595      | -0.024 |
| orf_v2_attn | dORF   | 0.348 | 0.340      | -0.008 |

Honest read: depth normalization costs a SMALL but real amount within-dataset (~0.006 to 0.024 median
Pearson, mostly ~0.01; baseline lncRNA is the one exception at +0.014). It is not the "no change" a
first guess would expect -- the raw absolute depth carried a little within-Fibroblast-useful signal (a
mild abundance prior on which transcripts have clean profiles) that normalizing away discards. This is
the expected price of depth-invariance and it is worth paying: a model that only works at one
sequencing depth is useless cross-dataset (the actual objective), and the depth x10 check proves the
normalized input is invariant by construction. The cost is ~1.5 pp of median Pearson to buy
transferability. (The 13-tx pc count difference 33257 -> 33244 is the chrM exclusion landing in the
normalized runs' eval, negligible on the median per Task 11.) Note: this normalizes the depth SCALE
only; cross-dataset/cross-tissue coverage SHAPE differences (3' bias, library prep) remain a separate
transfer concern, now being tested directly by the 8-tissue LOTO (job 35180395, cov_norm=global_mean).

## Task 11: mitochondrial gene exclusion -- DONE

Mitochondrial protein-coding genes were slipping into the universe. `Mt_rRNA` / `Mt_tRNA` had
been removed upstream by the ncRNA drop-list, but the 13 chrM mRNAs (MT-ND1/2/3/4/4L/5/6,
MT-CO1/2/3, MT-CYB, MT-ATP6/8, tx ENST00000361227.2 .. ENST00000362079.2) passed through as
`protein_coding`. They should not be there: mitochondrial mRNAs are translated by the mitoribosome
with a different genetic code, are leaderless (no 5'UTR, `utr5=0`, so no uORFs), and carry none of
the cytoplasmic 3-nt periodicity that the ORF-candidate track and the FM embeddings encode. They are
extreme abundance outliers (MT-CO1 TPM 13,722; MT-CYB 69,976 pooled P-sites) and the model predicts
them poorly (whole-tx profile Pearson median 0.124 for the 13 vs 0.618 for real pc genes) -- exactly
because the standard-code ORF track and nuclear-trained embeddings do not apply to them.

Fix: drop chrM at two points. (1) `dataset.excluded_tx()` reads `data/tx2biotype.tsv` and
`load_split` filters every train/val/test list, so no chrM tx enters training or eval -- this is the
effective gate and covers the existing Fibroblast pack and all future LOTO runs (shared split). (2)
`define_universe_and_fasta.py` skips `chrom in {"chrM"}` so future universe / FASTA / embedding
builds never include them. Confirmed: 13 usable MT pc genes dropped, zero chrM leak across all 5
folds; split sizes essentially unchanged (train ~22.0k, val/test ~7k).

Impact on the already-reported matrix (Tasks 6 to 10, all computed with the 13 present): negligible.
pc profile-Pearson median 0.6180 with MT vs 0.6181 without (13 of 33,257 pooled); uORF/dORF metrics
unaffected (MT mRNAs are leaderless, already excluded by the `utr5 >= 30` window); config rankings
unchanged. So those numbers stand as reported; the exclusion is active for every future run.

## Biotype / chromosome table

`data/tx2biotype.tsv` built from GENCODE v49 GTF: 507,365 transcripts
(`tx_id, gene_id, gene_name, chrom, strand, transcript_type, gene_type, length`; length
= summed exon lengths). About 2,285 of the 509,650 hd5-axis transcripts have no GTF
transcript-feature match and carry `NA`/`unknown` in the summary (843 of them have
nonzero signal). Reused for the summary join and the later held-out-chromosome split.

## Task 12: leave-one-tissue-out (cross-tissue transfer) -- DONE

The project objective: does the model, trained on several tissues, predict the per-nt Ribo-seq
profile of a tissue it never saw? Pilot: train on 8 Chothani tissues, hold out **Hepatocytes**
(design in `LOTO_PLAN.md`). The FM embeddings + ORF track are sequence-based / tissue-independent,
so the fixed 36,668-tx universe is reused with zero re-extraction; per-tissue expression enters
through the RNA-seq coverage input. Each tissue has its own packed target + coverage + depth
normalizer (`data/packed_<Tissue>/`, `cov_norm=global_mean` per tissue); training concatenates the
8 held-in tissues (`scripts/train_loto.py`, `dataset.loto_split`, tissue held out not chromosome),
capped at 6,000 tx/tissue for balance (deep tissues would otherwise dominate). chrM excluded
(Task 11). Held-out Hepatocytes is the *deepest* target (2.46e9 pooled P-sites) -- a clean test but
its depth flatters the numbers vs a shallower held-out tissue would.

Consolidated held-out Hepatocytes profile Pearson (`scripts/loto_table.py`, from pertx.tsv), against
the **matched** within-Fibroblast ceiling at the same normalization (Task 10 covnorm, NOT the raw
Task 9 numbers):

| backend | config | pc | lncRNA | uORF | dORF |
|---------|--------|------|--------|------|------|
| RiNALMo | baseline    | 0.565 | 0.361 | 0.459 | 0.209 |
| RiNALMo | orf_v2_attn | **0.611** | **0.438** | 0.534 | 0.234 |
| within-Fib (norm) | baseline    | 0.593 | 0.363 | 0.544 | 0.285 |
| within-Fib (norm) | orf_v2_attn | 0.609 | 0.387 | 0.595 | 0.340 |
| Orthrus | baseline    | 0.585 | 0.376 | 0.494 | 0.214 |
| Orthrus | orf_v2_attn | **0.621** | **0.428** | 0.530 | 0.239 |

Headline (RiNALMo, orf_v2_attn): **cross-tissue transfer is essentially lossless for pc (0.611 vs
0.609 within-tissue) and exceeds within-tissue on lncRNA (0.438 vs 0.387)** -- a model predicts an
unseen tissue's Ribo-seq profile about as well as it predicts its own. uORF (-0.06) and dORF (-0.11)
keep real gaps: 5'/3'UTR translation is more tissue-specific. The abundance/count head transfers very
well (count Pearson ~0.83). Architecture matters cross-tissue: **baseline transfers *worse* on pc
(0.565 vs 0.593 within-tissue, -0.03)** -- the ORF track + attention is what buys the near-lossless
transfer (+0.045 pc, +0.075 uORF over baseline), not a marginal gain.

Concentration robustness (`scripts/robustness_concentration_loto.py`), the credibility check that the
pc number is not spike-inflated: 75% of pc transcripts are *distributed* (<10% of P-sites on any one
nt), and orf_v2_attn's distributed-pc median is **0.610 = the all-pc 0.611**, so the headline is
genuine (0.610 distributed vs 0.609 within-tissue ceiling). lncRNA holds (distributed 0.441 >= all
0.438). uORF is partly spike-inflated (44.5% of uORF windows spiky): honest distributed uORF is 0.512,
below all-0.534, so the transfer gap to within-tissue (0.595) is real. orf_v2_attn beats baseline on
distributed transcripts everywhere (pc +0.040, lncRNA +0.078, uORF +0.052) -- the architecture edge is
not a spike artifact.

Infrastructure note: the first pilot (job 35180395) died mid-training with no traceback (exit 120) --
root cause was the 15 TB carpenterlab ceph group quota filling during the run, so checkpoint writes
failed (see memory `feedback_group_ceph_15tb_quota`). Fixed by reclaiming redundant intermediates and
hardening `train_loto.py` with per-epoch `last.pt` + `--resume`; the re-run (35183966) completed clean
on the same config. Orthrus arm (job 35224524, `--emb_backend orthrus`) COMPLETED and is backfilled
above: **the efficiency edge holds cross-tissue** -- Orthrus orf_v2_attn reaches pc 0.621 (edging
RiNALMo's 0.611) and lncRNA 0.428 at 40% the embedding dimension, with the same architecture pattern
(orf_v2_attn > baseline on pc: 0.621 vs 0.585). Task 7's within-tissue Orthrus advantage transfers.

## Task 13: localization eval -- does the model concentrate signal on translated ORFs? -- DONE

Profile Pearson scores the *shape* of the predicted profile but never asks the sharpest question this
project cares about: on a held-out transcript, does the model localize predicted Ribo-seq signal to
the ORFs that are genuinely translated -- especially non-canonical ones -- rather than to the many
untranslated candidate AUG ORFs sharing the same transcript? `scripts/eval_localization.py` scores
this against RiboCode ORF calls (`<Tissue>_collapsed.txt`; design in methods.md sec 5G). Every maximal
AUG..in-frame-stop ORF on a test transcript is a candidate (length >= 30 nt); a candidate is POSITIVE
if its start matches a RiboCode call. The orf_track marks candidate positions + frames but NOT which
translate, so the discrimination must come from the model. Two length-fair scores from the predicted
softmax profile: `density` (mean per-nt) and `frame0` (in-ORF frame-0 fraction, the periodicity
signature, ~1/3 with no translation). RiboCode here was run AUG-only, so this is translation
localization of AUG ORFs, **not** a non-AUG-start discovery test.

Within-tissue fold-0 (held-out chromosomes; **2,687 translated** [1,652 annotated + 1,035
non-canonical]). Internal in-frame alt-start candidates -- AUG ORFs re-initiating strictly inside an
annotated CDS in its own reading frame -- are collapsed by default (methods.md sec 5G): they are
5'-truncations of the same protein, not distinct ORFs. This drops **61,911** candidates from the
negative set (0 non-canonical positives; 144 annotated calls that land at an internal in-frame position
via a RiboCode-vs-`tx2cds` coordinate offset are kept), leaving **175,422 untranslated** candidate AUG
ORFs. dORFs (in-frame past the stop and off-frame overlapping), out-of-frame CDS overlaps, in-frame
N-terminal extensions, and all UTR / lncRNA ORFs are kept. AUROC separating translated from untranslated
candidates:

| config | stratum | pred frame0 | pred frame0 (len-ctrl) | pred density | ORF-length floor | obs frame0 (raw) |
|--------|---------|-------------|------------------------|--------------|------------------|------------------|
| baseline    | all           | 0.899 | 0.856 | 0.767 | 0.839 | 0.886 |
| orf_v2_attn | all           | 0.914 | 0.870 | 0.770 | 0.839 | 0.886 |
| baseline    | non-canonical | 0.778 | 0.770 | 0.675 | 0.616 | 0.864 |
| orf_v2_attn | non-canonical | 0.813 | **0.807** | 0.697 | 0.616 | 0.864 |

Headline: **the non-canonical stratum is the honest test, and the model reaches 94% of the observed
ceiling on it.** On non-canonical ORFs (short: median 102 nt) the **length-controlled predicted-frame0
AUROC is 0.807 (orf_v2_attn) / 0.770 (baseline)**, against a length-controlled observed-periodicity
ceiling of 0.854 -- so orf_v2_attn recovers **94.4% of the data ceiling** (baseline 90.1%), localizing
signal to translated non-canonical ORFs using learned periodicity, not length (the raw ORF-length floor
is only 0.616 and the length-controlled score sits well above it). orf_v2_attn beats baseline on every
non-canonical score (len-ctrl frame0 +0.037, density +0.022), the same architecture edge the
whole-transcript LOTO showed. Even the ALL stratum now shows real discrimination beyond length
(length-controlled AUROC 0.870 vs a 0.839 length floor): collapsing the internal in-frame CDS fragments
(below) removes the long, legitimately-periodic negatives that would otherwise make the ALL stratum a
pure length artifact.

Fidelity (predicted vs observed frame0 on translated ORFs) -- **the model reproduces 3-nt periodicity
it was never explicitly trained on** (the loss is multinomial NLL on raw per-nt counts, no periodicity
term):

| config | stratum | n | median pred frame0 | median obs frame0 | frame0 Pearson |
|--------|---------|---|--------------------|-------------------|----------------|
| orf_v2_attn | annotated     | 1652 | 0.812 | 0.801 | 0.416 |
| orf_v2_attn | non-canonical | 1035 | 0.553 | 0.706 | 0.495 |
| baseline    | annotated     | 1652 | 0.796 | 0.801 | 0.380 |
| baseline    | non-canonical | 1035 | 0.504 | 0.706 | 0.432 |

For canonical CDS the model matches the observed periodicity almost exactly (predicted 0.81 vs observed
0.80). For non-canonical ORFs it predicts periodicity well above the 1/3 null (0.55) but *undershoots*
observed (0.71) -- it knows they translate, less confidently. orf_v2_attn is closer than baseline
everywhere (all-positives frame0 Pearson 0.625 vs 0.589). By non-canonical type (orf_v2_attn median
predicted frame0): uORF 0.630, novel 0.596, dORF 0.458, Overlap_uORF 0.458, internal 0.177. uORF and
novel are strong; dORF moderate (matching the LOTO dORF gap); `internal` is low but *correct* -- internal
out-of-frame ORFs sit inside the CDS, so their own-frame periodicity is masked by the dominant CDS
frame (observed is likewise low). Figure: `scripts/make_localization_figure.py` (3 panels: the frame0
separation, the length-controlled non-canonical AUROC, the per-type fidelity). The fair upper bound is
the *length-controlled observed ceiling* -- how well the observed Ribo-seq periodicity itself separates
translated from untranslated within length strata: **0.854** on non-canonical ORFs, so orf_v2_attn's
0.807 reaches **94% of the data ceiling** (baseline 0.770 = 90%). Caveat: AUG-only truth set -- non-AUG
start discovery needs re-running the caller with near-cognate starts (future work).

**Why the untranslated-candidate periodicity is bimodal, and why the internal in-frame alt starts are
collapsed** (the frame0 violin in figure panel A has a large cluster above 0.7, not just the 1/3 null).
Classifying each of the 237,333 enumerated untranslated candidates by its relationship to the annotated
CDS (`cds_rel` column in `localization_orfs.tsv`) shows the set mixes two populations and the model is
correct on both. The high cluster (predicted frame0 > 0.7 = 30% of candidates) is **94% in-frame with
the annotated CDS** -- 87% internal in-frame alt starts (`cds_inframe`; predicted frame0 0.83, observed
0.81, 98% above 0.7) plus 7% canonical starts RiboCode did not confirm on that transcript. These are AUG
ORFs sharing a translated CDS's reading frame: RiboCode collapses each locus to one representative call,
so the internal fragments read as "untranslated" while actually sitting inside real periodic
translation, and the model correctly predicts them periodic. Because they are 5'-truncations of the same
protein rather than distinct ORFs, they are collapsed out of the discrimination metric (methods.md sec
5G) -- removing exactly the long, legitimately-periodic negatives that otherwise contaminate the negative
set. The low cluster at the null is 56% off-frame CDS overlap (CDS P-sites are off-frame relative to
these, so predicted and observed frame0 ~0.09) + 27% 3'UTR (near-zero signal) + 12% in-frame dORFs past
the stop (`dorf_inframe`; predicted 0.41 / observed 0.38 -- real but sparse, and KEPT as distinct ORFs).
5'UTR and lncRNA candidates sit intermediate (observed 0.62 / 0.60), consistent with genuine but sparser
non-canonical translation; in-frame N-terminal extensions (`ext_inframe`, 158 candidates) are highly
periodic (predicted 0.81 / observed 0.80) and are also KEPT, since a non-AUG extension of this kind is a
genuine alt-ORF (MYC / VEGFA; methods.md sec 5G). Figure `results/figures/untranslated_decomp.png`
(`scripts/make_untranslated_decomp_figure.py`).

Cross-tissue (LOTO, held-out Hepatocytes; **10,482 translated** [7,808 annotated + 2,674 non-canonical];
1,154,734 candidates enumerated, 314,037 internal in-frame alt starts collapsed, leaving **830,215
untranslated** -- deeper and better-powered than Fibroblast). **Translated-ORF localization transfers
losslessly -- in fact slightly better on the deeper held-out tissue:**

| config | stratum | pred frame0 (len-ctrl) cross-tissue | within-Fib (len-ctrl) | pred frame0 (raw) | ORF-length floor |
|--------|---------|-------------------------------------|-----------------------|-------------------|------------------|
| baseline (rinalmo)    | non-canonical | 0.774 | 0.770 | 0.781 | 0.604 |
| orf_v2_attn (rinalmo) | non-canonical | **0.827** | 0.807 | 0.830 | 0.604 |
| baseline (orthrus)    | non-canonical | 0.809 | n/a | 0.813 | 0.604 |
| orf_v2_attn (orthrus) | non-canonical | 0.822 | n/a | 0.824 | 0.604 |

The non-canonical length-controlled predicted-frame0 AUROC on the *unseen* tissue is 0.774 (baseline) /
0.827 (orf_v2_attn), matching or exceeding within-Fibroblast (0.770 / 0.807); the raw ORF-length floor
is 0.604 and the length-controlled score sits well above it. **orf_v2_attn's edge over baseline GROWS
cross-tissue** (+0.053 vs +0.037 within-tissue) -- the same signal as the whole-transcript LOTO, that the
ORF track + attention is what buys the transfer. Fidelity: the model reproduces the (even cleaner)
Hepatocytes periodicity on a tissue it never trained on -- annotated predicted frame0 0.824 vs observed
0.853; non-canonical predicted 0.618 vs observed 0.751 with frame0 Pearson 0.563 (orf_v2_attn, up from
within-tissue 0.495). Even on the length-confounded ALL stratum, orf_v2_attn's length-controlled AUROC is
0.896 -- above the 0.884 length floor -- while baseline (0.872) sits below it: the ORF track + attention
discriminates translated ORFs beyond length where the FM embedding + coverage alone cannot. Against the
length-controlled observed ceiling on the held-out tissue (**0.879** non-canonical), orf_v2_attn's 0.827
reaches **94% of the data ceiling** (baseline 0.774 = 88%) -- the model localizes non-canonical
translation on an unseen tissue nearly as well as the observed Ribo-seq data itself does. The absolute
non-canonical AUROC shifts with how the negative set is defined (uncollapsed vs internal-in-frame
collapsed), but the model-reaches-94%-of-the-observed-ceiling ratio is stable across both definitions and
both tissues, which is why the ceiling ratio is the headline.

**Orthrus backend (`--emb_backend orthrus`, LOTO Hepatocytes arm, job 35247400, completed 2026-07-12)**
completes the 4-way. On the non-canonical stratum Orthrus reaches length-controlled AUROC 0.809
(baseline) / 0.822 (orf_v2_attn) -- **92.0% / 93.5% of the 0.879 observed ceiling**. The backend pattern
differs from RiNALMo: Orthrus's baseline embedding alone already reaches 92% (vs RiNALMo baseline's 88%),
so the ORF-track + attention adds only **+0.013** on Orthrus versus **+0.053** on RiNALMo -- the richer
mature-mRNA (Mamba) embedding encodes more of the translation-localization signal up front, leaving less
for the architecture head to recover. Both attention configs converge near the ceiling (RiNALMo 0.827 =
94.1%, Orthrus 0.822 = 93.5%), with **RiNALMo + orf_v2_attn the single best**. Orthrus within-Fibroblast
fold-0 localization was not run (only the cross-tissue LOTO arm), so those cells are `n/a`. The
whole-transcript profile-Pearson LOTO 4-way is in `results/loto/loto_summary.tsv` (Orthrus orf_v2_attn pc
0.621 vs RiNALMo 0.611; lncRNA 0.428 vs 0.438 -- effectively tied, same story as localization: Orthrus and
RiNALMo are within noise on this task, both backends benefit from the ORF track + attention).

## Task 14: RiboCode drop-in test -- do the predicted profiles reproduce the ORF caller's own calls? -- DONE

The localization eval (Task 13) scores periodicity with a custom AUROC. A stronger, more concrete test:
feed the model's predicted per-nt profile into the ACTUAL ORF caller (RiboCode) in place of the
experimental Ribo-seq, and ask whether the ORF calls track the calls RiboCode makes on the real data.
Method in methods.md sec 5H. RiboCode's per-transcript P-site density is a `dict {tx_id: np.ndarray}`
consumed by `detectORF.main`; `scripts/ribocode_dropin.py` substitutes a predicted density and calls it
directly, reusing the GENCODE v49 annotation DB and RiboCode defaults (start ATG, min_AA 5, pval 0.05,
fdr_bh). Held-out Hepatocytes, orf_v2_attn RiNALMo LOTO model, model test transcripts. Three density
variants, identical caller + params, only the density differs: `real` (real pooled counts, reference),
`pred_obsdepth` (predicted shape x real per-tx depth, rounded -- isolates shape), `pred_preddepth`
(predicted shape x count-head predicted depth, rounded -- fully standalone, no experimental data).
Predicted densities are scaled to realistic depth and rounded to int (RiboCode hard-gates per-tx and
per-ORF frame0 P-site sums >= 5, and a smooth density with no zeros would look artificially periodic).
Calls compared by `(transcript_id, ORF_tstart)`: precision / recall / F1 vs the real-density calls
(overall + recall by ORF type), the real-vs-official overlap (harness validation), and confidence /
magnitude / type agreement on matched calls.

Held-out Hepatocytes, orf_v2_attn RiNALMo LOTO model, 33,918 test transcripts (0 skipped; every tx maps
1:1 to the annotation DB). **ORFs are floored at 90 nt (30 codons)** and calls are matched by **genomic
ORF locus `(gene_id, ORF_gstop)`, not by isoform.** RiboCode's per-locus collapse picks the longest ORF
per genomic stop (`detectORF.py:440-442`; density only breaks exact-start ties), so restricting the
transcript universe (34k test tx here vs the official 509k) merely reshuffles which isoform carries a
given call -- matching on the genomic locus is invariant to that. Calls at `pval_combined <= 0.05`,
ORF_length >= 90 nt, restricted to the 11,091 test-tx genes; predicted calls are additionally
enrichment-filtered at the default **0.5x uniform** (the chosen FP-minimizing operating point that
strips the spurious 3'UTR dORFs -- sec below; `compare_dropin_calls.py --min_enrichment`):

| comparison | nPred | nReal | match | precision | recall | F1 | -log10p Pearson | type concord |
|------------|-------|-------|-------|-----------|--------|-----|-----------------|--------------|
| real vs official *(harness check)*      | 13,179 | 16,555 | 13,099 | **0.994** | 0.791 | 0.881 | **1.00** | **1.00** |
| **pred (shape, real depth) vs real**    | 13,737 | 13,179 | 12,260 | 0.892 | **0.930** | 0.911 | 0.79 | 1.00 |
| **pred (standalone, no real data) vs real** | 13,534 | 13,179 | 12,226 | **0.903** | 0.928 | **0.915** | 0.82 | 1.00 |
| pred (shape) vs official                | 13,737 | 16,555 | 12,267 | 0.893 | 0.741 | 0.810 | 0.79 | 1.00 |

Pre-filter (no enrichment cut) the raw reproduction is pred-shape 0.815 / 0.946 / 0.876 and standalone
0.866 / 0.940 / 0.901 -- untuned, the model already recovers 94% of the real calls; the 0.5 floor trades
~1 point of recall for +4 precision by removing the softmax's 3'UTR leak. The 0.5 threshold was chosen on
this held-out Hepatocytes set and should be re-validated per tissue before the exact value is trusted.

**The `real vs official` precision of 0.994 proves the earlier isoform-level gap was purely
representative reassignment, not a harness error:** running detectORF on the real packed counts
reproduces 99.4% of the official calls at the genomic locus. Isoform-keyed `(transcript_id, ORF_tstart)`
the same rows read pred-shape 0.704 / 0.817 / 0.756, standalone 0.729 / 0.791 / 0.759, and real-vs-official
only 0.612 precision -- that view understates the model by penalizing collapse-representative
disagreements (the same genomic ORF carried by a different isoform in the two runs), so it is reported
here only for contrast. Separately, the 90 nt floor is load-bearing: without it the predicted variants
are much noisier (isoform-keyed pred-shape precision 0.554 / F1 0.651), because the sub-90 nt regime is
where the model over-calls (false-positive analysis below); the floor lifts precision ~+0.15 at no recall
cost.

**Headline: the model's predicted profile dropped into RiboCode reproduces the real ORF calls at F1 0.90
untuned (standalone: 94% recall at precision 0.87), and at the chosen 0.5 enrichment operating point it
is precision 0.90 / recall 0.93 / F1 0.915 -- the filter strips the spurious 3'UTR dORFs (below) for ~1
point of recall.** On every matched call RiboCode assigns the same ORF type (concordance 1.00) and a
strongly correlated confidence (pval_combined Pearson 0.79 / 0.82). The fully-standalone variant
(count-head predicted depth, no experimental data at all) is the stronger of the two throughout, because
its leaner depth over-calls less. The harness itself is validated at 0.994 real-vs-official precision, so
any residual gap is the model, not the injection.

Recall by ORF type (pred shape vs real, genomic locus, ORFs >= 90 nt, enrichment >= 0.5):

| ORF type | recovered / real | recall | (pre-filter) |
|----------|------------------|--------|--------------|
| annotated    | 10,886 / 10,947 | **0.994** | 0.995 |
| novel        |    330 / 377    | **0.875** | 0.891 |
| uORF         |    715 / 968    | **0.739** | 0.877 |
| Overlap_uORF |    281 / 474    | 0.593 | 0.605 |
| Overlap_dORF |     19 / 51     | 0.373 | 0.471 |
| dORF         |     18 / 149    | 0.121 | 0.503 |
| internal     |     11 / 213    | 0.052 | 0.056 |

Canonical CDS stays almost perfect (0.994) and novel strong (0.88); the enrichment floor is what costs
uORF recall (0.88 -> 0.74) and deliberately collapses genuine short dORF recall (0.50 -> 0.12), since
those dORFs are magnitude-indistinguishable from the spurious 3'UTR leak (below). internal stays a blind
spot (0.05, its out-of-frame ORFs sit in the CDS-masked frame the model cannot see), unaffected by the
filter.

**False positives bucket by length, and the residual class signal is spurious 3'UTR dORFs.** Before the
90 nt floor the predicted density (real depth) made 9,837 false positives (calls not in the real set):
median FP length 78 nt vs 1,003 nt for true positives, 64% of FPs <= 120 nt, FP rate falling
monotonically from 0.81 on the shortest ORFs to 0.14 on the longest, and within every ORF class the FPs
are the shorter members -- so length, not class identity, is the driver, and the floor removes ~64% of
them. **The dORF over-call is the clearest case.** RiboCode calls only 275 dORFs on the real data (0.9%
of 24,160 -- genuinely rare); the predicted density (real depth) calls 1,708, an ~8-10x over-call (the
standalone variant is milder at 673, its leaner predicted depth putting less mass over the gate).
Diagnosing the 1,638 dORF FPs: **37% sit in 3'UTR windows with ZERO real P-sites and 46% below RiboCode's
>= 5 gate** -- the softmax never zeros the 3'UTR, so scaling those tiny 3'UTR probabilities by a
well-expressed transcript's large depth and rounding yields callable counts where there is no real
signal; predicted frame0 there is 0.71 vs a real 0.50 (the model over-predicts 3'UTR phasing, 44% in the
CDS reading frame = echo enriched over the 33% null, scattered a median 142 nt past the stop, not just at
it). The 90 nt floor halves the dORF FPs (1,638 -> 878) and predicting depth instead of using the real
depth cuts the dORF calls 61% (1,708 -> 673); both together leave 421 standalone dORFs. **This over-call
is genuine, not a keying artifact:** by genomic locus the model still calls 928 dORF loci vs 149 real
(~6x), so genomic keying (which fixes the representative-isoform artifact and lifts every other class)
does NOT rescue dORFs -- they are real 3'UTR hallucinations. This is the sole remaining driver of the
0.82 precision (annotated/uORF/novel recall is 0.88-0.995). **Caveats.** (1) The dORF over-call caps
precision at 0.82; the count-head-depth (standalone) variant is the practical mitigation (precision
0.87), and a predicted-density or UTR-mass threshold would target it directly. (2) AUG-only truth set,
so this is translation-call concordance, not non-AUG discovery. Runners: `dump_pred_profiles.sbatch`
(job 35245771), `ribocode_dropin.sbatch` (35245772), then `compare_dropin_calls.py --key genomic`;
metrics in `results/loto/orf_v2_attn_holdout_Hepatocytes/dropin/dropin_metrics.json`.

**Minimizing the dORF false positives (2026-07-12).** The spurious dORFs are the diffuse
low-probability softmax leak into the 3'UTR: as predicted enrichment over uniform (mean predicted
probability x transcript length), CDS is 1.76x and uORF 1.62x, while genuine dORFs are 0.15x and
spurious dORF FPs 0.06x -- so a predicted-magnitude threshold removes them. Two thresholds were tested
(genomic locus, ORFs >= 90 nt; `ribocode_dropin.py --floor_mult` + `analyze_floor_sweep.py`, job
35246903):

| config | precision | recall | F1 | dORF FPs | uORF recall |
|--------|-----------|--------|-----|----------|-------------|
| baseline (obsdepth)                    | 0.815 | 0.946 | 0.876 | 853 | 0.877 |
| + per-position floor 1.0               | 0.873 | 0.936 | 0.903 | 204 | 0.800 |
| standalone (preddepth)                 | 0.866 | 0.940 | 0.901 | 359 | 0.798 |
| standalone + floor 1.0                 | 0.888 | 0.933 | 0.910 | 155 | 0.751 |
| **standalone + mean-density >= 0.5x (default)** | **0.903** | 0.928 | **0.915** | **60** | 0.693 |
| standalone + mean-density >= 0.7x      | 0.912 | 0.923 | 0.917 |  36 | 0.643 |
| standalone + mean-density >= uniform   | 0.921 | 0.906 | 0.914 |  27 | 0.554 |

- **Per-position density floor** (zero predicted positions below `floor_mult` x uniform, re-run
  RiboCode): a gentle denoise -- recall barely moves (0.94), F1 up to 0.91, dORF FPs cut ~76%
  (853 -> 204 / 359 -> 155), and it PRESERVES uORF recall (0.80) because a genuine uORF's frame0 peaks
  survive the floor even when its mean is below uniform.
- **Call-level mean-density filter** (keep calls whose mean predicted density >= tau x uniform): the
  aggressive lever -- cuts dORF FPs ~97% (853 -> 27) and lifts precision to 0.92, at a larger
  uORF-recall cost (0.55-0.64). It SUBSUMES the floor (floor + filter == filter alone).
- **Standalone (count-head predicted depth) dominates real depth** at every operating point (leaner
  depth -> fewer 3'UTR counts over RiboCode's >= 5 gate). annotated recall stays 0.985-0.999 throughout.

**Chosen default: standalone + mean-density >= 0.5x uniform** (`compare_dropin_calls.py --min_enrichment
0.5`, now the default) -- precision 0.90 / F1 0.915, spurious dORFs cut 853 -> 60 (-93%) while uORF recall
holds at 0.69; dial to 0.7-1.0 for more precision (0.91-0.92) at more uORF cost, or the per-position floor
for a gentler denoise. **The CDS-relative alternative was tested and does NOT help:** normalizing each
ORF's density to the transcript's own median or 75th-percentile density (instead of to uniform) is a
WORSE discriminator of genuine uORF vs spurious dORF (AUROC 0.78 / 0.84) than the plain uniform-enrichment
score (0.947) -- the uORF cost is a threshold-level choice, not a metric problem, so a lower uniform
threshold is the right lever. Figure `results/figures/dropin_floor_sweep.png`.

**Transfer check (2026-07-12): the 0.5 threshold holds on a second tissue.** Re-running the identical
pipeline on the within-Fibroblast fold-0 model (`results/improve/orf_v2_attn_f0`, a different, deep
tissue) + Fibroblast RiboCode calls: the harness validates (real-vs-official precision 0.991, matching
Hepatocytes' 0.994), and the enrichment sweep (standalone, genomic, >= 90 nt) tracks Hepatocytes -- >= 0.5
cuts dORF FPs ~87% (409 -> 52, cf. 359 -> 60 on Hepatocytes), F1 rises then plateaus from 0.3-1.0 (knee at
~0.3-0.5 in both), precision gains ~+0.10 from 0 -> 0.5, and uORF recall costs ~0.15-0.19, all matching.
Absolute level is ~6 points lower on Fibroblast (standalone enr>=0.5 precision 0.842 / recall 0.854 / F1
0.848 vs Hepatocytes 0.903 / 0.928 / 0.915) -- a genuine tissue/regime effect (Fibroblast over-calls more
at baseline, and this is a held-out-chromosome fold, not a held-out-tissue LOTO) -- but the tradeoff shape
and knee location are the same, so 0.5 is a robust operating point, not a Hepatocytes artifact. A pure
held-out-tissue LOTO replicate (rather than within-tissue) would be the last confirmation.

**Backend comparison (2026-07-13): the Orthrus drop-in matches RiNALMo, tied on novel ORF calling.** The
drop-in was re-run end-to-end on the Orthrus LOTO `orf_v2_attn` model (held-out Hepatocytes, identical
harness: dump predicted profiles job 35250402 -> `ribocode_dropin.py` 3 variants job 35250403 -> compare
job 35250404; genomic locus, >= 90 nt, 0.5x enrichment). The real-vs-official harness check is byte-identical
to RiNALMo's (precision 0.994 / recall 0.791 / F1 0.881 -- both consume the same real calls), so any residual
difference is the model, not the injection. Overall reproduction sits within noise of RiNALMo:

| comparison | backend | precision | recall | F1 |
|------------|---------|-----------|--------|-----|
| pred (shape, real depth) vs real         | orthrus | 0.892 | 0.935 | 0.913 |
| pred (shape, real depth) vs real         | rinalmo | 0.892 | 0.930 | 0.911 |
| pred (standalone, no real data) vs real  | orthrus | 0.903 | 0.930 | 0.916 |
| pred (standalone, no real data) vs real  | rinalmo | 0.903 | 0.928 | 0.915 |

Recall by ORF type (pred shape, real depth, genomic locus, >= 90 nt, enrichment >= 0.5) shows the backends
trade classes to a net wash:

| ORF type | Orthrus | RiNALMo | delta (orth - rina) |
|----------|---------|---------|---------------------|
| annotated    | 0.995           | 0.994           | +0.000 |
| uORF         | **0.807**       | 0.739           | **+0.068** |
| Overlap_uORF | 0.618           | 0.593           | +0.025 |
| novel        | 0.862 (325/377) | **0.875** (330/377) | **-0.013** |
| Overlap_dORF | 0.353           | 0.373           | -0.020 |
| dORF         | 0.094           | 0.121           | -0.027 |
| internal     | 0.023           | 0.052           | -0.028 |

**On novel ORF calling the backends are a dead heat: Orthrus recovers 325/377 vs RiNALMo 330/377, a 5-ORF
difference out of 377.** Orthrus's one clear edge is uORFs (+0.068) -- the mature-mRNA Mamba embedding
front-loads 5'UTR translation signal -- while RiNALMo is marginally ahead on the sparse, hard classes
(novel, dORF, internal, all rare and near the noise floor). Aggregate F1 is tied at ~0.91, the same
within-noise verdict the localization eval (Task 13) reached: on the concrete ORF-calling test the FM
backbone choice does not move novel-ORF recovery. Metrics in
`results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes/dropin/dropin_metrics.json`.

## Task 15: one-hot control -- how much lift do the FM embeddings actually provide? -- DONE

The whole pipeline conditions on a frozen RNA foundation-model embedding as the sequence input. The
sharpest test of whether that buys anything is to replace it with a raw `(L,4)` A/C/G/T one-hot and hold
everything else fixed (same dilated-CNN + transformer, same ORF track, same RNA coverage, same
gene-disjoint splits) -- so the delta is exactly the value of the FM pretraining over raw sequence.
Implemented as an `onehot` backend in `dataset.py` (computed on the fly; the universe is the rinalmo n
orthrus intersection so the LOTO splits are byte-identical to the FM runs). Trained the same
held-out-Hepatocytes LOTO in both configs (`baseline` = no ORF track; `orf_v2_attn` = + ORF track v2 + the
2-layer transformer), then ran the identical three evals: profile Pearson (Task 10), non-canonical
localization AUROC (Task 13), and the RiboCode drop-in concordance (Task 14).

**Headline: the FM embeddings provide NO lift. One-hot matches or beats both RiNALMo and Orthrus on every
metric, and is often slightly better.** The `orf_v2_attn` comparison (the apples-to-apples config):

| metric | onehot | rinalmo | orthrus | best |
|--------|--------|---------|---------|------|
| profile Pearson pc                 | **0.640** | 0.611 | 0.621 | onehot |
| profile Pearson lncRNA             | **0.447** | 0.438 | 0.428 | onehot |
| profile Pearson uORF               | **0.559** | 0.534 | 0.530 | onehot |
| localization non-canon AUROC (len-ctrl) | **0.832** (94.6% of ceiling) | 0.827 | 0.822 | onehot |
| localization novel-type AUROC      | **0.874** | 0.871 | 0.865 | onehot |
| drop-in standalone F1              | **0.923** | 0.915 | 0.916 | onehot |
| drop-in precision / recall        | 0.925 / 0.920 | 0.903 / 0.928 | 0.903 / 0.930 | onehot (precision) |
| drop-in novel recall              | 0.862 (325/377) | **0.875 (330/377)** | 0.862 (325/377) | rinalmo (+5 ORFs) |

One-hot wins profile Pearson on **every** stratum in **both** configs (baseline pc 0.587 vs 0.585/0.565;
lncRNA 0.405 vs 0.376/0.361), wins the non-canonical localization AUROC in both configs (baseline 0.815 vs
0.809/0.774), and wins the drop-in F1 (0.923, driven by the best precision, 0.925). The **only** cell where
an FM is nominally ahead is drop-in novel recall, where RiNALMo recovers 330/377 vs one-hot's 325/377 -- a
5-ORF difference, within noise, and one-hot ties Orthrus exactly there. So across profile shape, translated-
ORF discrimination, and ORF-call concordance, there is no metric on which the frozen FM embeddings
meaningfully beat a raw one-hot.

**Interpretation.** The embeddings are *frozen*, i.e. a fixed nonlinear compression of the sequence
optimized for the FM's own pretraining objective (masked-LM / contrastive), whereas one-hot preserves exact
per-nt codon and frame identity losslessly. With this much supervised signal (8 tissues x per-nt counts),
the dilated CNN + transformer just learns the task-relevant features directly from raw sequence, and the FM
compression only discards information. The lncRNA stratum is the tell: it is exactly where FM pretraining is
supposed to help most (sparse, "needs language"), yet one-hot's margin there is among the largest
(orf_v2_attn 0.447 vs 0.428/0.438; baseline 0.405 vs 0.376/0.361). Notably one-hot's advantage *grows* with
the orf_v2_attn architecture (pc margin over Orthrus +0.002 baseline -> +0.019 orf_v2_attn): the FM already
encodes some of what the ORF track + attention supply, so those additions help the FM less, while one-hot
starts with nothing and ends ahead. **Caveats.** (1) This is the frozen-embedding regime; a *fine-tuned* FM
backbone is a different experiment (impractical for the 650M/512-d models on this budget, and the point here
is the deployed pipeline uses frozen embeddings). (2) The one nominal FM win (novel recall +5 ORFs) is
noise-level. (3) AUG-only truth set throughout, as elsewhere.

**Practical consequence:** for per-nt Ribo-seq P-site prediction, the model can drop the FM-embedding
dependency entirely -- one-hot is as good or better, ~5M params either way, no per-token embedding
extraction step (the slowest and most storage-heavy part of the pipeline), and it runs on CPU. Runs:
`train_loto_onehot.sbatch` (jobs 35251300_0/_1), evals `eval_localization_onehot_cpu.sbatch` (35310632),
`dump_pred_profiles_cpu.sbatch` -> `ribocode_dropin.sbatch` -> compare (35310629-31). The dump + localization
were run on CPU (dump 3h36m vs ~15min GPU) because the GPU partition was 48/48 booked cluster-wide; the tiny
model makes CPU a viable fallback. Metrics in `results/loto/{baseline,orf_v2_attn}_onehot_holdout_Hepatocytes/`.

---

## Task 16: cross-study (human Ruiz-Orera) + cross-species (mouse Wang) held-out transfer (2026-07-14, LANDED 2026-07-16)

The strictest test of the profile's transferability: apply a Chothani-trained orf_v2_attn model,
without retraining, to two fully independent datasets and score the three evals (profile Pearson,
localization non-canonical AUROC, RiboCode drop-in) against each dataset's OWN RiboCode calls.

- **Cross-study (human):** Ruiz-Orera 2024 iPSC-cardiomyocyte (PRJEB65856; 5 Ribo + 5 RNA), a cell
  type absent from the 9-tissue training panel, same GENCODE v49 annotation. Reuses the Fibroblast
  universe, so all three backends (one-hot / Orthrus / RiNALMo) are scored. Pack: 36,668 tx, 33,185
  scorable; global-mean depth 743 (~10x shallower than Fibroblast -- the depth-normalized coverage
  input is what makes the transfer fair).
- **Cross-species (mouse):** Wang 2021 P42 liver (GSE94982; 2 Ribo + 2 RNA), GENCODE vM38. One-hot
  only, since the one-hot control (Task 15) needs no species-specific FM embeddings -- the human-trained
  one-hot model is applied directly to mouse sequence. Mouse universe 221,835 tx (65,936 pc + 155,899
  lncRNA), 41,096 scorable; global-mean depth 43.9. Mouse orf_track_v2 + tx2cds built from vM38.

Infrastructure built this session (methods.md 5K): `heldout/build_heldout_pack.py` (pool per-SRR ->
pack, by-id join, no 3.7 GB intermediates), `heldout/build_mouse_universe.py`, `--pack/--fasta` on
`build_orf_track.py`, `--gtf/--out` on `build_tx2cds.py`, `--heldout` mode on `dump_pred_profiles.py`
+ `eval_localization.py`, `heldout/eval_heldout_{predict,dropin}.sbatch`, `heldout/assemble_heldout.py`.
The trained model's `cov_norm=global_mean` uses each held-out pack's own mean depth (the transfer
mechanism). Both drop-in variants are reported: `pred_obsdepth` (predicted shape at real depth =
transferable-profile test) and `pred_preddepth` (fully standalone, also tests the depth-dependent count
head, expected to transfer worse cross-dataset).

Original CPU chain (35314728 ...) hit short-wall timeouts; re-run as self-contained GPU jobs
35335044-47 (see infrastructure note below). Results landed 2026-07-15.

### Cross-study (human Ruiz-Orera) -- LANDED

| backend | profile r (pc) | profile r (lnc) | drop-in F1 | precision | recall | novel-ORF recall | count-head F1 gap |
|---------|---------------|-----------------|-----------|-----------|--------|------------------|-------------------|
| one-hot | 0.4248        | 0.3350          | 0.931     | 0.937     | 0.925  | 0.797 (244/306)  | +0.002            |
| Orthrus | 0.4138        | 0.3626          | 0.925     | 0.917     | 0.932  | 0.801 (245/306)  | -0.001            |
| RiNALMo | 0.4178        | 0.3651          | 0.924     | 0.919     | 0.929  | 0.814 (249/306)  | -0.002            |

(n_pc = 32,400; n_lnc = 772 scorable held-out transcripts; drop-in over n_real = 12,871 real-density calls.)

- **The profile transfers.** One-hot pc profile Pearson 0.425 ties/beats the FMs (Orthrus 0.414, RiNALMo
  0.418); the one-hot >= FM finding (Task 15) holds on fully independent data. FMs edge one-hot only on
  lncRNA profile and novel-ORF recall (marginal).
- **Drop-in F1 ~0.93 cross-study**, one-hot best (0.931). Notably the predicted-profile F1 (0.93) EXCEEDS
  the real-held-out-profile-vs-official F1 (0.878, same for all backends): the model's denoised profile
  reproduces the deep official RiboCode calls BETTER than the shallow real Ruiz-Orera data (depth 743)
  does. The Ribo-seq refinement use case validates itself on independent data.
  **Checkpoint caveat (added 2026-08-08):** these three backends were run on
  `orf_v2_*_onehot_holdout_Hepatocytes`, not on the shipping union models. On the released models the
  0.931 is a `pred_obsdepth` number and holds (0.929 attn / 0.934 mamba4), but the *standalone*
  (`pred_preddepth`) arm falls to 0.876/0.880. See "Released-model re-run (2026-08-08)" below. The
  backend-comparison conclusion in this bullet is unaffected -- it is a within-checkpoint comparison.
- **Localization (full set, all 3 backends) -- LANDED 2026-07-16.** 827,551 candidate ORFs (10,071
  translated positives / 817,480 negatives; `min_orf_nt` 30, in-frame CDS collapsed). Predicted-frame0
  AUROC (does the predicted 3-nt periodicity discriminate a translated ORF from a candidate?), by ORF
  category, backends near-identical:

  | category | n_pos | pred_frame0 (onehot / orthrus / rinalmo) | obs_frame0 ceiling (Ruiz-Orera's own Ribo-seq) |
  |----------|-------|------------------------------------------|------------------------------------------------|
  | all          | 10,071 | 0.945 / 0.946 / **0.946** | 0.908 |
  | annotated    | 7,668  | 0.976 / 0.979 / 0.978     | 0.914 |
  | non-canonical| 2,403  | 0.847 / 0.838 / 0.845     | 0.887 |
  | uORF         | 1,463  | 0.909 / 0.905 / 0.911     | 0.894 |
  | novel (pred_density) | 184 | 0.924 / 0.935 / **0.936** | 0.642 |

  Three headline points. (1) **The predicted profile beats the study's own measured periodicity** for
  discriminating translated ORFs: all-category pred_frame0 0.945-0.946 EXCEEDS the obs_frame0 ceiling of
  0.908 computed from Ruiz-Orera's actual (shallow, ~10x under Chothani) Ribo-seq -- the "predict beats
  measure on shallow data" result made concrete on a fully independent human dataset. Same pattern on
  annotated (0.976-0.979 vs 0.914) and novel-lncRNA-type ORFs (density 0.92-0.94 vs 0.64). (2) The prelim
  0.752 non-canonical number was a small-shard artefact (n_pos 81); on the full set (n_pos 2,403) it is
  **0.847**, and it is the one category where the model sits just below the observed ceiling (0.887) rather
  than above it -- non-canonical ORFs are where measured signal, when you have it, still adds. (3) **Backend
  is irrelevant here:** one-hot, Orthrus, and RiNALMo agree to within 0.01 on every category, so the FM
  embeddings buy nothing over one-hot for cross-study localization (consistent with Task 15). Metrics:
  `results/heldout/human_ruizorera/{onehot,orthrus,rinalmo}/localization_metrics_human_ruizorera.json`.

### Count-head transferability -- RESOLVED

The `pred_obsdepth` vs `pred_preddepth` drop-in F1 gap is ~0 across all three backends (+0.002, -0.001,
-0.002). The predicted depth (count head, regressed to the Chothani reference) reproduces RiboCode calls
as well as the REAL observed depth. So the fully standalone `pred_preddepth` path (predicted shape x
predicted depth, needs no Ribo-seq at all) is validated, and the absolute count head transfers cross-study
for the ORF-calling use case. Consequence: the TE / per-million count-head refactor discussed in
`design_count_magnitude_transferability.md` is NOT needed for single-dataset deployment; it would only
help multi-dataset training coherence. This settles the transferability question empirically (the earlier
concern that the count head would not transfer was wrong for the calling use case).

### Cross-species (mouse Wang liver) -- LANDED (drop-in), 2026-07-15

The mouse genomic-key compare was blocked by a hardcoded human tx->gene map in
`compare_dropin_calls.py`: it always loaded the human v49 `tx2biotype.tsv`, so all 41,096 mouse ENSMUST
test tx missed the map, `keep_genes` came out empty, and `--key genomic` filtered out 100% of calls
("test genes: 0", n_real 0). Fixed 2026-07-15: added a `--tx2gene` argument (default = human) plus a loud
guard that aborts when 0 tx match the map (so this can't silently recur), built `data/mouse_tx2biotype.tsv`
from the vM38 GTF (278,326 tx, via the now-parameterized `build_tx2biotype.py --gtf/--out`), and reran.
Now `test genes: 11,617`.

> **SUPERSEDED 2026-08-08 -- see "Released-model re-run" immediately below.** The table and bullets in
> this subsection came from `orf_v2_attn_onehot_holdout_Hepatocytes`, a pre-nokozak / pre-mm1 /
> pre-union checkpoint, on the pre-decontamination Janich universe. They are kept as the historical
> record. The shape conclusion survived the re-run; the standalone-depth number did not.

Mouse drop-in, one-hot, genomic key (min_len 90, pred enrichment >= 0.5):

| comparison | nPred | nReal | F1 | precision | recall |
|---|---|---|---|---|---|
| real vs official (harness validation) | 13,579 | 13,992 | 0.981 | 0.995 | 0.966 |
| pred_obsdepth vs real (PRIMARY)       | 13,365 | 13,579 | 0.929 | 0.936 | 0.922 |
| pred_preddepth vs real (standalone)   | 13,914 | 13,579 | 0.919 | 0.908 | 0.930 |
| pred_obsdepth vs official             | 13,365 | 13,992 | 0.915 | 0.937 | 0.895 |

- **The human-trained model transfers across the species boundary.** Applied without retraining to mouse
  liver, the one-hot model recovers 92-93% of RiboCode's own mouse calls from predicted profiles alone
  (pred_obsdepth F1 0.929) -- essentially the same as the human cross-study F1 (0.931). It recovers 99.5% of
  annotated ORFs (11,897/11,962) and 76% of novel ORFs (440/579); the non-canonical strata are harder
  (uORF recall 0.208, dORF 0.120), the same pattern as human.
- **Count-head resolution holds cross-species.** pred_preddepth F1 0.919 vs pred_obsdepth 0.929 (gap
  -0.010) -- the fully standalone predicted-depth path is nearly as good as borrowing the real depth, even
  across the human->mouse boundary. This mirrors the human cross-study gap ~0 and extends the
  count-head-transfers conclusion to a second species.
- **Harness validated.** real_vs_official F1 0.981 confirms the mouse pipeline (vM38 RiboCode annotation +
  P-site calling) reproduces the official mouse calls near-perfectly, so the ~0.92 predicted numbers are a
  real model result, not a broken harness.
- Transcript-key cross-check (universe-sensitive, no gene map) is systematically lower as expected
  (real_vs_official 0.838, pred_obsdepth 0.816, pred_preddepth 0.772): the difference is the
  collapse-representative isoform reshuffle that genomic keying neutralizes, which is exactly why genomic is
  the headline key here and for human. Written to `dropin/transcript_key/dropin_metrics.json`.

Update (2026-07-16): the human full-set localization AUROC was reported in Task 16/17 (predicted-frame0
AUROC 0.945, beating the observed-profile ceiling 0.908). The cross-species headline became the drop-in
ORF-calling F1 (mouse Wang 0.929, Task 16/19) rather than a mouse profile-Pearson shape number; that
per-nt mouse profile-Pearson assembly remains the one un-collated secondary metric (low priority -- the
drop-in F1 is the stronger cross-species statement).

### Released-model re-run (2026-08-08): shape held, standalone depth did not

Every cross-study drop-in number quoted above and in the human cross-study section was measured on
`orf_v2_attn_onehot_holdout_Hepatocytes`. The project ships `*_union_noBrain_nokozak_mm1_*`. Directory
names encode the DATASET, not the CHECKPOINT, so the staleness was invisible until the npz `meta` field
was read. Re-dumped both released models over all four held-out sets
(`scripts/dump_liver_released.sbatch`, `scripts/ruizorera_released_dumps.sbatch`) and re-scored through
the same `compare_dropin_calls.py`, genomic key.

| dataset | arm | attn F1 | mamba4 F1 | previously reported |
|---|---|--:|--:|--:|
| Ruiz-Orera (human iPSC-CM) | pred_obsdepth  | 0.929 | 0.934 | 0.931 |
| Ruiz-Orera (human iPSC-CM) | pred_preddepth | 0.876 | 0.880 | 0.930 |
| Wang liver (mouse)         | pred_obsdepth  | 0.923 | 0.927 | 0.929 |
| Wang liver (mouse)         | pred_preddepth | 0.867 | 0.867 | 0.919 |
| Janich liver (mouse)       | pred_obsdepth  | 0.919 | 0.926 | not run on released |
| Janich liver (mouse)       | pred_preddepth | 0.898 | 0.893 | not run on released |
| GSE243134 liver (mouse)    | pred_obsdepth  | 0.919 | 0.921 | not run on released |
| GSE243134 liver (mouse)    | pred_preddepth | 0.895 | 0.896 | not run on released |

- **The shape claim is intact.** `pred_obsdepth` moves by at most 0.006 in either direction on both the
  human and the mouse set that had a prior number. The profile the model predicts is as good on the
  released checkpoints as it was on the old one, and it is now confirmed on two additional independent
  mouse-liver studies that were never scored before.
- **The standalone claim was overstated.** `pred_preddepth` drops 0.052-0.054 on the two datasets with a
  prior number. The mechanism is visible in the counts: the standalone arm predicts MORE ORFs than the
  reference has (Ruiz-Orera 14,641 vs 13,147 real), so the loss is entirely precision (0.942 -> 0.832 on
  Ruiz-Orera) while recall actually rises. The count head over-calls on the union universe -- expected,
  since union is 84,472 tx against a much smaller earlier universe, and the extra transcripts are
  low-expression ones where predicted depth is least constrained.
- **"Count-head resolution holds cross-species" needs qualifying.** The obsdepth-to-preddepth gap is no
  longer ~0.010; it is -0.024 to -0.062 (mean -0.039 across the six mouse arms). The count head still
  transfers, but it costs real precision, and the Ribo-seq-free Poisson calibration (below) is the lever
  for it rather than an optional refinement.
- **attn and mamba4 are interchangeable on every one of these**, max separation 0.007 F1. Consistent with
  the 3-seed union comparison.

Per-dataset tables with class-level recall, and the checkpoint string behind every row, are generated
into the tutorial by `tutorial/make_heldout_{human,mouse}.py`.

### Infrastructure note (the timeout re-run)

The original held-out chain sharded the predicted-profile dump across the `short` partition (1 h wall)
with afterok-chained finalize. RiNALMo's O(L^2) transformer attention on long transcripts blew the 55-min
wall (6/8 shards TIMEOUT), and `afterok` let a single timed-out shard poison the whole finalize (all four
finalize jobs went DependencyNeverSatisfied). Re-run as one self-contained GPU job per (dataset, backend)
(`heldout/dump_finalize_gpu.sbatch`): the un-sharded dump on a GPU finishes in minutes and writes
`pred_profiles.npz` directly (no merge), then runs the 3 drop-in density variants + compare in the same
job. This removed all array/merge/afterok fragility. Jobs 35335044-47.


## Task 17: mixer + input-modality + capacity ablations, 9-fold LOTO, Mamba variant (2026-07-15)

Controlled single-variable studies off the `orf_v2_attn_onehot` config, all on the SAME Hepatocytes
hold-out (except 17D, which varies the held-out tissue). Metric below is the CONVERGED held-out
TEST-set score (median per-transcript profile Pearson and median predicted-profile 3-nt periodicity
over the 33,918 test transcripts = 32,869 protein-coding + 1,049 lncRNA), read from each run's
`test_metrics.json` -- not the mid-flight validation Pearson used in earlier drafts of this section.
Reference baseline `orf_v2_attn_onehot` (transformer, both inputs, 2 attn layers, ch 256): protein-coding
Pearson **0.639**, periodicity **+0.391**, lncRNA Pearson **0.447**. 17A/17B/17D-landed are converged;
17C (Mamba) and the remaining 17D folds are still marked [PRELIM].

### 17A. Input-modality ablation -- RNA-seq drives MAGNITUDE, sequence drives SHAPE (LANDED)

| input_mode | pc Pearson | pc periodicity | lncRNA Pearson |
|------------|------------|----------------|----------------|
| both (baseline)        | 0.639 | +0.391 | 0.447 |
| emb one-hot (seq only) | 0.579 | +0.471 | 0.390 |
| emb RiNALMo (seq only) | 0.546 | +0.487 | 0.346 |
| cov (RNA-seq only)     | 0.122 | -0.011 | 0.088 |

The profile head (shape) is a scale-free multinomial and the count head (magnitude) reads
log1p(total coverage), so the architecture was designed for RNA-seq -> magnitude and sequence -> shape.
This ablation confirms that split cross-tissue for the first time (previously only within-tissue, fold 0):

- Sequence-only recovers 91% of the both-modality profile Pearson (0.579 of 0.639) and REPRODUCES the
  3-nt periodicity -- in fact sharpens it (+0.471 vs the baseline's +0.391), because the coverage channel
  contributes a smooth envelope that slightly dilutes fine periodicity.
- RNA-seq-only collapses to 0.122 with periodicity DESTROYED (-0.011): a smooth magnitude signal cannot
  express the periodic shape at all.
- Adding RNA-seq on top of sequence lifts whole-transcript pc Pearson by +0.060 (0.579 -> 0.639) -- it
  marks where coverage exists and firms up the broad envelope -- but the periodicity itself is sequence-borne.

So the shape (and all of the periodicity RiboCode reads) comes from sequence; RNA-seq supplies magnitude
plus a modest whole-transcript-envelope lift. This is the mirror image of Translatomer's abundance-weighted
envelope metric (RNA-seq-only 0.731, sequence +5.3%); the difference is the metric (within-CDS periodic
shape here vs a 65 kb genomic-window envelope there). one-hot >= RiNALMo (0.579 vs 0.546), consistent with
the FM-no-lift finding in Task 15.

### 17B. Capacity -- attention depth trades whole-tx Pearson for periodicity; width is neutral (LANDED)

| config | pc Pearson | pc periodicity | lncRNA Pearson |
|--------|------------|----------------|----------------|
| baseline (2 attn layers, ch 256) | 0.639 | +0.391 | 0.447 |
| attn = 4 layers                  | 0.603 | +0.430 | 0.423 |
| channels = 384                   | 0.633 | +0.418 | 0.441 |

On convergence, NEITHER deeper nor wider beats the baseline on whole-transcript pc Pearson. Deeper
(4 attn layers) actually HURTS Pearson (0.603 < 0.639) while sharpening periodicity (+0.430 > +0.391) --
a real depth-vs-periodicity tradeoff, not a free lift. Wider (384 ch) is within noise on Pearson
(0.633 ~ 0.639) and slightly better on periodicity. Net: the baseline 2-attn/256-ch body is at the
capacity sweet spot for whole-tx shape. (This REVISES the earlier mid-flight read of this section, which
had reported "attn=4 gives a small lift 0.502 > 0.495" off a noisy validation Pearson; the converged
held-out test metric reverses the Pearson ordering.)

**Does the +0.039 periodicity help ORF calling? Drop-in test of attn=4 (2026-07-15, LANDED).** The
4-layer variant was pushed through the identical RiboCode drop-in chain (dump -> detectORF x3 -> compare,
genomic key) -- a rerun of the exact baseline drop-in with ONLY depth changed (jobs 35357991-993):

| metric (pred_obsdepth vs real) | baseline (2L) | attn=4 (4L) | delta |
|---|---|---|---|
| drop-in F1              | 0.922 | 0.904 | -0.018 |
| precision              | 0.913 | 0.871 | -0.042 |
| recall                 | 0.932 | 0.940 | +0.008 |
| non-canonical recall   | 0.614 | 0.678 | +0.064 |
| uORF recall            | 0.689 | 0.809 | +0.120 |
| Overlap_uORF recall    | 0.673 | 0.732 | +0.059 |

The periodicity gain IS real for calling, and it lands exactly where predicted -- the non-canonical strata.
Deeper attention recovers +0.064 more non-canonical ORFs overall and +0.120 more uORFs, because uORF / dORF
/ novel calls hinge on clean 3-nt periodicity rather than the whole-tx envelope, so the sharper periodicity
pays off there. BUT it costs precision (-0.042: the deeper model also over-calls more), so the AGGREGATE F1
is slightly LOWER (0.904 vs 0.922). Verdict: whole-tx Pearson is a good proxy for aggregate F1 (both favor
2-layer), but it MASKS a real recall/precision tradeoff on the non-canonical strata. The 4-layer is the
better uORF / non-canonical DISCOVERY model (recover more, then filter on precision); the 2-layer is the
better-calibrated general caller. So "2-layer is the capacity sweet spot" holds for aggregate F1 and
whole-tx shape, but should be QUALIFIED: for uORF / non-canonical discovery specifically, deeper attention
is worth the precision cost.

### 17C. Mixer swap -- bidirectional Mamba (adapted from the seq2ribo polisher) (LANDED)

A bidirectional Mamba body (`BiMambaBlock` = forward + reversed Mamba, summed; `MambaBody` stack) replaces
the 2 transformer layers -- a controlled single-variable swap, trained in the Orthrus SIF (mamba_ssm
1.2.0.post1). Adapted from seq2ribo's polisher (Kaynar & Kingsford 2026), made BIDIRECTIONAL because that
model leans on its sTASEP simulation prior for downstream context and this model has no simulator (see
`manuscript/related_work.md`). O(L) vs the transformer's O(L^2). `mamba_onehot` (5,768,962 params) CONVERGED
(job 35349301, early-stopped near epoch 9). Held-out Hepatocytes test-set medians vs the transformer anchor:

| mixer (2 layers, ch 256, both inputs, one-hot) | pc Pearson | pc periodicity | lncRNA Pearson |
|---|---|---|---|
| Transformer (baseline)      | 0.639 | +0.391 | 0.447 |
| bidirectional Mamba         | 0.599 | +0.458 | 0.413 |

The bidirectional Mamba lands BELOW the transformer on whole-tx pc Pearson (0.599 vs 0.639, -0.040) but
ABOVE it on periodicity (+0.458 vs +0.391) -- the same Pearson-for-periodicity tradeoff as deeper attention
(17B: attn=4 gave 0.603 / +0.430). It is competitive, not better, on this hold-out; the mid-flight
validation read predicted exactly this (Mamba plateaued at val pc Pearson 0.492 @ epoch 3 vs the
transformer's converged 0.505 @ epoch 9 on the same val split). The O(L) vs O(L^2) scaling is the standing
structural advantage, so the Mamba mixer is worth keeping as an option for very long transcripts where the
transformer's O(L^2) attention is the inference bottleneck, but it is not the default. (The 1-epoch smoke's
periodicity +0.623 was an undertraining artifact; it settled to +0.458.) Code: `model.py`
BiMambaBlock/MambaBody + `train_loto.py --mixer mamba`; all
model-reconstruction sites (eval_extra, dump_pred_profiles, eval_localization, plot_example_uorf) read
`mixer` from `args.json` (transformer default, backward compatible). Smoke run confirmed the CUDA
forward/backward/eval work end-to-end.

**Mamba drop-in test (2026-07-16): the periodicity-for-precision tradeoff is a GENERAL property of
stronger mixers, not a transformer-depth quirk.** The Mamba variant was pushed through the identical
RiboCode drop-in chain (dump in the Orthrus SIF -> detectORF x3 -> compare). It shows the SAME signature as
attn=4 (17B) -- two independent "stronger global mixer" architectures, same result:

| pred_obsdepth vs real | baseline (2L transf) | attn=4 | Mamba |
|---|---|---|---|
| drop-in F1            | 0.922 | 0.904 | 0.906 |
| precision            | 0.913 | 0.871 | 0.880 |
| recall               | 0.932 | 0.940 | 0.933 |
| non-canonical recall | 0.614 | 0.678 | 0.638 |
| uORF recall          | 0.689 | 0.809 | 0.761 |

Both the deeper-attention and the state-space mixer trade aggregate F1 (via precision) for non-canonical /
uORF recall relative to the 2-layer baseline: Mamba lifts uORF recall +0.072 (0.689 -> 0.761) and
non-canonical +0.024, at precision -0.033 and F1 -0.016. So "sharper periodicity from a stronger mixer
buys non-canonical ORF recovery at a precision cost" now holds across two unrelated architectures -- a
general property of the profile model, not a one-off. The 2-layer transformer remains the best-calibrated
general caller; deeper-attention or Mamba are the better uORF / non-canonical DISCOVERY front-ends.

### 17D. Full 9-tissue leave-one-tissue-out (one-hot) (LANDED -- 9/9 folds)

Each of the 9 Chothani tissues held out once (previously only Hepatocytes), one-hot `orf_v2_attn`,
converged held-out test-set medians. The `period_OBS` column is the 3-nt periodicity of the REAL held-out
Ribo-seq -- i.e. the target's own quality (independent of the model):

| held-out tissue | pc Pearson | pc period_pred | pc period_OBS | lncRNA Pearson |
|-----------------|------------|----------------|---------------|----------------|
| Hepatocytes | 0.638 | 0.380 | 0.158 | 0.444 |
| Fibroblast  | 0.629 | 0.340 | 0.222 | 0.473 |
| HCAEC       | 0.570 | 0.415 | 0.171 | 0.436 |
| Fat         | 0.562 | 0.299 | 0.260 | 0.416 |
| ES          | 0.561 | 0.375 | 0.202 | 0.446 |
| VSMC        | 0.490 | 0.442 | 0.110 | 0.352 |
| HA_EC       | 0.447 | 0.452 | 0.098 | 0.338 |
| HUVEC       | 0.413 | 0.414 | 0.110 | 0.347 |
| Brain       | 0.234 | 0.431 | 0.044 | 0.153 |

Mean pc Pearson over all 9 folds is **0.505** (**0.539** excluding the Brain outlier). The headline:
**the cross-tissue pc-Pearson spread (0.23 to 0.64) is largely a held-out-TARGET data-quality ceiling, not
a model-generalization gradient.** The model's own output quality (`period_pred`) is roughly constant
across tissues (0.30 to 0.45 -- it predicts a cleanly periodic profile everywhere); what varies is how
periodic the OBSERVED held-out Ribo-seq is (`period_OBS` 0.044 to 0.260), and pc Pearson broadly tracks
`period_OBS` (strong at the low end, noisier at the top where depth also matters -- Hepatocytes is deep so
it scores high at moderate OBS). Brain is the extreme: observed periodicity 0.044 (near noise -- the
shallowest, noisiest tissue in the panel), so even a good periodic prediction cannot correlate highly
against an essentially aperiodic target, capping pc Pearson at 0.234. This is the same shallow-data ceiling
seen cross-study in Task 16: on the full-set Ruiz-Orera localization the predicted-frame0 AUROC (0.945 all,
0.976 annotated) actually EXCEEDS that study's own observed-periodicity ceiling (0.908 / 0.914), because the
external Ribo-seq is ~10x shallower than the model's denoised profile. So the LOTO result reads better as
"the model predicts a consistent periodic profile across every held-out tissue; the score reflects target
quality" than as "the model generalizes unevenly across tissues."

Jobs: all 9 LOTO folds + ablations + mamba converged (35348241 A100 + 35348242 A5500). Anchor baseline:
`results/loto/orf_v2_attn_onehot_holdout_Hepatocytes`.


## Task 18: depth crossover -- at what sequencing depth does PREDICTING beat MEASURING? (2026-07-16, LANDED)

Makes the Task 17D "profile-Pearson dispersion is a target-quality ceiling" claim concrete and
quantitative, and directly demonstrates the Task 16 denoising finding: below a crossover depth, predicting
the Ribo-seq profile from RNA-seq (no experiment) recovers MORE of the deep ORF-call truth than a real
experiment at that depth. Hepatocytes held out (`orf_v2_attn_onehot`); full test-set depth 6.30e8 P-sites
(33,918 tx). Figure: `figures/depth_crossover/depth_crossover.png` (+ `FIGURE_DATA_INPUTS.md`).

Method: binomially thin the REAL Hepatocytes per-nt P-site profile to a fraction f of depth (each footprint
kept i.i.d. w.p. f; `ribocode_dropin.py --subsample f --seed s`, `subsample_depth.sbatch` array), run
RiboCode `detectORF` on the thinned counts, and score the calls against the deep OFFICIAL Hepatocytes calls
(genomic locus key, 90 nt, 0.5x enrichment). The measurement curve F1(f) is compared to the depth-INDEPENDENT
prediction line `pred_preddepth` (predicted shape x count-head predicted depth; uses RNA-seq + sequence,
ZERO Ribo-seq), also vs the deep official. 8 depths x 3 seeds + a low-depth tail (`subsample_depth_curve.py`).

| test-set P-sites (f) | measured F1 | measured non-canon recall |
|---|---|---|
| 6.30e8 (1.0)   | 0.881 | 0.57 |
| 2.21e8 (0.35)  | 0.870 | 0.51 |
| 1.26e8 (0.2)   | 0.860 | 0.46 |
| 6.30e7 (0.1)   | 0.844 | 0.37 |
| 3.15e7 (0.05)  | 0.821 | (~0.32) |
| 2.21e7 (0.035) | 0.805 | 0.28 |
| 1.26e7 (0.02)  | 0.779 | 0.23 |
| 6.30e6 (0.01)  | 0.733 | 0.16 |
| 3.15e6 (0.005) | 0.671 | 0.09 |
| 1.26e6 (0.002) | 0.545 | 0.04 |
| 3.15e5 (0.0005)| 0.271 | 0.00 |
| **predict, NO Ribo-seq** | **0.818** | **0.323** |

**Crossover: ~2.9e7 test-set P-sites (f~0.047).** Below ~29 million P-sites on these transcripts,
predicting the profile from RNA-seq alone recovers more of the deep ORF-call truth (F1 0.818) than
measuring Ribo-seq at that depth. At 6.3e6 P-sites (a normal-depth library on this test set) measuring
gives F1 0.733 vs predicting 0.818 -- an 8.5-point gap in the prediction's favor; at 3.15e6 it is 0.671
vs 0.818. The non-canonical panel crosses at a similar depth (~4e7). (`pred_obsdepth`, predicted shape at
full real depth, sits at 0.819 -- essentially on the standalone line, confirming the count head's predicted
depth is as good as the real depth here, consistent with Task 16.)

This is the concrete form of the ceiling argument: a profile-Pearson or F1 number measured against shallow
Ribo-seq is not a model ceiling, it is a MEASUREMENT ceiling -- and below ~29M P-sites the measurement is
the weaker of the two estimates of the deep truth. Refined grid pinned the crossover at 29M (the coarse
pilot's 56M was interpolated across a gap; the curve is convex near full depth). Jobs 35360209 (pilot),
35364855 (refined 8x3).


## Task 19: non-AUG (ATG + CTG) drop-in -- does the model help call CTG-initiated ORFs, or only ATG? (2026-07-16, LANDED)

> **CHECKPOINT NOTE (added 2026-08-08).** This whole section is on
> `orf_v2_*_onehot_holdout_Hepatocytes` and the ATG+CTG ORF track -- pre-nokozak / pre-mm1 /
> pre-union. It was not re-run: the CTG track is a separate ORF track from the one the released
> models trained on (`--kozak none`, ATG-only), so re-running is a new experiment rather than a
> refresh, and it sits on the OPEN-elective list. The comparisons here are internally consistent
> (all arms share one checkpoint), so the conclusions stand as stated. Do not quote the absolute
> numbers -- including the "ATG F1 0.929" reference point below -- next to released-model numbers;
> on the released models the corresponding Wang figure is 0.923 (attn) / 0.927 (mamba4).

Every RiboCode call in Tasks 14-18 was ATG-only (verified: 100.0% of real and predicted calls, all
categories incl. the non-annotated uORF/dORF/novel classes, start with ATG). So the "non-canonical"
categories test non-canonical POSITION, all ATG-initiated; the v2 non-AUG-graded start channel had only
ever been validated on profile SHAPE (uORF Pearson, Task 9c), never on non-AUG DETECTION. This task
re-runs the Hepatocytes held-out drop-in (one-hot `orf_v2_attn`, reusing the existing `pred_profiles.npz`,
no re-dump) with RiboCode's `ALTERNATIVE_START_CODON_LIST=["CTG"]` (== the `RiboCode -A CTG` CLI) on all
three density variants, so real-vs-predicted ORF calling can be compared on CTG ORFs. CTG is FALLBACK-ONLY
in RiboCode (`orf_finder.orf_find`: an alt-start opens an ORF only where its in-frame stop has no ATG), so
CTG calls are a DISJOINT addition on top of the unchanged ATG calls -- confirmed empirically below.

**Start-codon composition of the call sets** (genomic loci, 90 nt / raw-pval / 0.5x-enrichment filters):

| set | n calls | ATG | CTG |
|-----|--------:|----:|----:|
| real (CTG-aware truth) | 14,622 | 13,141 (89.9%) | 1,481 (10.1%) |
| pred_obsdepth | 14,416 | 13,423 (93.1%) | 993 (6.9%) |
| pred_preddepth (standalone) | 13,846 | 13,087 (94.5%) | 759 (5.5%) |

Real CTG ORFs are **92% uORFs** (920 uORF + 446 Overlap_uORF of 1,481; then 59 novel, 26 dORF, 19 internal,
8 Overlap_dORF, 3 annotated) -- exactly where near-cognate CUG initiation is expected biologically.

**Predicted-vs-real ORF calling, stratified by start codon** (genomic-locus key, same filters as Task 14):

| start | variant | n real | n pred | match | precision | recall | F1 |
|-------|---------|-------:|-------:|------:|----------:|-------:|---:|
| ATG | pred_obsdepth | 13,141 | 13,423 | 12,257 | 0.913 | 0.933 | **0.923** |
| ATG | pred_preddepth | 13,141 | 13,087 | 12,103 | 0.925 | 0.921 | **0.923** |
| CTG | pred_obsdepth | 1,481 | 993 | 507 | 0.511 | 0.342 | **0.410** |
| CTG | pred_preddepth | 1,481 | 759 | 420 | 0.553 | 0.284 | **0.375** |

Findings:
1. **Method is sound.** The stratified ATG F1 0.923 reproduces the ATG-only drop-in exactly, so enabling
   CTG did NOT perturb the ATG calls -- the fallback-only semantics hold. And of the 1,481 real CTG loci,
   the model recovered 521 at the locus (locus-recall 0.352) of which 507 it also called as CTG, so there
   is essentially no codon-swapping: the model genuinely misses the other ~2/3, it does not mislabel them.
2. **CTG detection is real but much weaker than ATG: F1 0.41 vs 0.92.** The predicted profile recovers
   only ~1/3 of real CTG ORFs (recall 0.342) at ~half precision (0.511). Still, 507 recovered CTG ORFs
   (almost all uORFs) that the model was never explicitly supervised to call is the FIRST direct evidence
   the pipeline has any non-AUG detection ability -- the detection-level correlate of the v2 start
   channel's uORF-Pearson gain (Task 9c).
3. **The model UNDER-proposes CTG:** 993 CTG calls vs 1,481 real (6.9% vs 10.1% of its calls). It is
   conservative on near-cognate starts, consistent with being trained on a periodicity-filtered P-site
   target where CTG-uORF signal is weaker and noisier than CDS signal.
4. **The CTG limitation is shape/detection, not depth.** Standalone `pred_preddepth` CTG F1 0.375 is close
   to `pred_obsdepth` 0.410, so the count head's predicted depth is not the bottleneck (consistent with
   Tasks 16/18); the gap is whether the predicted profile concentrates enough in-frame signal at the CTG.
5. **Deployment implication.** Adding CTG drags the COMBINED drop-in F1 to 0.880 (from ATG-only 0.923),
   because the noisy CTG tail (~10% of calls at F1 0.41) dilutes it. So keep RiboCode ATG-only when the
   target is canonical ORFs; enable CTG only when non-AUG uORFs are specifically wanted and the softer
   truth is acceptable -- it recovers ~500 real CTG uORFs at ~half precision.

CAVEAT: RiboCode's CTG calls are a SOFTER TRUTH than ATG (near-cognate initiation is inherently noisier
and less validated), so recall 0.34 is partly the truth being less reproducible, not only model miss.
Scripts: `ribocode_dropin.py --alt_start_codons CTG`, `ribocode_dropin_ctg.sbatch`, `compare_dropin_ctg.py`
(start-codon lookup from the RiboCode annotation FASTA), `ctg_compare.sbatch`. Metrics:
`results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/dropin_ctg/{dropin_metrics,dropin_ctg_metrics}.json`.
Jobs 35368963 (drop-in array), 35368993 (compare).

### Does the CTG result hold up elsewhere? -- 9 runs across backend / dataset / species / architecture (2026-07-16)

Re-ran the ATG+CTG drop-in on 8 more runs that already had a dumped `pred_profiles.npz` (no re-dump),
spanning four axes. `pred_obsdepth` vs real, genomic-locus key, same filters. F1saC = standalone
(`pred_preddepth`) CTG F1; %CTGr = CTG share of the real calls; %uORF = uORF share of real CTG ORFs.

| run | ATG F1 | nRealC | nPredC | match | CTG prec | CTG rec | CTG F1 | F1saC | %CTGr | %uORF |
|-----|-------:|-------:|-------:|------:|---------:|--------:|-------:|------:|------:|------:|
| Hep onehot (baseline)      | 0.923 | 1481 |  993 | 507 | 0.511 | 0.342 | **0.410** | 0.375 | 10.1% | 92% |
| Hep orthrus                | 0.914 | 1481 | 1122 | 490 | 0.437 | 0.331 | 0.376 | 0.360 | 10.1% | 92% |
| Hep rinalmo                | 0.911 | 1481 | 1160 | 543 | 0.468 | 0.367 | 0.411 | 0.378 | 10.1% | 92% |
| Ruiz-Orera onehot          | 0.931 | 1713 |  687 | 374 | 0.544 | 0.218 | 0.312 | 0.262 | 11.8% | 93% |
| Ruiz-Orera orthrus         | 0.925 | 1713 |  828 | 411 | 0.496 | 0.240 | 0.323 | 0.305 | 11.8% | 93% |
| Ruiz-Orera rinalmo         | 0.924 | 1713 |  904 | 483 | 0.534 | 0.282 | 0.369 | 0.336 | 11.8% | 93% |
| attn4 onehot (Hep)         | 0.905 | 1481 | 1831 | 691 | 0.377 | 0.467 | 0.417 | 0.403 | 10.1% | 92% |
| mamba onehot (Hep)         | 0.906 | 1481 | 1223 | 547 | 0.447 | 0.369 | 0.405 | 0.360 | 10.1% | 92% |
| mouse Wang liver onehot    | 0.929 |  606 |  233 |  49 | 0.210 | 0.081 | 0.117 | 0.136 |  4.3% | 84% |

**ATG F1 stays 0.905-0.931 everywhere** -- the sanity check holds across all 9 (enabling CTG never perturbs
the ATG calls). **CTG F1 ranges 0.117-0.417 (mean 0.349)**, i.e. CTG detection is real-but-weak everywhere
and much below the ~0.92 ATG line. By axis:

1. **Backend is irrelevant for CTG too:** within-Hepatocytes CTG F1 0.376-0.411 (onehot ~ rinalmo > orthrus
   by a hair), the same backend-agnosticism seen on every other metric.
2. **Independent human dataset (Ruiz-Orera) drops CTG F1 to 0.31-0.37, entirely via RECALL** (0.22-0.28 vs
   0.33-0.37 within-tissue) while precision HOLDS (~0.50-0.54). Ruiz-Orera is ~10x shallower, so the real
   CTG truth is sparser/noisier and the predicted profile has less to lock onto -- the model still proposes
   precise CTG calls, just fewer of them. Composition is stable (11.8% CTG, 93% uORF).
3. **Cross-species mouse collapses to CTG F1 0.117 (recall 0.081, 49/606).** The hard triple of cross-species
   + shallow depth (mouse Wang mean depth ~44, the shallowest) + a smaller/sparser CTG truth (only 4.3% of
   calls, 606 loci). Canonical transfer is fine (ATG F1 0.929), but non-AUG uORFs are less positionally
   conserved and the shallow noisy target leaves almost nothing to recover. Small-n, treat as a floor.
4. **A stronger mixer lifts CTG RECALL, the same periodicity-for-precision tradeoff seen on ATG uORFs.**
   attn4 pushes CTG recall 0.342 -> **0.467** (+0.125) by proposing far more CTG ORFs (1831 vs 993, now
   OVER-proposing vs 1481 real) and recovering more (691 vs 507), but precision falls 0.511 -> 0.377, so net
   CTG F1 barely moves (0.417). mamba is milder (recall 0.369). This is the first evidence that the
   architectures which boosted uORF recall (Tasks 9/17) improve non-AUG DETECTION recall specifically -- the
   v2 start-channel + attention synergy operates on CTG ORFs, but the extra sensitivity is bought with false
   positives, not free.

Bottom line: the ~0.41 CTG F1 is a within-tissue, deep-target ceiling. It is backend-agnostic, degrades
under distribution shift (independent dataset, and especially cross-species) primarily through RECALL, and
is uniformly ~0.5x the ATG F1. The model's non-AUG ability is genuine but fragile, and a stronger mixer
trades precision to raise its recall. All caveats from the Hepatocytes result stand, amplified on the
shallow held-outs: RiboCode's CTG truth is soft, and shallower data makes it softer, so the recall drops
are part model, part truth quality. Aggregate: `results/ctg_across_runs.json`; per-run
`results/<run>/dropin_ctg/dropin_ctg_metrics.json`; `aggregate_ctg.py`. Jobs 35369472-35369487.

### Task 19b: all-near-cognate drop-in (2026-07-17) -- RiboCode's fallback pre-emption collapses it to CTG

Extended the drop-in to ALL 9 near-cognates (`--alt_start_codons CTG,GTG,TTG,ACG,ATA,ATT,ATC,AAG,AGG`,
`ribocode_dropin_ctg.sbatch OUTSUB=dropin_allalt`; new `compare_dropin_allalt.py` stratifies by codon) on
the base / mamba / attn4 Hepatocytes models. Result: the real call set is ONLY ATG (13,141) + CTG (1,481)
-- **no GTG/TTG/ACG/etc. at all**. RiboCode's `ALTERNATIVE_START_CODON_LIST` is fallback-only AND
priority-ordered, so CTG (first in the list) pre-empts every other near-cognate; adding the other 8 codons
yields essentially zero additional calls. So "all alt ORFs" via RiboCode == ATG+CTG in practice, and the
per-codon numbers reproduce Task 19 exactly (base ATG F1 0.923 / CTG 0.410; attn4 lifts CTG recall
0.342->0.467, F1 0.417; mamba CTG 0.405). **Implication for the proteogenomics DB (proteogenomics/):** the
model-selected novel-ORF search database CANNOT be built from RiboCode's alt-start caller (it will only
ever surface CTG). The DB must instead ENUMERATE candidate ORFs directly (all near-cognate start ..
in-frame stop) and score each with the model's predicted profile -- the caller is the bottleneck, not the
model. Jobs 35499198-35499200 (drop-in) + 35499374-35499376 (compare); metrics in each run's
`dropin_allalt/dropin_allalt_metrics.json`.

## Task 20: Kozak start-context ablation -- remove / empirical / learned (2026-07-16, LANDED)

Motivated by disliking the hand-picked Kozak heuristic in the ORF-track start channel
(`kz = 0.5*[purine at -3] + 0.5*[G at +4]`). Four one-hot arms differing ONLY in start-propensity
channel 3 (occupancy channels 0-2 + stop channel 4 byte-identical), all trained on Fibroblast and tested
on held-out Hepatocytes (`orf_v2_attn` config, identical protocol: 40 epochs, patience 8, budget 16000),
answer three questions: **Q1** remove the heuristic (V1 no-Kozak), **Q2** replace it with an
empirically-fit PWM (V2), **Q3** let the model learn the context weights (V3, a `Conv1d(4->1, k=10)` gate
over the one-hot slice). Full setup in methods.md 5N; plan + pre-registered prediction in `KOZAK_PLAN.md`.
Consolidated `results/kozak/kozak_summary.{md,json}`; per-arm under `results/kozak/<arm>_onehot_fib2hep/`.

| metric | V0 heuristic | V1 no-Kozak | V2 empirical | V3 learned |
|--------|:---:|:---:|:---:|:---:|
| **CTG non-AUG drop-in F1** (pred profile @ obs depth) | 0.340 | **0.369** | 0.343 | 0.320 |
| CTG recall | 0.296 | 0.302 | 0.287 | 0.238 |
| CTG precision | 0.401 | 0.475 | 0.428 | **0.489** |
| n pred CTG loci (real = 1481) | 1092 | 942 | 994 | 720 |
| non-canonical AUROC (len-ctrl) | 0.839 | **0.848** | 0.836 | 0.833 |
| uORF AUROC (len-ctrl) | 0.921 | 0.923 | 0.917 | 0.916 |
| uORF frame0 Pearson (fidelity) | 0.461 | **0.503** | 0.460 | 0.472 |
| ATG drop-in F1 | 0.907 | 0.923 | 0.906 | 0.914 |
| annotated-ORF AUROC (len-ctrl) | 0.912 | 0.912 | 0.910 | 0.915 |
| pc whole-tx profile Pearson | 0.555 | **0.587** | 0.551 | 0.565 |
| pc 5'UTR profile Pearson | 0.501 | **0.536** | 0.494 | 0.503 |
| **best Fibroblast-val Pearson** (in-distribution) | 0.6458 | 0.6418 | 0.6431 | 0.6424 |

**Interpretive anchor -- the four arms are in-distribution equivalent.** Best Fibroblast validation
Pearson ties at 0.6418-0.6458 (spread 0.004), and V1 (no-Kozak) is actually the *lowest* there. So none of
the cross-tissue test differences below come from one arm being a better-fit model; they are pure
held-out-tissue generalization on a soft-truth task, single seed per arm. Read the deltas as directional,
not decisive.

**Q1 (remove the heuristic): removing it is free, and directionally best.** V1 tops nearly every held-out
metric -- CTG F1 0.369 vs 0.340 (+0.029, via precision 0.475 vs 0.401 at equal recall), uORF frame0
Pearson 0.503 vs 0.461, non-canonical AUROC 0.848, whole-pc 0.587, 5'UTR 0.536, even ATG F1 0.923 (the
Kozak factor multiplies ATG starts too, so dropping it un-penalizes ATG starts with poor context). The
consistency of the V1 win across ~7 semi-independent metrics argues it is a real, if small, effect rather
than noise. The heuristic was a mild net *negative*: it injects a canonical-CDS prior that mildly
mis-weights the non-canonical starts the project targets.

**Q2 (empirical PWM): does not rescue it.** V2 CTG F1 0.343 == V0 0.340 (+0.003), uORF AUROC 0.917 vs
0.921 -- a wash, and below V1. Fitting the Kozak context *better* does not help, because the problem is not
the fit quality but that imposing a canonical-CDS start-context prior on ~92%-uORF CTG starts is the wrong
move (pre-registered prediction confirmed: the CDS PWM floors 15-25% of uORFs to 0.5x; see `KOZAK_PLAN.md`
and `kozak_context_alt_orfs.py`).

**Q3 (learn it): the model CAN, which is exactly why the explicit gate is redundant.** V3's learned 4x10
kernel (`figures/kozak/learned_vs_empirical_vs_heuristic.png`), never shown the annotation, rediscovers
the canonical -3 purine (A+G weight +0.759 vs C+T -0.808), correctly down-weights the weak +4 position
(G +0.128 vs +0.029 others; empirical +4 G is only 51.5%), and **matches the empirical PWM at Pearson r =
0.807** over the 7 context positions -- two independent routes (Ribo-seq supervision vs annotated-start
counting) converging on the same matrix. But *behaviorally* the learned gate is the WORST detector: it is
the most conservative (720 CTG calls, precision 0.489 but recall 0.238), netting the lowest CTG F1 (0.320).
The one-hot backbone already carries -3/+4 in its receptive field and uses start context implicitly (that
is why V1, with no explicit start-context channel at all, does best), so bolting on an *explicit* Kozak
gate -- heuristic, empirical, or learned -- is at best redundant and at worst (the hard sigmoid gate)
suppresses recall.

**Synthesis + recommendation.** The Kozak start-context factor does essentially no useful work in this
model. In-distribution all four arms tie (val 0.642-0.646); on cross-tissue non-canonical detection the
only directional signal is that *removing* the explicit factor (V1) is never worse and modestly best. The
mechanistic reason is Q3: the sequence backbone reconstructs Kozak on its own (r=0.807), so an explicit
start-context channel is redundant. **Recommendation: drop the heuristic -- build the ORF track with
`--kozak none` (the V1 track).** It removes a hand-tuned literature prior that mildly mis-serves the
non-canonical targets, simplifies the track, and directionally improves held-out detection. Neither the
empirical PWM nor the learnable gate earns its added complexity.

**Implemented (2026-07-16, default-flip only, no deployment retrain).** `build_orf_track.py` now defaults
`--kozak none` (was `heuristic`); the `orf_track()` function default and the module docstring were updated
to match. The `heuristic` and `pwm` paths remain available (opt-in) for reproducibility. IMPORTANT
consistency note: the existing deployment `orf_v2_attn` model was TRAINED on the heuristic track
`data/packed/orf_track_v2.npy`, which is left untouched -- so the deployed model keeps using its heuristic
track (no train/inference mismatch). Only NEW `ext`-mode track builds are no-Kozak by default. A deployment
rebuild (retrain orf_v2_attn on a no-Kozak train-on-8 track) is deferred; recommend a 3-seed confirmation
of the V1 edge first (the cross-tissue deltas are single-seed).

**Caveats.** Single seed per arm on a soft-truth task (RiboCode CTG calls); the cross-tissue deltas
(CTG F1 spread 0.049) are modest and would be worth a 3-seed confirmation before hardening the V1
recommendation into a deployment rebuild. One train/test tissue pair (Fib->Hep). The in-distribution
val-Pearson tie is the robust part; the V1 test-set edge is the suggestive part. Jobs 35403140-35403146
(v0/v2/v3), 35411708-35411710 (v1); orchestration `run_kozak_eval_chain.sh`; driver produced the final
`kozak_summary.md`.

## Task 21: Posture-B multimap sensitivity check (2026-07-17, LANDED)

The pre-registered definitive version of the isoform-multimap check (methods.md 5O). Posture A counts
each footprint on every within-gene isoform (phantom-signal risk on low-expressed siblings); posture B
commits each footprint to one representative isoform. Key equivalence: for a gene's highest-expressed
isoform, posture-A counts already equal posture-B counts, so posture B = "restrict train+eval to the
highest-expressed isoform per gene" with NO target re-derivation (`make_representative_tx.py` ->
`data/fibroblast_representative_tx.txt`, 12,485 genes = 34% of packed tx; `dataset.representative_tx()`
via env `RIBO_REPRESENTATIVE_TX`). Airtight 2x2: {baseline, orf_v2_attn} x {A full, B representative},
one-hot, fold-0, identical code/seed, only the universe differs (`train_postureB.sbatch`, job 35498336).
Compare pc/lncRNA/uORF/dORF pooled medians A vs B. Expected to CONFIRM (Task 6 single-isoform-gene proxy).

| config | pc | lncRNA | uORF | dORF | count Pearson | n test pc |
|--------|----|--------|------|------|---------------|-----------|
| baseline    A (full)           | 0.638 | 0.364 | 0.622 | 0.329 | 0.806 | 6,698 |
| baseline    B (representative) | 0.623 | 0.346 | 0.618 | 0.259 | **0.823** | 2,202 |
| orf_v2_attn A (full)           | 0.647 | 0.394 | 0.626 | 0.329 | 0.808 | 6,698 |
| orf_v2_attn B (representative) | 0.635 | 0.371 | **0.655** | 0.309 | **0.823** | 2,202 |

**Result: the multimap posture does NOT distort the conclusions -- CONFIRMED as pre-registered.** Posture B
(one representative isoform per gene, the clean posture-B target) lands within ~0.01-0.02 of posture A on
pc/lncRNA (B-A: baseline pc -0.015 / lncRNA -0.018; orf_v2_attn pc -0.012 / lncRNA -0.023), and the two
headline conclusions survive intact: (1) the architecture edge holds on B -- orf_v2_attn beats baseline on
every class (pc 0.635 vs 0.623, lncRNA 0.371 vs 0.346, uORF 0.655 vs 0.618); (2) magnitude is unaffected.
Two things even IMPROVE on the clean representative universe: count Pearson rises to 0.823 (vs 0.806-0.808;
representatives have cleaner expression signal) and orf_v2_attn's uORF is the best of all four (0.655). So
the ~18x isoform-multimap inflation is genuinely absorbed by per-transcript normalization, and phantom
signal on low-expressed siblings is not materially inflating posture A. The small B<A on pc/lncRNA is
consistent with either mild phantom-easing in A or representatives (highest-expressed, often longer
canonical isoforms) being intrinsically harder. Caveat: A and B are scored on DIFFERENT test sets (full
6,698 vs representative 2,202), so this mixes training- and eval-composition; a common-subset cross-eval
would isolate them, but the practical multimap question -- does posture A distort the story -- is answered
NO. The pre-registered heavy check is closed; posture A stands. Jobs 35498336 (all 4 arms);
`train_postureB.sbatch`, `make_representative_tx.py`, `dataset.representative_tx()`.

## Task 22: replicate-concordance ceiling -- how much of the profile is even predictable? (2026-07-17, LANDED)

Motivation: architecture changes (inputs, FM vs one-hot, attn depth, mamba, width, Kozak) have all
plateaued at pc profile Pearson ~0.61-0.65. Is that the noise floor (near the ceiling, architecture not
the bottleneck) or is there reproducible headroom the tested models miss? To find out, split the 32
Fibroblast Ribo-seq samples into two independent 16-sample half-pools, compute the SAME per-transcript
profile Pearson the model reports (whole-tx + 5'UTR/3'UTR windows, by biotype) between the halves, over 5
random splits, then Spearman-Brown correct each half-depth r to full-pool depth (`replicate_concordance.py`,
job 35535899; `results/replicate_concordance.json`).

| class | replicate ceiling (SB, full) | raw r (half-depth) | model | model / ceiling |
|-------|:---:|:---:|:---:|:---:|
| pc whole-tx | **0.952** | 0.909 | 0.638 | 67% |
| lncRNA whole-tx | **0.890** | 0.802 | 0.364 | 41% |
| uORF (5'UTR) | **0.953** | 0.910 | 0.622 | 65% |
| dORF (3'UTR) | **0.749** | 0.599 | 0.329 | 44% |

(n = 33K pc / 1.2K lncRNA / 22K uORF-window / 11K dORF-window scorable tx per split; SB = 2r/(1+r).)

**Result: the profile is HIGHLY reproducible, and the model captures only ~2/3 of it (pc) to ~2/5
(lncRNA) -- so profile Pearson is NOT saturated; there is large, real reproducible headroom.** The
architecture-sweep plateau is therefore NOT a noise-floor ceiling. The model nails the frame / 3-nt
periodicity (most of its 0.64), but misses the reproducible position-specific magnitude modulation
(elongation / pausing -- which positions pile up P-sites), which is reproducible at ~0.95 and no tested
architecture (one-hot/FM/attn/mamba/width) captures.

**Reconciliation -- two metrics, two stories (this is the key nuance):**
- **ORF localization / calling** (the project's actual deliverable): the model IS at/above the ceiling --
  Ruiz-Orera pred_frame0 AUROC 0.945 beats the observed-profile ceiling 0.908 (Task 16), plus the depth
  crossover (Task 18). The frame signal is what calling needs, and the model has it near-perfectly.
- **Per-nt profile SHAPE** (elongation dynamics): 41-67% of the ceiling -- a stated LIMITATION.

**Honest limitation (paper framing).** This model is optimized for translation LOCALIZATION, which it
achieves at the replicate-reproducibility ceiling. It does NOT capture the finer position-specific
elongation dynamics (the reproducible profile-shape residual), which is the explicit target of dedicated
ribosome-density models (seq2ribo, Riboformer, RiboMIMO, RiboNN's density head).

**Per-codon (elongation-only) resolution + seq2ribo comparison (2026-07-17).** To separate periodicity
from elongation, the concordance was re-run on the CDS collapsed to in-frame codons (periodicity removed):

| per-codon CDS Shape r (elongation only) | value | note |
|-----------------------------------------|:-----:|------|
| replicate ceiling (SB, full-depth) | **0.950** | half-depth 0.905; n=31.6K pc |
| this model (orf_v2_attn fold-0) | **0.526** | ~55% of ceiling; per-nt sanity 0.614 |
| seq2ribo (Kaynar & Kingsford 2026, published "Shape r") | **0.05-0.19** | its own within-transcript metric |

Three findings that overturn the naive read: (1) **the headroom is NOT periodicity** -- the per-codon
elongation profile is *itself* reproducible at ~0.95, so the codon-to-codon pausing pattern is real
biology, not noise; the per-nt 0.95 was not an artifact of the deterministic 3-nt frame. (2) **This model
captures ~55% of the reproducible elongation signal (0.53/0.95)** -- moderate, with real remaining headroom.
(3) **It substantially OUTPERFORMS the purpose-built seq2ribo on within-transcript shape** (0.53 vs its
published Shape r 0.05-0.19), because this model is trained on the within-transcript profile (multinomial
NLL) while seq2ribo optimizes cross-transcript magnitude -- seq2ribo's advertised 0.92 is "Elemwise r"
(pooled across transcripts ~ TE/expression), NOT within-transcript shape. So seq2ribo is NOT the tool that
closes this gap; if anything this model is closer. Caveat: seq2ribo's 0.05-0.19 is on its GWIPS-viz A-site
data (domain + A-vs-P-site + median-vs-mean differ), so the cross-study number is directional; but since its
within-transcript shape is weak even on its home turf, the direction is robust.

**Corrected honest limitation (paper framing).** This model achieves translation LOCALIZATION at the
replicate ceiling (the deliverable) AND captures roughly half the reproducible codon-level ELONGATION
profile (0.53/0.95) -- better than the dedicated SOTA sequence model -- while ~half of that reproducible
signal remains unpredicted by ANY current sequence model. That residual (codon-context-specific ribosome
pausing) is a real, shared, open problem, not a tuning gap this architecture family can close. Decision:
cite seq2ribo's published Shape r rather than stand it up domain-matched (its home-turf number already
makes the point; an out-of-domain run on Fibroblast would only be <= that). `replicate_concordance.py`
(pc_cds_codon), jobs 35535899 + 35539171; per-codon model number from
`results/improve/orf_v2_attn_f0/dropin/pred_profiles.npz`.

## Task 25: mm1 vs mm20 RNA-seq coverage posture -- fast proxy (2026-07-21, LANDED)

Question the user raised ("run train + predict on uniquely-mapped RNA and compare"): the deployed RNA-seq
COVERAGE input is built with STAR `--outFilterMultimapNmax 20` (mm20). Does single-mapper (mm1) coverage
change the model? Retraining on mm1 would need re-downloading + re-aligning the 57 Chothani FASTQs (the
training FASTQ were deleted) -- expensive. Fast proxy instead: align one 30M-read DoHH2 RNA-seq subsample
BOTH ways (mm1 and mm20), build per-nt transcript coverage each way (`align_proxy_mm.sbatch`,
`--outFilterMultimapNmax {1,20}`), and compare per-transcript on the DoHH2 expressed universe
(`compare_coverage_mm.py`).

| metric (DoHH2 universe, 30,950 tx with mm20 depth >= 50) | value |
|---|---|
| median per-tx Pearson(mm1, mm20) | **1.0000** |
| tx with Pearson >= 0.999 | **85.7%** |
| tx with Pearson < 0.95 | 5.3% |
| total depth ratio mm1/mm20 | 0.9236 (mm20 has ~8% more reads, all multimappers) |
| multimap-inflated tx (mm20 > 1.5x mm1) | 983 (3.2%) |

**Result: the multimap posture does NOT materially change the model input -- a full mm1 retrain is
unwarranted.** Coverage is identical (Pearson >= 0.999) on 85.7% of expressed transcripts and only 5.3%
differ at all; the divergent tail is exactly the paralog / repeat / multimap-heavy loci (the 10 lowest-
Pearson transcripts have mm20/mm1 depth ratios 7-145x -- histone clusters, paralog families, repeat-
embedded transcripts). Because the model is per-transcript-normalized and sees identical input on ~95% of
transcripts, its outputs are insensitive to the posture. This is consistent with Task 21 (posture-B: the
~18x isoform-multimap inflation is absorbed by per-transcript normalization; conclusions posture-invariant)
and confirms the standing note that `--outFilterMultimapNmax 20` is the correct RNA-seq-COVERAGE posture
(the single-mapper rule is Ribo-seq-only, for P-site periodicity). The Chothani mm1 re-download was
therefore stopped (ENA was throttling it to ~2 days anyway). Figure: `figures/B7_multimap_posture/`.
Proxy aligns: jobs on `align_proxy_mm.sbatch`; comparison `compare_coverage_mm.py` ->
`proteogenomics/data/mm_proxy/coverage_mm_{comparison.txt,pertx.tsv}`.

## Figure set for the brief communication (2026-07-21, building)

Reusable per-figure generators under `figures/<name>/` (each: `make_<name>.py` + `FIGURE_DATA_INPUTS.md`
with a `## Regenerate` command + PDF/PNG). Master index: `figures/README.md`. Built this session (all read
result JSON/TSV so they regenerate on retrain): A1 localization-beats-ceiling, A2 RiboCode drop-in F1 0.923,
A3 even 9-fold LOTO generalization (profile-Pearson spread is held-out target quality, r=0.81 vs period_obs;
count transfer even), B4 replicate ceiling, B5 translation signal exceeds an expression-only baseline
(+0.12 length-ctrl), B7 mm1/mm20 posture (Task 25), C10 input saliency (start codon among the most salient
nt in the transcript; the model reads the start context), D11 vs seq2ribo. Prior PNG figures adopted:
prediction_examples, kozak (C9), depth_crossover. GATED on the immunopeptidome MS consolidation (Fig 2):
B8 class-vs-global FDR, D12 CPAT/CPC2 baseline, D13 discovery forest plot, D14 novel-antigen case study.

### mm1 adopted as the standard + strict-uniform re-run (2026-07-21, user directive)

The user chose STRICT UNIFORM mm1: since the proxy proves mm1 == mm20 on the coverage input, unique mappers
become the standard everywhere (smaller intermediates, faster STAR, one posture shared with Ribo-seq).
Actions: (1) switched --outFilterMultimapNmax 20->1 in all 8 RNA-seq COVERAGE scripts + methods.md sec 2;
(2) re-running the completed Fig 2 immunopeptidome datasets end-to-end on mm1 coverage so nothing in the
figure is mixed-posture. Fully SLURM-native afterok chain per line (no head-node babysitting):
align(mm1, SKIP_SALMON=1 since the universe is salmon-derived and unchanged) -> rebuild the model-input
pack from the new coverage -> GPU predict -> enumerate/score ORFs -> model+null DBs -> HLA search ->
MS2Rescore -> class-FDR compare (-> db_comparison_rescored_mm1.md). Driver: rerun_mm1_line.sh <LINE>
<HELDOUT> <SRR>; new SLURM steps post_predict_line.sbatch + compare_line.sbatch. Launched DoHH2
(35864420..23), SUDHL4 (35864424..27), HBL1 (35864428..31). The prior SU-DHL-4 mm20 predict (35864407)
was cancelled. A549 (tryptic; RNA FASTQ had been deleted) is re-downloading its 4 ENCODE FASTQs to re-run
mm1 too. The deployed model itself is NOT retrained (mm20-trained but proxy-equivalent on mm1 input; the
Chothani training FASTQ are gone). Expectation per Task 25: the mm1 results match the mm20 results within
noise; the value is uniform provenance, not a changed answer.

A549 (tryptic) mm1 re-run also launched (35864436-39) via rerun_mm1_a549.sh: its canonical "fresh universe"
lineage reads the ENCODE FASTQ (ENCSR000CON, re-downloaded, gzip 4/4 OK) projected onto A549's own universe;
chain align_a549_encode(mm1) -> pack human_a549_fresh -> dump_a549_fresh -> enumerate -> db_fresh -> tryptic
msf_full -> ms2rescore_fanout -> class-aware compare (-> _mm1.md, mm20 preserved as _mm20.md). All 4 Fig 2
datasets (DoHH2, SU-DHL-4, HBL-1, A549) now re-running strictly on mm1 coverage.

### mm1 strict-uniform re-run: outcome + a stale-rescore bug + a mokapot-stochasticity finding (2026-07-21)

The 3 immunopeptidomes (HBL-1, DoHH2, SU-DHL-4) were re-run end-to-end on mm1. Two things surfaced:

1. **Stale-rescore bug (fixed).** ms2rescore_line_fanout.sbatch has a resume guard (skip if psms.tsv
   exists). The old mm20 rescore PSMs were not cleared, so the first mm1 pass SKIPPED rescoring (1 s/task)
   and the compare read stale mm20 PSMs -> falsely byte-identical mm1==mm20. Caught by the suspicious exact
   identity. Fixed: post_predict_line.sbatch + post_predict_a549.sbatch now clear stale rescore/search PSMs
   before submitting; re-ran rescore+compare clean.

2. **The DB is posture-invariant; peptide counts are mokapot noise.** The deterministic, posture-affected
   artifact is the model's novel-ORF DB, and its counts are IDENTICAL mm1 vs mm20 to the digit (model
   72,439 / 60,814 / 67,670; null 220,190 / 190,619 / 205,041). The downstream peptide discovery counts
   wobble (e.g. DoHH2 model 17 vs 18 on two runs of the SAME mm1 search) -- MS2Rescore/mokapot is UNSEEDED
   and stochastic on the sparse novel class. So mm1 and mm20 give the same DB and the same discovery within
   rescoring noise -> the posture is safe (confirms Task 25 end-to-end).

**Process caveat:** compare_rescored.py hardcodes db_comparison_rescored.md and compare_line.sbatch tees the
same run to _mm1.md, and rescore/*.psms.tsv is a shared dir -- so the clean mm20 peptide baseline was
overwritten during the re-run. It was not needed for the conclusion (DB-count identity + the determinism
test carry it), but the compare/rescore I/O should be made posture-namespaced before any future A/B.

**NEW Fig 2 rigor item (independent of posture): the Fig 2b discovery counts have mokapot run-to-run
variance (+/- ~1-2 on these sparse counts).** Before finalizing Fig 2b, either seed mokapot for
reproducibility OR report the counts as a mean +/- SD over N rescoring runs (error bars on the forest plot).

### Strict-uniform mm1 re-run COMPLETE (2026-07-21/22) -- posture-invariance confirmed on all 4 Fig 2 datasets

All 4 immunopeptidome/proteome datasets re-run end-to-end on mm1 coverage. Results match the preserved mm20
baselines within noise, confirming the multimap posture does not change the Fig 2 conclusions:
- HBL-1 / DoHH2 / SU-DHL-4 (immunopeptidome): model novel-ORF DB counts IDENTICAL mm1 vs mm20 (deterministic);
  peptide discovery matches within mokapot stochasticity (characterized separately, D13 / stabilized).
- A549 (tryptic): model DB 98,157 (mm1) vs 98,047 (mm20) = +0.1%; ALL discovery numbers identical -- raw
  hyperscore class-FDR model 10 / null 13 (2.31x by rate), global-rescore 1/6, class-aware UNTRAINABLE (tryptic
  novel class too sparse). Fix applied: compare_a549.sbatch now uses the ms2rescore env (compare_rescored_
  classaware.py imports mokapot, absent from cas12a). mm20 preserved as *_classaware_mm20.md.

Net: strict-uniform mm1 achieved; the deployed model is NOT retrained (proxy + full re-run both show
equivalence); mm1 is now the standard coverage posture (smaller intermediates, faster STAR, consistent with
Ribo-seq). Queue clear.

## De novo ORF calling: no-Kozak mm1 retrain, count-head calibration, and the reproducibility ceiling (2026-07-24)

Detailed running plan + all numbers: `docs/count_head_calibration_checks.md`. Headlines below.

### Final model (mm1 + no-Kozak + Brain-drop) is localization-equivalent to the deployed model
Retrained the one-hot `orf_v2_attn` holding out Hepatocytes, on 7 tissues (Brain dropped as noise), mm1
coverage, no-Kozak ORF track (`orf_track_v2_nokozak.npy`; Kozak heuristic never set per Task 20). Run
`results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes`, best val_pearson 0.527.
Localization on the human Hepatocytes LOTO holdout: `pred_frame0` AUROC **0.941** vs deployed 0.940 (both
beat the observed ceiling 0.907). The three changes (mm1, no-Kozak, drop-Brain) are localization-neutral --
no regression, no obvious gain on that metric. All three tie within ~0.004 (seed noise).

### But localization can't see the real problem: the standalone model OVER-CALLS non-canonical ORFs
Standing evaluation now runs the RiboCode drop-in (`real` / `pred_obsdepth` = shape at real depth /
`pred_preddepth` = standalone) on predicted profiles, per dataset, both arms (see below). On Hepatocytes the
shape is fine (`pred_obsdepth` F1 ~0.80) but `pred_preddepth` over-calls: precision 0.71 (no-Kozak) / 0.74
(deployed). Over-call is concentrated in short non-canonical ORFs and is DEPTH-dependent (worst on shallow
data). Root cause is NOT count-scale: it is SMOOTHNESS -- the model predicts a noiseless *rate*, RiboCode's
frame test is tuned for noisy integer counts, so spurious weak ORFs look cleanly periodic. Proof: `pred_obsdepth`
(perfect per-tx counts) still has novel precision only 0.31; per-tx isotonic calibration cannot help.

### Count-head calibration: Poisson injection wins; the CDS-anchored stringency dial
Added `--pred_scale theta` (scale predicted effective depth) + `--pred_poisson` (sample Poisson(rate) instead
of rounding) to `ribocode_dropin.py`. Check 4, matched CDS recall, Hepatocytes: **Poisson injection is the
best lever** (+0.10 novel precision vs plain depth-scale at 90% CDS recall; p-value tightening WORSE -- the
smooth density gives spurious ORFs deceptively good p-values). Recommended de novo recipe:
`Poisson-sample(predicted_rate x theta) -> RiboCode`, theta = a **CDS-anchored** operating point (CDS is the
one class we trust; CDS precision ~0.93 across the whole sweep). The **dial** (human Hepatocytes, 3-seed
stable): theta=0.05 -> ~92% CDS recall, **226 novel @0.58, 660 uORF @0.73**; tighten for fewer/higher-precision
novels. STANDING RULE: every ORF-calling test reports BOTH arms (standard pred_preddepth theta=1 + downsample+
Poisson at CDS-anchored theta). Sweep TSVs preserved as supplemental-figure data.

### O2: the between-experiment reproducibility ceiling (the honest de novo bar)
Two independent mouse-liver Ribo-seq datasets (Wang, Janich; both pass periodicity QC -- 76-93% / 85-92%
frame-0) agree only at **F1 0.675 ALL / 0.735 CDS / 0.505 non-canonical** (`o2_between_dataset.py`, coord-key
matched). So non-canonical ORFs are only ~50% reproducible even between real experiments -- the model must be
judged against ~0.5, NOT 1.0. A predictor reaching ~0.5 non-canonical F1 is "as informative as a second real
experiment."

**Model vs ceiling (Wang no-Kozak, 2026-07-24) -- CDS reaches the ceiling once calibrated.** Standard arm
(theta=1) de novo: CDS F1 0.621 (0.84 of ceiling), non-canonical 0.067. The CDS shortfall is pure over-calling
(the model makes 21,315 CDS calls vs ~14,700 real -> precision 0.51 despite recall 0.786 > ceiling's 0.756).
The **Poisson arm fixes it**: at theta=0.02 the calibrated model makes 13,251 CDS calls, precision 0.76, and
**CDS F1 = 0.744, edging PAST the 0.735 between-experiment ceiling** -- de novo (sequence + RNA-seq, no Ribo)
predicts liver CDS ORF calls as well as a second real Ribo-seq experiment, matching the with-Ribo arm (0.759).
Non-canonical improves ~3x under Poisson (0.067 -> 0.195 at theta=0.05) but stays far below its 0.505 ceiling
-- the FP-filter + peaky model target that gap, not the depth dial. `results/o2_liver/o2_{,poisson_}summary.txt`.
Janich symmetric model-vs-ceiling still queued (dump 35890103).

### Transfer + the class-specific verdict (Check 5, GSE120762, deployed model)
The dial KNOB transfers: theta=0.05 -> ~90% CDS recall on human liver, mouse BMDM, shallow AND deep. But
non-canonical precision at an operating point is NOT a portable number -- it depends on the OBSERVED reference
depth. Deep-reference follow-up (GSE-NT calls vs deep GSE-LPS observed) resolves it CLASS-SPECIFICALLY: CDS
solid (0.83-0.89); uORF 0.23->0.46 (shallow reference unfairly penalized it -- half the "misses" are real);
novel 0.07->0.11 vs a 0.54 baseline (model novels validate at only ~20% of the real-novel rate) = **novel is
GENUINELY over-called, not an artifact**. Net thesis: CDS + (corrected) uORF calls are defensible; novel-ORF
calling needs a robust caller, and must be scored vs a DEEP independent reference (never a shallow same-dataset
one). This motivated O3.

### O3 (forward goal): robust novel-ORF caller beyond RiboCode's frame test
Move past the single univariate periodicity test (what the smooth density fools). Two tracks: co-opt an ML
caller runnable on a substituted density (RibORF SVM / RP-BP / DeepRibo), or a lightweight per-ORF FP-filter
classifier trained on reproducibility labels. Survey verdict (`docs/o3_orf_caller_miniplan.md`): NO wholesale
co-opt -- every periodicity-based caller (RiboCode/Ribo-TISH/Ribotricer/RP-BP/ORFquant) inherits the same bug
on a smooth predicted density, and the discriminating signal is read-distribution SHAPE, not periodicity. So
keep RiboCode for candidate generation + build a custom FP-filter.

**O3 track 1 -- FP-filter (`scripts/heldout/{orf_features,fp_filter}.py`, 2026-07-24):**
- **Phase 1 refuted the central premise.** Extracting shape features (codon uniformity PME, Gini/CV, 5' ramp,
  3'-of-stop drop, in-frame fraction, length) from the model's predicted density and labelling by the deep
  observed calls: PME alone gives AUC ~0.50 for novel. The model over-smooths REAL ORFs too, so the shape
  info that separates real from spurious is largely absent from its output. Best multivariate 5-fold CV =
  logreg = RandomForest = **AUC 0.669** (novel) -- a linear, modest signal.
- **Trained FP-filter (per-class balanced logreg, CV):** novel base rate 0.235 -> AUC 0.668; lift 1.5x at
  50% recall. It leaves CDS untouched (annotated AUC 0.591, precision ~0.95 held).
- **On the CDS-anchored framework (in-sample):** at fixed CDS recall 0.956 (theta=1.0 standard arm), the
  filter lifts novel precision 0.294 -> 0.398 (+0.10) and uORF 0.425 -> 0.540 (+0.12) -- its niche is raising
  non-canonical precision WITHOUT the CDS-recall cost that lowering theta pays (theta=0.05 reaches 0.590 novel
  but drops CDS recall to 0.925). Complementary levers. Caveat: Hepatocytes-trained + applied -> optimistic
  upper bound; cross-dataset (Wang/Janich) validation pending GPU.
- **Verdict:** a real but MODEST post-hoc lever, ceiling-capped by the model's smoothness.

**O3 track 2 -- anti-smoothing model (built + training in parallel, 2026-07-24):** the higher-ceiling fix.
`model.py::profile_entropy_gap` adds `--peakiness_weight * relu(H(pred) - H(obs))` to the profile loss so the
model stops hedging with smooth profiles (the multinomial NLL alone tolerates a smooth prediction). Deep pooled
targets (~6.6 P-sites/nt) mean H(obs) reflects real 3-nt periodicity, not shot noise. `train_loto_noBrain_peaky.sbatch`
(PK=0.2) is an exact A/B vs the current no-Kozak mm1 model (same config/holdout + the term), queued 35890525
behind the O2 dumps. A/B when it lands: does feature-separation AUC clear 0.669, does novel over-calling drop
at matched CDS recall, does val_pearson hold. See Check 8 + the mini-plan.

## Task 26: 3-seed union comparison -- attn vs mamba4, settled (2026-07-31, LANDED)

Which architecture ships as the primary model. Both are trained on the union universe with the
no-Kozak ORF track (`--kozak none`) and mm1 coverage; only the sequence mixer differs. Three seeds
each, scored on held-out test median per-transcript Pearson.

| seed | `orf_v2_attn` (dilated CNN + 2 transformer layers) | `orf_v2_mamba4` (dilated CNN + 4 Bi-Mamba blocks) |
|---|--:|--:|
| 0 | 0.6585 | 0.6851 |
| 1 | 0.6603 | 0.6753 |
| 2 | 0.6598 | 0.6794 |
| **mean** | **0.6595** | **0.6799** |
| spread | 0.6585 - 0.6603 | 0.6753 - 0.6851 |

**Gap = +0.0204 in favour of mamba4, and the seed ranges do not overlap: mamba4's worst seed (0.6753)
beats attn's best (0.6603).** At n=3 per architecture this is a clean separation, not a seed artefact.
Note the two architectures differ markedly in seed stability -- attn's spread is 0.0018 while mamba4's
is 0.0098, i.e. mamba4 is ~5x noisier across seeds. It wins anyway, but a single-seed mamba4 number
should not be quoted without the spread.

Params: attn 5,071,106; mamba4 7,521,026. Mamba4 is GPU-only (mamba-ssm CUDA kernels).

**Both models ship** (user decision, 2026-07-31): mamba4 is the primary model and the headline result;
attn moves to supplemental for the paper but stays maintained as the released inference path for
CPU-bound users, since it is the only one of the two that runs without a GPU. Both are covered by the
architecture figure pair (`figures/arch_attn/arch_rinalmo.png`, `arch_rinalmo_mamba.png`) and the
tutorial's architecture page.

Provenance note: attn seed1 initially hit the 24 h wall at epoch 23 and was relaunched with a 48 h
limit and `--resume`. The resumed run was explicitly verified to restart from epoch 23 with its prior
best (`val_pearson=0.6015`) rather than silently reinitialising from epoch 0, so seed1's 0.6603 is a
genuine continuation, not a short run.

## Task 46: macrophage proteogenomics -- model-selected vs null ORF DB, 12 populations (2026-07-30, LANDED)

Does restricting the MS search database to model-selected ORFs beat a size-matched null selection?
12 mouse macrophage populations, 214 mzML. Novel = peptide maps ONLY to model/null-called ORFs and is
scored against `REV_nuORF|` decoys only (class-specific FDR at 1%; a global FDR inflates non-canonical
discovery ~10-13x and must not be used here). Full per-population table:
`proteogenomics/data/macrophage_tissue/macro_model_vs_null.md`.

| metric | canonical baseline | model DB | null DB |
|---|--:|--:|--:|
| novel peptides (1% class FDR) | -- | **440** | 374 |
| canonical PSMs | 9,961,958 | 9,430,641 | 8,481,926 |
| **PC-churn** (vs baseline) | -- | **-531,317** | **-1,480,032** |
| net PSM (canonical + novel targets) | -- | 9,975,276 | 10,074,790 |

> **RETRACTED 2026-08-01 -- the churn and net-PSM columns above are NOT FDR-filtered.** See Task 52.
> `macro_churn_aggregate.analyze()` applied the FDR loop only to the novel-peptide column; `canon` and
> `total_t` counted every rank-1 PSM. Unfiltered, those columns measure how many spectra a database
> ABSORBED, which scales with database size. Corrected findings:
>
> - **The -1,480,032 PC-churn does not survive FDR.** Canonical PSMs at 1% FDR are 335,548 /
>   335,804 / 335,392 / 336,284 for canonical-only / model@0.5 / model@0.74 / null (BMDM) -- a
>   **0.27% spread**. Across 12 populations the shift is -0.6% to +1.0%. The displaced PSMs were
>   sub-threshold matches, never confident identifications. The unfiltered number overstated the
>   effect 10-25x.
> - **The original sentence "adding poorly-chosen ORFs actively costs previously-confident canonical
>   identifications" is WRONG** and is withdrawn. Those identifications were not confident.
> - **The null's net-PSM "win" was an absorption artefact**: 1,592,864 unfiltered novel rank-1 PSMs
>   yielding only 374 FDR-surviving peptides (4,259 raw matches per real peptide, vs 399 for the
>   CDS-anchored arm). Its novel gain was ~93% accounted for by its canonical loss -- it was
>   relabelling canonical matches as junk-novel.
>
> **What survives unchanged** is the novel-peptide axis, where every arm uses the same decoy class,
> procedure and denominator: model 440 vs null 374, and the model wins the per-population discovery
> rate in all 12 (2.94x to 8.05x). That comparison is untouched by the artefact. The claim from this
> experiment must be a discovery-rate claim, not a churn claim and not a net-PSM claim.

Caveat on provenance: the macrophage MODEL database was built with
`orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes` -- the pre-union deployed checkpoint, a
7-tissue LOTO model (train `Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC`, hold out Hepatocytes) applied
cross-species to mouse. (An earlier draft of this note called it "the deployed fibroblast checkpoint";
that was wrong. Fibroblast is one of seven training tissues and the name of the default universe FASTA,
not the model's training set.) It is NOT the union-universe attn/mamba4 pair selected in Task 26, so
these numbers are a lower bound on the shipped architectures; the rebuild is Task 50.

## Task 52: the DB-tradeoff metric was measuring database size (2026-08-01, LANDED)

Triggered by a question that should not have had an answer: how can the CDS-anchored macrophage DB
report MORE novel peptides and LESS canonical churn than the null, yet a LOWER net PSM?

**Cause.** `macro_churn_aggregate.analyze()` FDR-filtered only the novel-peptide column:

```python
canon   = sum(1 for (_,_,k) in sel if k == "canon_t")             # NO FDR
total_t = sum(1 for (_,_,k) in sel if k in ("canon_t","novel_t")) # NO FDR
nov     = ...                                                     # <- FDR applied only here
```

Unfiltered, `canon` and `total_t` count how many spectra a database ABSORBED. A larger database
absorbs more by chance, so both columns scale with database size rather than with correctness.

**Decomposition (12 populations).** Splitting net PSM into its parts explains the paradox exactly:

| arm | canon PSM | novel_t PSM | net PSM | novel pep @1% | raw PSM per real peptide |
|---|--:|--:|--:|--:|--:|
| deployed attn @0.5 | 9,430,641 | 544,635 | 9,975,276 | 440 | 1,238 |
| mamba4 @0.5 | 9,307,176 | 687,819 | 9,994,995 | 446 | 1,542 |
| **mamba4 @0.74 anchored** | 9,804,842 | **184,016** | 9,988,858 | **461** | **399** |
| null (2.36M seqs) | 8,481,926 | **1,592,864** | 10,074,790 | 374 | **4,259** |

The null's entire net-PSM lead is 1,592,864 junk matches that do not survive FDR. Its novel gain
(+1,592,864) is ~93% accounted for by its canonical loss (-1,480,032): it is mostly RELABELLING
canonical matches as junk-novel, then being credited for them by a column that cannot tell.

**The churn claim does not survive FDR either.** Canonical PSMs retained at 1% FDR, BMDM:

| db | canon_t raw | other_d raw | FDR thresh | **canon_t kept @1%** |
|---|--:|--:|--:|--:|
| canonical-only | 849,994 | 361,892 | 18.517 | **335,548** |
| model @0.5 | 811,734 | 325,175 | 18.509 | **335,804** |
| model @0.74 | 836,911 | 349,132 | 18.522 | **335,392** |
| null | 737,895 | 255,138 | 18.484 | **336,284** |

A **0.27% spread** (-0.6% to +1.0% across all 12 populations) against an unfiltered churn of -5.3% to
-14.9%. The displaced PSMs were sub-threshold matches, never confident IDs. **Task 46's statement that
poorly-chosen ORFs "actively cost previously-confident canonical identifications" is withdrawn.**

**A second, subtler trap.** Naively fixing this by applying class-specific FDR to the canonical side
produces an impossible result -- the null gaining +3,317 canonical peptides over a canonical-only
search. Junk cannot create real canonical IDs. Mechanism, visible in the table above: a large novel
space soaks up SPURIOUS matches preferentially, so canonical decoys are cannibalised faster than
final-recipe targets (`other_d` -29.5% vs `canon_t` -13.2%), deflating the estimated canonical FDR and
relaxing its threshold. So the unfiltered and the naive class-specific-FDR versions are BOTH biased
toward the larger database, by different mechanisms.

**Resolution (now a standing rule).** Every column in a DB-comparison table must be FDR-filtered by
the same procedure; never place a filtered column beside an unfiltered one. Canonical/net columns use
GLOBAL 1% FDR (keeps the target-decoy competition intact); class-specific FDR is for the novel column
only, where the denominator is the point. Report discovery rate per 100k DB sequences, since raw
counts are not comparable across databases of different size. Reference impl for the corrected
computation: `proteogenomics/scripts/net_psm_fdr.py`.

**What was unaffected.** The novel-peptide comparison -- same decoy class, same procedure, same
denominator across every arm -- stands: **461 (CDS-anchored) > 446 (mamba4@0.5) > 440 (deployed@0.5)
> 374 (null)**, at 277.6 vs 15.9 discoveries per 100k sequences (17.5x). The macrophage conclusion is
a discovery-rate result, not a churn result.

## Task 53: `pgx` -- RiboCode-called search databases replace the f0 threshold (2026-08-01, BMDM LANDED)

Retires the `pred_frame0`-threshold database build. ORFs are now called by Poisson-calibrated
RiboCode on the model's predicted signal, N-terminal extensions are tested separately, and both
naive nulls are built from the same expressed universe. Pipeline: `proteogenomics/scripts/pgx/`;
methodology in `methods.md` 5S; plan and full numbers in `docs/proteogenomics_pipeline_miniplan.md`.

### Why the f0 threshold had to go

`f0` is a ratio over the ORF interval, so prepending a zero-signal upstream region changes neither
numerator nor denominator: it is INVARIANT to N-terminal extension. Measured on BMDM over 12,719
stop-codon groups, **96.2%** of 51,646 upstream non-AUG ORFs sat within 0.05 of the true CDS's
`f0`, median |delta f0| = **0.0026**. There was also no Poisson step, no periodicity test and no
significance filter anywhere in that path.

### BMDM headline (mamba4 union, theta* = 0.05, tryptic, 18 fractions)

| arm | novel PSMs | novel pept | novel seqs w/PSM | DB novel seqs | GENCODE PSMs | dPSM | dPept | ncStart |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| gencode (baseline) | 0 | 0 | 0 | 0 | 335,548 | +0 | +0 | 0 |
| **model (RiboCode-called)** | 121 | 40 | 20 | **1,533** | 335,187 | **-361** | -61 | **3** |
| null_atg (naive AUG) | 177 | 50 | 104 | 235,686 | 325,254 | **-10,294** | -2,048 | 0 |

Novel columns at 1% class-specific FDR; GENCODE columns and dPSM at 1% global FDR. The GENCODE
baseline reproduces the legacy final-recipe arm EXACTLY (335,548 both), so the two builds are on
identical footing and only the novel content differs.

- **Database size:** 154x smaller (1,533 vs 235,686 novel sequences).
- **Canonical cost:** 28x lower (-361 vs -10,294 PSMs; 0.11% vs 3.07% of baseline).
- **Discovery density:** 26.1 vs 0.21 novel peptides per 1,000 DB sequences, a **123x** difference.
- **Net trade for using the null instead:** +10 novel peptides for -1,987 canonical peptides.

### The model recovers peptides the 154x larger null cannot report

Peptide overlap: 23 shared, **17 model-only**, 27 null-only. Every model-only peptide is accounted
for mechanistically:

| cause | n | explanation |
|---|--:|---|
| absent from the null's sequence space | 3 | non-AUG start; unreachable by AUG-only enumeration at ANY threshold |
| present in the null DB but below its FDR threshold | 13 | search-space inflation raises the class-specific hyperscore cut from **20.20 to 28.10** (+7.91); all 13 clear the model's cut |
| present but not seen as rank-1 | 1 | lost to rank-1 competition in the larger space |

So a bigger database does not merely cost canonical identifications, it **actively destroys
non-canonical discovery** it structurally contains: 13 real peptides sit in the null's own FASTA
and cannot be reported at 1% FDR because the null's own size raised the bar.

### Extensions (the new, separately-tested entry class)

RiboCode cannot produce these: `orf_finder.orf_find` sets `alt_flag = 0` as soon as an in-frame ATG
precedes the stop, so a near-cognate extension of an annotated CDS is unreachable at any parameter
setting. `pgx.extensions` applies RiboCode's own frame statistic to the **extension region alone**
(candidate start -> annotated CDS start). BMDM: 30,243 candidate starts -> 884 testable ->
**564 pass** q<=0.05 (CTG 152, GTG 95, AAG 72, AGG 67, TTG 46, ACG 37, ATC 34, ATT 27, ATG 18,
ATA 16). Of the 884 that reached the test the retired whole-ORF `f0 >= 0.5` criterion would have
admitted **884 of 884**, and by inheritance essentially all 30,243: a 54x reduction on evidence.
RiboCode itself labelled 236 `annotated` calls as extensions but only 18 ATG extensions survive the
region-specific test, so its most-5'-start choice is unsupported ~92% of the time.

### Comparison with the retired build (same spectra, same FDR treatment)

| build | novel pept | DB novel seqs | dPSM | pept per 1,000 seqs |
|---|--:|--:|--:|--:|
| legacy model (f0 >= 0.5) | 46 | 13,835 | -2,526 | 3.3 |
| legacy null (all ATG candidates) | 41 | 187,043 | -9,114 | 0.2 |
| **pgx model (RiboCode-called)** | 40 | **1,533** | **-361** | **26.1** |

Honest reading: absolute novel yield is slightly LOWER than the retired build (40 vs 46 peptides).
The gains are in efficiency and canonical preservation -- 9x fewer sequences, 7x lower canonical
cost, 8x higher discovery density -- plus a statistical basis the f0 cut never had, and 3 peptides
no AUG-only pipeline can reach.

### Four defects found, three of them in this work

All caught by controlled comparison, not by inspection; none crashed anything.

1. **Absolute CDS-recall anchor unreachable.** Recall of all annotated CDS plateaus at 0.533
   (10,236 of 19,190), so a 0.90 target can never be met. Re-anchored on the achievable CDS set;
   theta* = 0.05 then reproduced the Check 5 operating point independently.
2. **Null missing an ORF class the model had.** RiboCode maps out-of-frame CDS overlaps to
   Overlap_uORF/Overlap_dORF (novel); the null enumerator called them `internal` and dropped them.
   117 of 1,533 model sequences escaped a 2.3M-sequence null, all ATG. After aligning the
   classifiers: 0 escape null_nc, and 546/546 (100%) escaping null_atg do so by start codon.
3. **Enzyme override.** `search_enzyme_nocut_1` was set to `""`, silently switching trypsin to
   trypsin/P. Detected because the GENCODE-only baseline read 352,307 against the legacy arm's
   335,548 (+5.0%) on byte-identical databases. Fixed to `P`; the baseline now matches exactly.
4. **`num_slices` is not an MSFragger 4.2 parameter** (logged as "Unknown parameters"); removed
   rather than left as a warning-generating no-op. MSFragger sizes slices from `-Xmx`.

### Status

`null_nc` (2.5M-target near-cognate null) still searching; its arm will be appended. Open decision:
`--biotype-field` transcript_type (current default, matches all prior universes) vs gene_type
(adds NMD / retained_intron isoforms, where non-canonical ORFs concentrate; A549 45,166 vs 41,139).

### Task 53 addendum: five-arm tables, both shipping models (2026-08-01)

Both calling arms now get their own database and search by default (`pgx.run --db-arms`, default
`standard,poisson`), so the benefit of calibration is visible rather than asserted. BMDM, tryptic,
18 fractions; `gencode`, `null_atg`, `null_nc` are shared between the two models (identical
databases, searched once via `--shared-search-root`).

| arm | novel PSMs | novel pept | seqs w/PSM | DB novel seqs | GENCODE PSMs | dPSM | dPept | ncStart |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| gencode | 0 | 0 | 0 | 0 | 335,548 | +0 | +0 | 0 |
| mamba4 standard (theta=1) | 86 | 31 | 21 | 10,454 | 334,790 | -758 | -138 | 3 |
| **mamba4 poisson (theta*=0.05)** | 121 | **40** | 20 | 1,533 | 335,187 | **-361** | -61 | 3 |
| attn standard (theta=1) | 77 | 30 | 16 | 5,991 | 334,963 | -585 | -122 | 0 |
| **attn poisson (theta*=0.05)** | 88 | 33 | 11 | 1,244 | 335,158 | -390 | -72 | 2 |
| null_atg | 177 | 50 | 104 | 235,686 | 325,254 | -10,294 | -2,048 | 0 |
| null_nc | 657 | 177 | 760 | 2,512,388 | 319,544 | -16,004 | -4,634 | 9 |

**Calibration wins on every axis, in BOTH architectures.** The Poisson arm finds more peptides than
theta=1 (40 vs 31 for mamba4, 33 vs 30 for attn) from 6.8x / 4.8x FEWER sequences and at roughly
half the canonical cost. This was not the expected direction: relaxing calibration was supposed to
trade precision for recall. Instead the extra 8,921 (mamba4) / 4,747 (attn) sequences that theta=1
admits are actively harmful -- they inflate the class-specific FDR threshold enough to bury genuine
peptides the smaller database reports, the same mechanism that costs `null_nc` 16 of the model's own
peptides, at smaller scale.

This closes the question left open by the four-arm table: the model's lower recall against
`null_nc` is **not** an artifact of over-conservative calibration, because loosening theta makes it
worse. Whatever the near-cognate null finds, the caller does not reach it by relaxing stringency.

**mamba4 vs attn: more ORFs, not better ORFs.** mamba4 finds 21% more novel peptides (40 vs 33) at
slightly lower canonical cost (-361 vs -390), so it is the better arm overall. But per-sequence
discovery density is indistinguishable -- **26.1 vs 26.5 peptides per 1,000 DB sequences** -- so the
advantage comes from calling MORE ORFs (1,533 vs 1,244), not from better per-ORF discrimination.
State it that way; this dataset does not show mamba4 is the better discriminator.

## Tasks 60 / 63 / 66: the comparators, the QC panels, and the first ground-truth check (2026-08-07, LANDED)

Detail and full tables live in `proteogenomics/results.md`; this is the model-level summary.

### Task 60 -- CPAT / CPC2 close the "beats a naive enumeration" gap (decision D3)

The honest comparator. CPAT and CPC2 shrink the search space by the same order of magnitude as the
model, from sequence composition alone, and both arms enumerate from the identical candidate pool as
`null_atg` (pinned by `test_coding_potential_pool_matches_null_arm`), so only the selection rule
differs. Novel peptides, 1% class-specific FDR, released mamba4, frozen params:

| dataset | substrate | CPAT | CPC2 | model theta=1 | model Poisson |
|---|---|--:|--:|--:|--:|
| A549 | tryptic | **25** | 23 | 9 | 11 |
| HBL-1 | HLA-I | 3 | 10 | 8 | **14** |
| SU-DHL-4 | HLA-I | 7 | 14 | **22** | 10 |
| DoHH2 | HLA-I | 4 | 8 | **32** | 14 |

**The result splits by assay and the split is the finding.** CPAT/CPC2 win the tryptic whole
proteome; the model wins all three immunopeptidomes. On discovery efficiency (peptides per 1,000 DB
sequences) the Poisson arm leads all four (4.58 / 7.62 / 6.17 / 7.16 vs 0.03-0.11 for the null).
Coding-potential tools score a TRANSCRIPT's composition, so they keep long codon-biased sequences --
what a tryptic digest samples well and the least novel population available. The model scores
per-nucleotide translation from the sample's own RNA-seq, so it keeps short non-canonical
cell-type-specific ORFs, which is what HLA-I presentation samples. Figure `figures/D12_cpat_cpc2/`.

### Task 66 -- B721.221: predicted vs MEASURED ORF calls, the first non-null benchmark

Every earlier comparison scored the model against a null or another selection rule. B721.221 has
both arms in one line: Sarkizova RNA-seq drives the prediction, Ouspenskaia Ribo-seq (327 M unique
footprints) gives the measured calls. The model never sees the Ribo-seq. Genomic keying, restricted
to the model's 11,527-gene space (21,974 measured calls genome-wide -> 16,331 in scope).

| arm | precision | recall | F1 | recall canonical | recall non-canonical |
|---|--:|--:|--:|--:|--:|
| mamba4 theta=1 | 0.780 | 0.580 | 0.665 | 0.754 | 0.241 |
| mamba4 Poisson | 0.887 | 0.506 | 0.644 | 0.714 | 0.099 |
| attn theta=1 | 0.761 | 0.592 | **0.666** | 0.759 | 0.267 |
| attn Poisson | **0.918** | 0.491 | 0.640 | 0.705 | 0.074 |

Three readings. (1) **The two-arm framework behaves as designed on real data**: Poisson buys
precision (0.78 -> 0.89, 0.76 -> 0.92) and pays recall, now confirmed against measured translation
rather than against a null -- the strongest independent support the calibration dial has. (2) mamba4
and attn are indistinguishable AGAIN (F1 0.665 vs 0.666), a third replication after the macrophage
cross-subtype run and the 4-dataset MS panel. (3) **Canonical recall (0.70-0.76) far exceeds
non-canonical recall (0.07-0.27)** -- the class the proteogenomics work depends on is the one the
model recovers worst. That belongs in the manuscript, not in a footnote.

**The split-half ceiling LANDED and it reframes reading (3).** RiboCode run independently on two
disjoint halves of the same Ribo-seq (split by HLA allele so each half holds complete libraries;
A = 18,907 calls, B = 24,538) reproduces its own calls at **F1 0.889**, canonical recall 0.973, and
non-canonical recall **0.498**, scored by the same code on a comparable reference (n_ref 15,703 vs
the model's 16,331). As a fraction of that ceiling the model reaches:

| arm | F1 | recall canon | recall non-canon |
|---|--:|--:|--:|
| attn theta=1 | **75%** | **78%** | **54%** |
| mamba4 theta=1 | 75% | 77% | 48% |
| Poisson arms | 72% | 72-73% | 15-20% |

So two independent measurements of the SAME cells agree on only half of each other's non-canonical
calls. The model's 0.24-0.27 non-canonical recall sits against a denominator of ~0.50, not 1.0: it
recovers about **half the non-canonical signal a replicate Ribo-seq experiment would**, not a quarter
of perfect. The gap is real and stays in the manuscript, but stated against the right denominator.
On canonical ORFs the model reaches 77-78% of the assay's self-agreement using no Ribo-seq at all,
and **75% of the ceiling on overall F1** is the fair headline for a prediction made from RNA-seq
alone against a reference built from 327 M measured footprints.

**FIGURE (added 2026-08-14): `figures/A4_b721_ground_truth`.** This result had no figure for a
week despite being the project's only non-null benchmark. Three panels: predicted-vs-measured
precision/recall/F1, recall split by class against the per-class ceiling, and every model number as a
percentage of the split-half ceiling. The panel-C framing is the one to quote -- panel A's 0.24
non-canonical recall against an implied denominator of 1.0 is the misreading the ceiling exists to
prevent.

### Task 63 -- three supplemental QC panels, all built from bytes already on disk

- **`figures/S_riboseq_qc/`** (closes D5): all 16 Ribo-seq library groups the study uses -- 9
  Chothani training tissues, 5 held-out, 2 proteogenomics cell lines -- on read-length, 3-nt
  periodicity and P-site-offset consistency. **Every library clears the bar**: f0 74.8-87.9% against
  a 33% floor, footprints 28-31 nt, offset 12 nt in most. Two honest observations: HBL-1 is the
  weakest on two axes at once (34 nt mode, 128,641 P-sites at CDS, lowest f0) which is a reason to
  weight its proteogenomics result below the others; B721.221 is by far the deepest (18.8 M), which
  is what makes it usable as the measured-translation reference. Built from RiboCode `metaplots`
  output rather than ribotish, which would have required re-aligning a dozen studies to regenerate
  quantities already on disk.
- **`figures/S_fdr_rigor/`** (replaces the gated B8): **a global 1% FDR fails to control the
  non-canonical class in BOTH directions.** 37.5x over-report on A549 (tryptic, 65,405 canonical
  peptides set the cut at hyperscore 19.5, inside the novel decoy bulk; the class-specific cut is
  30.3) down to a 0.2x UNDER-report on DoHH2. The error tracks canonical-class size, not database
  size -- the four `null_atg` arms span 242k-369k sequences and give 37.5x / 1.9x / 0.8x / 2.2x.
  This corrects the project's earlier "global inflates ~10x" shorthand: it is not merely
  anti-conservative, and in the immunopeptidomes it discards real identifications.
- **`figures/S_rescore_substrate/`** (closes rigor point 4): MS2Rescore helps HLA-I and hurts
  tryptic because mokapot is semi-supervised and needs a learnable target class. The intuitive
  metric points the WRONG way (A549 has 15.1% novel-touching PSMs vs 4.0-4.7% for the
  immunopeptidomes); the metric that matters is CONFIDENT novel targets per 1,000 novel-touching
  PSMs, where A549 sits at 0.16 against 21.9-100.7. Measured, not assumed -- the guess was wrong.

---

## Released-model held-out re-run, both calling arms (2026-08-08, LANDED)

Three things landed together here, because they turned out to be the same problem: the held-out
numbers, the calibration arm, and the scoring convention had each drifted independently.

### 1. Every held-out number is now on the shipping models

Run directories are named for the DATASET, not the CHECKPOINT, so a stale checkpoint is invisible from
the path -- it lives only in the npz `meta` field. Four held-out sets were re-dumped on both released
models and re-scored: Ruiz-Orera (human), Wang / Janich-decontaminated / GSE243134 (mouse liver). The
detailed tables with class-level recall are generated into the tutorial by
`tutorial/make_heldout_{human,mouse}.py`, which prints the source checkpoint beside every row so this
class of staleness cannot recur silently. Headline: shape (`pred_obsdepth`) held to within 0.006;
standalone (`pred_preddepth`) fell 0.052-0.054 on the two sets that had a prior number. See
"Released-model re-run (2026-08-08)" above for the mechanism.

### 2. Both calling arms, on a Ribo-seq-FREE anchor

`scripts/liver_calibrate_freeanchor.sbatch` runs `pgx.calibrate --cds-recall 0.90 --recall-mode
relative` per dump. The anchor is annotated-CDS relative recall -- reachable without any Ribo-seq for
the query sample -- rather than a match to an observed between-experiment ceiling, which would require
the very data the standalone claim says is unnecessary.

`scripts/score_released_two_arm.py`, genomic key, min_len 90, pred enrichment >= 0.5:

| dump | arm | nPred | nRef | P | R | F1 | ncP | ncR | ncF1 |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| attn/Wang | theta=1 | 9,873 | 8,552 | 0.809 | 0.934 | 0.867 | 0.229 | 0.454 | 0.305 |
| attn/Wang | Poisson th=0.02 | 7,998 | 8,552 | 0.917 | 0.858 | **0.887** | 0.254 | 0.127 | 0.170 |
| mamba4/Wang | theta=1 | 9,832 | 8,552 | 0.810 | 0.932 | 0.867 | 0.227 | 0.441 | 0.299 |
| mamba4/Wang | Poisson th=0.02 | 7,869 | 8,552 | 0.918 | 0.844 | **0.879** | 0.271 | 0.135 | 0.180 |
| attn/Janich | theta=1 | 10,263 | 9,850 | 0.880 | 0.917 | **0.898** | 0.467 | 0.530 | 0.497 |
| attn/Janich | Poisson th=0.05 | 8,367 | 9,850 | 0.960 | 0.816 | 0.882 | 0.570 | 0.176 | 0.269 |
| mamba4/Janich | theta=1 | 9,731 | 9,850 | 0.899 | 0.888 | **0.893** | 0.493 | 0.446 | 0.468 |
| mamba4/Janich | Poisson th=0.1 | 8,267 | 9,850 | 0.948 | 0.796 | 0.865 | 0.516 | 0.162 | 0.246 |
| attn/GSE243134 | theta=1 | 12,200 | 11,761 | 0.879 | 0.912 | **0.895** | 0.461 | 0.483 | 0.472 |
| attn/GSE243134 | Poisson th=0.05 | 10,093 | 11,761 | 0.950 | 0.815 | 0.878 | 0.543 | 0.134 | 0.216 |
| mamba4/GSE243134 | theta=1 | 12,177 | 11,761 | 0.881 | 0.912 | **0.896** | 0.466 | 0.481 | 0.473 |
| mamba4/GSE243134 | Poisson th=0.05 | 10,238 | 11,761 | 0.946 | 0.823 | 0.880 | 0.516 | 0.139 | 0.218 |
| attn/Ruiz-Orera | theta=1 | 14,641 | 13,147 | 0.832 | 0.926 | **0.876** | 0.385 | 0.551 | 0.453 |
| attn/Ruiz-Orera | Poisson th=0.05 | 11,329 | 13,147 | 0.930 | 0.801 | 0.861 | 0.448 | 0.154 | 0.229 |
| mamba4/Ruiz-Orera | theta=1 | 14,614 | 13,147 | 0.836 | 0.929 | **0.880** | 0.402 | 0.567 | 0.470 |
| mamba4/Ruiz-Orera | Poisson th=0.05 | 11,464 | 13,147 | 0.926 | 0.808 | 0.863 | 0.459 | 0.167 | 0.245 |

- **Poisson is a precision lever, not a free win.** Overall precision rises on all eight dumps (+0.07
  to +0.11, up to 0.960), but overall F1 improves on only ONE of the four datasets -- Wang, the one
  where theta=1 was over-calling badly (precision 0.809). On Janich, GSE243134 and Ruiz-Orera the
  calibrated arm is 0.013-0.028 F1 WORSE. The earlier project shorthand that Poisson is "the winning
  calibration lever" was generalised from the Wang-like case; it does not hold at 3 of 4 datasets.
- **It costs non-canonical recall by roughly 3x, every time.** ncR falls 0.454 -> 0.127, 0.530 ->
  0.176, 0.483 -> 0.134, 0.551 -> 0.154. Non-canonical F1 drops on all eight, typically by half.
  Non-canonical *precision* does improve (0.467 -> 0.570 on attn/Janich), so the dial is real -- it
  simply sits at a point tuned for annotated CDS, which is what the anchor is made of.
- **Practical reading:** use the calibrated arm when reporting canonical-CDS agreement or building a
  proteogenomics database where precision dominates; use theta=1 when the goal is non-canonical
  discovery and a downstream filter (the O3 FP-filter) will do the pruning. Report both, always.
- attn and mamba4 remain interchangeable: max separation 0.008 F1 at theta=1, and the two models even
  pick different theta* on Janich (0.05 vs 0.1) while landing within 0.017 F1 of each other.

### 3. One scoring convention, enforced by imports rather than by care

The first version of `score_released_two_arm.py` loaded the collapsed files directly and reported
attn/Wang theta=1 at F1 0.725, against 0.867 from the established path -- same dump, same arm. The
difference was entirely convention: `compare_dropin_calls.py` restricts to the dump's test
transcripts, drops ORFs under 90 nt, and requires predicted calls to reach 0.5x uniform mean density,
and none of that was reimplemented. Rather than copy the three constants, the loader was extracted
from `compare_dropin_calls.main()` into `compare_dropin_calls.build_loader()` and imported. Copied
constants drift; an import cannot.

### 4. The same staleness was in five of the six main Figure 1 panels

Finding the held-out numbers stale prompted a scan of EVERY figure generator rather than the ones that
seemed likely (feedback_blast_radius_scan_by_artifact: enumerate the artifact, not the datasets you
remember). Five of the six Fig 1 panels read `orf_v2_attn_onehot_holdout_Hepatocytes`. Only
`arch_attn` was on a released model.

| panel | state before | now |
|---|---|---|
| A1 localization ceiling | stale | rebuilt from `eval_localization_released.sbatch` |
| A2 drop-in | stale, F1 0.923 | released, F1 **0.909**; asserts the npz `meta` tag before plotting |
| A3 LOTO spread | pre-union 9-fold | **unchanged, labelled** -- no union 9-fold exists, redoing it is 9 retrains |
| B4 replicate ceiling | hardcoded pre-union dict | reads the released run's `extra_metrics.json` |
| B5 expression independence | stale | rebuilt from the released `localization_metrics.json` |
| C10 saliency | stale ckpt + Kozak track | released ckpt + union no-Kozak track (both had to move together) |
| prediction_examples | stale | released |

C10's conclusion survived the move (start-codon saliency rank 3/2036 uORF and 10/7914 lncRNA on the
released config, against 8 and 15 before -- it got slightly stronger).

**A second-order error, worth recording because the guard looked right.** The first attempt at the
released localization eval passed `data/packed/orf_track_v2_nokozak.npy` and asserted `kozak: none`.
The assertion PASSED -- that file is kozak=none. It is also the FIBROBLAST-universe track (979,634,968
bytes) rather than the union one (2,126,471,848). The guard tested the correct property on the wrong
axis, and both tracks are self-consistent with their own packs, so a mis-tracked run would have
produced a plausible number rather than an error. Both jobs were cancelled and resubmitted with the
pack suffix, ORF track and one-hot FASTA moved together, plus a byte-size assertion against the
training track. The same defect was present in C10 and was fixed there too.

The durable rule, now in `figures/README.md`: a generator that quotes a model number either asserts
the checkpoint tag or reads it from a file and prints it. Hardcoded model numbers are how this
happened, and directory names cannot be trusted because they encode the dataset, not the checkpoint.

Two smaller fixes in the same pass:

- `compare_dropin_calls.py --official` is now optional. It was `required=True`, but only Wang and
  GSE243134 have a call set from RiboCode over the full observed data outside the pack; Janich does
  not. The two `*_vs_official` rows are harness validation, not the primary metric, so omitting
  official now drops those rows instead of blocking the whole comparison.
- **Figure B4 compared two different quantities.** Its third bar put a per-codon CDS *elongation*
  ceiling (0.9501) against a *periodicity* score (0.526); `replicate_concordance.py` records the model
  side of that stratum as `None` / "our per-codon TBD", so it had never been computed and the derived
  "55% of ceiling" was meaningless. Replaced with the two strata the concordance script was written to
  pair with the eval (`pc_uorf5` / `pc_dorf3`). The model values are also no longer hardcoded: B4 now
  reads them from the released run's `extra_metrics.json` and writes the run name into
  `B4_values.json`. Released numbers: 70% (pc whole-tx), 51% (lncRNA), 61% (uORF 5'UTR), 34% (dORF
  3'UTR) of the reproducible ceiling. The dORF ceiling is itself only 0.749, so part of the weak dORF
  result is the assay, not the model.

## Task 61: the no-RNA-seq ablation on the DEPLOYED recipe (2026-08-08, LANDED)

`scripts/train_union_inputablation.sbatch`, two arms, cloned line-for-line from
`train_loto_union.sbatch` with `--input_mode` the only difference. The `both` arm is the deployed run
itself, not a re-train, so the contrast carries no extra seed noise. Held-out Hepatocytes, n=70,883
scored transcripts in all three.

| input_mode | pc profile r | lncRNA profile r | pc count r | pc periodicity (pred) |
|---|--:|--:|--:|--:|
| **both** (deployed) | 0.6699 | 0.4583 | 0.8990 | 0.3359 |
| **emb** -- sequence + ORF track, RNA-seq ZEROED | 0.6601 | 0.4462 | **0.7063** | 0.3704 |
| **cov** -- RNA-seq only, sequence + ORF ZEROED | 0.1226 | 0.0838 | 0.8078 | -0.0155 |

The double dissociation from Task 17A reproduces on the union recipe, but **sharper, and it revises
what the cell-type-specificity claim should say**:

- **Profile SHAPE is essentially sequence-borne.** Zeroing RNA-seq costs 0.0098 of pc profile Pearson:
  sequence-only recovers **98.5%** of the deployed model's shape (0.6601 of 0.6699). On the pre-union
  recipe that figure was 91%; the bigger universe made the model *more* sequence-driven, not less.
  Periodicity is not merely preserved without RNA-seq, it sharpens (0.3704 vs 0.3359) -- the coverage
  channel contributes a smooth envelope that slightly dilutes fine periodicity.
- **Sequence alone cannot do magnitude.** Count Pearson falls 0.8990 -> 0.7063 (**-0.193**) when RNA-seq
  is removed, and RNA-seq ALONE reaches 0.8078 -- better per-transcript depth than sequence alone.
- **RNA-seq alone cannot do shape at all.** 0.1226 pc profile Pearson with periodicity destroyed
  (-0.0155): a smooth magnitude signal cannot express a periodic one.

### What this means for the cell-type-specificity claim

The ablation was set up to test whether the model is a cell-type-invariant sequence prior in disguise.
The answer is not a clean yes or no, and the manuscript wording has to follow the measurement:

**The model's cell-type specificity runs through the count head, not the profile shape.** Handed a
different cell type's RNA-seq, the predicted SHAPE over a given transcript barely moves (98.5% of it is
recoverable with no RNA-seq at all), while the predicted DEPTH moves a lot (-0.193 count Pearson
without it). That is still a genuine cell-type-specific mechanism for ORF CALLING -- a transcript not
expressed in the query cell type gets low predicted depth, fails the caller's significance test, and is
not called -- but it is a claim about which ORFs clear threshold, not about the shape of the signal
over them.

So: "the model predicts cell-type-specific translation" is supportable; "the model predicts a
cell-type-specific profile shape" is not, and should not be written. The proteogenomics argument
depends on the former, so it stands; it should cite the count-head number (-0.193), not a shape number.

Recorded as the direct evidence for FIGURES_PLAN item 1. Both ablation runs are at
`results/ablation/orf_v2_attn_onehot_union_nokozak_input_{emb,cov}_holdout_Hepatocytes`.

## Task 68: does better RNA-seq change ORF-call precision and recall? (2026-08-08, LANDED)

Ribo-seq held FIXED across all three arms (the same 5 Janich P-site hd5, read out of arm a's
provenance so the arms cannot drift); only the RNA-seq changes. Design and the PRJEB34766/PRJEB86747
comparator reasoning are in methods.md 5S. Scored by `scripts/score_rna_quality_factorial.py` through
the same `compare_dropin_calls.build_loader` filtering as every other drop-in number.

| arm | model | variant | nTx | nPred | nRef | P | R | F1 | ncP | ncR | ncF1 |
|---|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| a Janich/Janich | attn | obsdepth | 14,887 | 9,518 | 9,850 | 0.935 | 0.904 | 0.919 | 0.599 | 0.493 | 0.541 |
| a Janich/Janich | attn | preddepth | 14,887 | 10,263 | 9,850 | 0.880 | 0.917 | 0.898 | 0.467 | 0.530 | 0.497 |
| a Janich/Janich | mamba4 | obsdepth | 14,887 | 9,476 | 9,850 | 0.944 | 0.908 | 0.926 | 0.643 | 0.499 | 0.562 |
| a Janich/Janich | mamba4 | preddepth | 14,887 | 9,731 | 9,850 | 0.899 | 0.888 | 0.893 | 0.493 | 0.446 | 0.468 |
| b Janich/34766 | attn | obsdepth | 14,887 | 9,315 | 9,850 | 0.952 | 0.900 | 0.925 | 0.677 | 0.471 | 0.555 |
| b Janich/34766 | attn | preddepth | 14,887 | 10,127 | 9,850 | 0.891 | 0.916 | 0.903 | 0.497 | 0.519 | 0.508 |
| b Janich/34766 | mamba4 | obsdepth | 14,887 | 9,290 | 9,850 | 0.955 | 0.901 | 0.927 | 0.693 | 0.461 | 0.553 |
| b Janich/34766 | mamba4 | preddepth | 14,887 | 10,112 | 9,850 | 0.891 | 0.915 | 0.903 | 0.498 | 0.513 | 0.506 |
| c 34766/34766 | attn | obsdepth | 16,474 | 10,866 | 11,547 | 0.953 | 0.897 | 0.924 | 0.673 | 0.457 | 0.544 |
| c 34766/34766 | attn | preddepth | 16,474 | 11,840 | 11,547 | 0.896 | 0.919 | 0.907 | 0.511 | 0.530 | 0.520 |
| c 34766/34766 | mamba4 | obsdepth | 16,474 | 10,863 | 11,547 | 0.956 | 0.899 | 0.927 | 0.686 | 0.453 | 0.546 |
| c 34766/34766 | mamba4 | preddepth | 16,474 | 11,768 | 11,547 | 0.900 | 0.917 | 0.908 | 0.521 | 0.519 | 0.520 |

### The COVERAGE path (a vs b) buys precision, and it is small

Universe held byte-identical, so this is a clean like-for-like contrast. All four model-by-variant
cells move the same way on precision and barely at all on F1:

| model | variant | dF1 | dPrecision | dRecall | d nc-F1 |
|---|---|--:|--:|--:|--:|
| attn | obsdepth | +0.006 | +0.017 | -0.004 | +0.014 |
| attn | preddepth | +0.005 | +0.010 | -0.001 | +0.011 |
| mamba4 | obsdepth | +0.001 | +0.011 | -0.007 | -0.008 |
| mamba4 | preddepth | +0.010 | -0.007 | +0.027 | +0.037 |

**Swapping to 15x deeper, paired, poly(A) RNA-seq coverage is worth +0.001 to +0.010 F1.** That is a
real but modest effect, and it is concentrated in precision (+0.010 to +0.017 on three of four cells).
The one cell that gains through recall instead, mamba4 preddepth, is also the one that started worst
(0.888 recall, the lowest of the eight arm-a numbers), so it had the most to recover.

This is the result the B6 input ablation predicts. B6 showed the coverage channel drives MAGNITUDE and
not shape (sequence alone recovers 98.5% of profile Pearson; the count head loses 0.193 without
RNA-seq). If shape is nearly RNA-independent, then improving the coverage track cannot move the shape
much -- and it does not: the obsdepth arm, which borrows real depth and therefore isolates shape,
moves +0.001 to +0.006.

### The UNIVERSE path is where the volume is

| | arm b | arm c | change |
|---|--:|--:|--:|
| scored transcripts | 14,887 | 16,474 | +1,587 (+10.7%) |
| reference ORF calls | 9,850 | 11,547 | +1,697 (**+17.2%**) |
| model predictions (preddepth, mamba4) | 10,112 | 11,768 | +1,656 |

Better RNA-seq reveals **17% more real ORFs** while adding only 10.7% more transcripts -- the
transcripts a deeper poly(A) library rescues from below the TPM >= 1 cut are ORF-denser than average.
The model keeps pace: +1,656 predictions against +1,697 new references. Arm c's F1 is NOT comparable
to arm b's (the reference set changed with the universe), which the scorer prints as a separate block
so the delta is not read as an effect size.

### The two models converge once the RNA is good

Full swap, a -> c, standalone arm:

| model | F1 | non-canonical F1 |
|---|---|---|
| attn | 0.898 -> 0.907 (+0.009) | 0.497 -> **0.520** (+0.023) |
| mamba4 | 0.893 -> 0.908 (+0.015) | 0.468 -> **0.520** (+0.052) |

Both land on the same non-canonical F1 to three decimals, having started 0.029 apart. On this axis
**RNA quality is a larger lever than the architecture choice**, and mamba4's larger gain is recovery,
not superiority. This is worth stating next to decision D1b, whose own caveat is that mamba4's
profile-Pearson edge does not carry into downstream tasks: here neither model's identity survives a
change of RNA-seq.

### Practical reading

- If the goal is non-canonical discovery, spend effort on the RNA-seq library before the architecture:
  paired poly(A) at depth adds 17% more findable ORFs and +0.023 to +0.052 non-canonical F1.
- The gain arrives mostly through the UNIVERSE (what is scoreable at all), not through the coverage
  track the model reads. A cheaper library that still clears TPM >= 1 on the same transcripts would
  capture most of the coverage-path benefit, which is small.
- Do not expect better RNA-seq to improve the predicted profile SHAPE. It does not, and B6 explains
  why.

## Tasks 62 and 64: the immunopeptidome panel is closed at 5, and the main figures are assembled (2026-08-09)

### Task 62 SUPERSEDED at 5 of 6 datasets

Closed by user instruction rather than completed. Locked decision D2 asked for 6 immunopeptidome
datasets; **5 were delivered** (A549 tryptic, plus HBL-1 / SU-DHL-4 / DoHH2 / THP-1 HLA-I), all on the
released models with frozen search parameters. Five is enough for the Fig 2a forest plot to demonstrate
a reproducible pattern rather than a one-off: 5 datasets x 2 models x 2 arms = 20 points, every one
above 1.0, median 93.6x.

The two not delivered are blocked on **data availability, not effort**, and both were checked:

- **B721.221 MS** -- MassIVE FTP port 21 is filtered from prism; three unblock routes are recorded in
  DATA_PROVENANCE.md. B721 is not wasted: it supplies the drop-in GROUND-TRUTH check (F1 0.666 = 75% of
  the 0.889 split-half reproducibility ceiling), which is the more valuable use of that dataset anyway,
  since it is the only place in the project where predicted ORF calls meet measured ones.
- **A mouse immunopeptidome** -- candidates rejected 2026-08-07: PXD008733 (Schuster murine tissue map)
  has 39 raw files but NO matched RNA-seq; Rospo CT26 MMRd has ENA RNA-seq but NO raw MS accession.
  Matched RNA-seq is non-negotiable because it IS the model input.

**Consequence, and it must be stated in the manuscript rather than left implicit:** Figure 2 has no
cross-species MS arm. FIGURES_PLAN always listed that as a stretch goal. The cross-species claim rests
on Ribo-seq -- the mouse Wang / Janich / GSE243134 drop-ins -- and not on mass spectrometry.

### Task 64: main figures assembled

`figures/main/make_main_figures.py` composites the per-panel PDFs vector-preserving (pypdf; panels are
placed as scaled PDF pages, not re-rastered) and writes a manifest naming every source panel and its
mtime. Layout is a data structure, so re-specifying a figure is a one-list edit.

- **Figure 1 (validity)**, 7.20 x 10.19 in, 5 panels: (a) observed-vs-predicted profile exemplars,
  (b) drop-in ORF calling F1 0.913, (c) localization AUROC vs the observed Ribo-seq ceiling, (d) depth
  crossover, (e) 9-fold LOTO spread.
- **Figure 2 (utility)**, 7.20 x 8.67 in, 3 panels: (a) discovery-density forest across 5 datasets,
  (b) model vs CPAT/CPC2 split by assay, (c) DB-design tradeoff.

Deviations from the plan, all deliberate: the plan's Fig 1d ("one-hot ~= token-emb") is omitted because
no standalone panel was ever built for it and its plan numbers are pre-union; A2 takes that slot. Fig 2
has no D14 highlighted-PSM panel, which was gated behind task 62.

Five known limitations are written into `figures/main/FIGURE_DATA_INPUTS.md` rather than left for a
reader to discover -- the most important being that **panel 1e is still on the pre-union recipe** (no
union 9-fold exists; building one is nine retrains) and that **Figure 1 at 10.2 in exceeds one printed
page**, so either the format grows or a panel moves to supplemental. That is a layout decision, not a
data one, and is left open.

While assembling, one more instance of the D1b problem surfaced and was fixed: `prediction_examples`
(panel 1a) was still on attn. It only READS the dump, so switching it to mamba4 needed no GPU despite
mamba4 inference requiring CUDA. It now takes `FIG_MODEL` like the other four.

## Mouse-liver 3x3 factorial (COMPLETE, 2026-08-10)

> **RESCORED CANONICALLY 2026-08-13 -- conclusions unchanged.** The original numbers were computed on
> packs built without `--alignEndsType EndToEnd` and without the ncRNA + cross-gene filter, both
> mandatory under `docs/PIPELINE_POLICY.md`. Everything below was rebuilt from re-aligned reads:
> 40 Ribo alignments -> 3 pools -> 9 packs -> 18 dumps -> rescore. Final-recipe pooling removes 3.5-4.9%
> of P-sites (janich 110.0M -> 106.1M, wang 78.3M -> 74.5M, gse243134 150.5M -> 143.2M). The
> off-recipe tables are archived at `results/_archive_offrecipe_2026_08_13/` for comparison; current
> tables are `results/mouse_liver_3x3_canon/scored/`.
>
> What moved: observed call sets fell ~4-6% (see below); `pred_preddepth` F1 fell ~0.01.
> What did not: every qualitative conclusion, including the headline matched-vs-mismatched RNA effect
> (mean delta +0.0030 off-recipe -> **+0.0035** canonical, negative in the same 2 of 12 cells).

3 Ribo-seq datasets x 3 RNA-seq inputs x 2 models, all pooled from BAMs through ONE pipeline onto one
shared universe (22,974 tx). 18/18 cells, both prediction arms each. Genomic keying
(gene_id, ORF_gstop), `compare_dropin_calls.build_loader` filtering. Tables:
`results/mouse_liver_3x3_canon/scored/`.

Observed call sets (post-filter), canonical vs the archived off-recipe values:

| Ribo dataset | off-recipe | canonical | change |
|---|--:|--:|--:|
| janich    | 12,804 | 12,289 | -4.0% |
| gse243134 | 13,145 | 12,584 | -4.3% |
| wang      | 12,104 | 11,646 | -3.8% |

Fewer calls is the expected direction: EndToEnd refuses the soft-clipped alignments that previously
contributed spurious P-sites, and the cross-gene filter drops reads that were being counted against
transcripts of other genes. Both remove signal that should never have been there.

**Construction self-check PASSES.** `real` calls are identical across each Ribo row whatever RNA sat
in the pack, and `pred_preddepth` is identical down each RNA column. Both follow from the design and
are asserted, not assumed.

### A. Ribo-vs-Ribo ceiling (no model involved)

| reference | compared | precision | recall | F1 all-class | F1 annotated |
|---|---|---|---|---|---|
| janich | gse243134 | 0.941 | 0.966 | 0.954 | 0.995 |
| janich | wang | 0.971 | 0.918 | 0.944 | 0.993 |
| gse243134 | wang | 0.975 | 0.898 | 0.935 | 0.992 |

Two real experiments agree at F1 0.935-0.954 all-class but 0.992-0.995 on annotated ORFs. Essentially
all cross-experiment disagreement is in non-canonical ORFs.

### B. Model vs observed (n=18 cells per arm) -- CANONICAL

| arm | F1 all-class | F1 annotated |
|---|---|---|
| pred_obsdepth (predicted shape, reference depth) | 0.897-0.925 | 0.982-0.990 |
| pred_preddepth (standalone, no Ribo-seq at inference) | 0.855-0.871 | 0.966-0.980 |

On annotated ORFs the model sits at the experimental ceiling (0.982-0.990 vs 0.991-0.995). The fully
standalone arm is within ~0.02 of two real experiments agreeing with each other. Arm structure is
consistent throughout: obsdepth precision-heavy (P~0.93 / R~0.90), preddepth recall-heavy
(P~0.82 / R~0.91) -- the documented standalone over-calling.

### C. Does matching the RNA input to the Ribo dataset help? No. -- CANONICAL

Mean delta **+0.0035 F1**, range **-0.0042 to +0.0070**, positive in 10 of 12 rows, every
|delta| <= 0.007. Negative for mamba4/Wang in both arms.

| model | ribo reference | arm | F1 matched | F1 mismatched (mean) | delta |
|---|---|---|---|---|---|
| mamba4 | gse243134 | pred_preddepth | 0.8669 | 0.8599 | +0.0070 |
| mamba4 | gse243134 | pred_obsdepth | 0.9126 | 0.9060 | +0.0066 |
| mamba4 | janich | pred_obsdepth | 0.9192 | 0.9127 | +0.0065 |
| attn | wang | pred_obsdepth | 0.9193 | 0.9140 | +0.0053 |
| attn | wang | pred_preddepth | 0.8665 | 0.8616 | +0.0049 |
| attn | gse243134 | pred_obsdepth | 0.9057 | 0.9012 | +0.0045 |
| attn | janich | pred_preddepth | 0.8688 | 0.8649 | +0.0038 |
| attn | gse243134 | pred_preddepth | 0.8657 | 0.8624 | +0.0033 |
| attn | janich | pred_obsdepth | 0.9105 | 0.9075 | +0.0030 |
| mamba4 | janich | pred_preddepth | 0.8678 | 0.8650 | +0.0028 |
| mamba4 | wang | pred_obsdepth | 0.9221 | 0.9240 | -0.0019 |
| mamba4 | wang | pred_preddepth | 0.8634 | 0.8676 | -0.0042 |

*Source: `results/mouse_liver_3x3_canon/scored/`. The off-recipe version of both tables (mean
+0.0030, range -0.0057 to +0.0079, arm F1 0.895-0.926 / 0.864-0.880) is archived at
`results/_archive_offrecipe_2026_08_13/mouse_liver_3x3/scored/`. The conclusion is identical; only
the third decimal moved.*

### Limitation

metaplots re-derives periodic read lengths per sample, so the rebuilt Janich pack differs ~1.7% from
its historical one. These nine cells are internally consistent but NOT comparable to pre-2026-08-10
Janich numbers, which are archived at `results/_archive_pre_2026_08_10_janich_pack/`.

## Human ORF calls: THP-1 and CAR-T (2026-08-11)

Two human Ribo-seq datasets pulled, processed and ORF-called end to end. Each is scored
against RiboCode run on its OWN observed P-sites, on its OWN expression-restricted universe,
so the transcript space is identical on both sides by construction. Genomic keying,
`build_loader` filtering, both prediction arms. Tables: `results/human_orf_calls/scored/`.

| dataset | universe | model | arm | P | R | F1 | n_ref |
|---|--:|---|---|--:|--:|--:|--:|
| cart | 23,887 tx | attn | `pred_obsdepth` | 0.934 | 0.921 | 0.927 | 8,856 |
| cart | 23,887 tx | attn | `pred_preddepth` | 0.888 | 0.922 | 0.905 | 8,856 |
| cart | 23,887 tx | mamba4 | `pred_obsdepth` | 0.932 | 0.926 | 0.929 | 8,856 |
| cart | 23,887 tx | mamba4 | `pred_preddepth` | 0.890 | 0.924 | 0.907 | 8,856 |
| gse208041 | 39,611 tx | attn | `pred_obsdepth` | 0.953 | 0.857 | 0.903 | 14,304 |
| gse208041 | 39,611 tx | attn | `pred_preddepth` | 0.824 | 0.879 | 0.851 | 14,304 |
| gse208041 | 39,611 tx | mamba4 | `pred_obsdepth` | 0.954 | 0.861 | 0.905 | 14,304 |
| gse208041 | 39,611 tx | mamba4 | `pred_preddepth` | 0.821 | 0.880 | 0.850 | 14,304 |

**CAR-T scores higher than THP-1 on both arms** (F1 0.929 vs 0.905 obsdepth, 0.907 vs 0.851
standalone), despite a smaller universe and shallower Ribo-seq. The standalone gap is the
striking part: CAR-T loses only 0.022 F1 going fully Ribo-seq-free, THP-1 loses 0.054.

### Per ORF class (`pred_obsdepth`, both models)

| class | GSE208041 F1 | CAR-T F1 | n_ref (208041 / CAR-T) |
|---|--:|--:|--:|
| annotated | 0.991-0.992 | 0.991-0.993 | 10,654 / 7,399 |
| uORF | 0.635-0.636 | 0.679-0.685 | 1,191 / 579 |
| novel | 0.650-0.666 | 0.615-0.636 | 1,347 / 297 |
| internal | 0.154-0.234 | 0.186-0.213 | 138 / 87 |

**Human internal ORFs are NOT at the floor.** They score 0.154-0.241 here against 0.000-0.037
in the mouse liver panel. The near-total failure on internal ORFs is specific to that mouse
panel, not a general property of the model, and the mouse figure should not be generalised.

### GSE39561: closed as a negative result

The third dataset has NO usable Ribo-seq signal and is excluded from the tables above.
Measured frame concentration against annotated CDS starts is 40.7% pooled (best read length
54.5%) against a 33.3% random baseline; GSE208041 reaches 91.2% at its dominant length. The
cause is ragged nuclease digestion: inserts smear across 26-32 nt with no single length above
18%, where GSE208041 concentrates 51.6% at 28 nt. With no reproducible 5'-end-to-P-site
offset, RiboCode cannot assign P-sites, and its refusal is correct rather than a threshold to
override. Adapter and trim were verified correct beforehand (poly-A tail, no ligated adapter
in the reads at all), and the library is not noise -- 40.7% over a 33.3% baseline means real
footprints are present, just not resolvable to a frame.

It was built as a coverage-only pack and does produce standalone predictions, but those are
IDENTICAL to GSE208041's (verified by call-set comparison) because it borrows that dataset's
RNA and universe. It is not an independent third measurement and must not be reported as one.
Its remaining value is as a demonstration that the model predicts where RiboCode cannot.


## Canonical redo of the interpretability chain (2026-08-13, LANDED)

The three "what did the model learn" analyses (Phases 3 and 4) were first run on off-recipe packs.
After the final-recipe rebuild they were re-run end to end on final-recipe substrate, writing to `*_canon`
paths alongside the originals so every finding could be compared rather than silently replaced.
Expectations were stated in advance, so the comparison is a check and not a rationalisation:
saliency should be unchanged (its gradient depends on sequence, ORF track and coverage, none of
which moved), ablation signs should hold, and codon was genuinely uncertain because its empirical
side rebuilds from shifted P-sites.

**Codon occupancy: unchanged, and that is the result.** Medians are identical to three decimals
(within-species ceiling 0.944, across-species 0.269, mouse 0.263, human 0.548). This is not a
skipped rerun: the canonical and off-recipe tables are different files, with `max|diff| = 0.0385`
per codon. They correlate at **r = 0.99952**. Per-transcript normalization to the CDS mean cancels a
roughly uniform ~5% change in P-site count, so the codon-occupancy vector is insensitive to the
alignment recipe. Worth stating in the manuscript: it means the codon result does not depend on
which of the two recipes a reader would have used.

The substantive codon finding therefore stands. Mouse (transcripts never seen in training) reaches
r = 0.263 against a cross-species empirical ceiling of 0.269 -- **98% of the achievable ceiling**.
Human reaches 0.548, higher in absolute terms but on transcripts that were 99.7% present in
training, so the two numbers measure different things and should not be read as a ranking.

**Saliency: unchanged, as predicted.** Class ordering by median `start_rank` is identical, with uORF
sharpest and annotated CDS weakest:

| ORF class | canonical start_rank | off-recipe | in_orf_ratio (canon) | seq_frac (canon) |
|---|--:|--:|--:|--:|
| uORF         |  35 |  33 | 9.59 | 0.635 |
| Overlap_uORF |  36 |  38 | 7.13 | 0.613 |
| internal     |  74 |  77 | 6.07 | 0.596 |
| novel        |  81 |  69 | 2.63 | 0.538 |
| dORF         |  84 |  86 | 3.48 | 0.816 |
| Overlap_dORF |  99 |  89 | 3.49 | 0.787 |
| annotated    | 160 | 154 | 1.28 | 0.602 |

**Channel ablation: signs hold, 57 of 60 comparisons agree.** The three disagreements are all cells
where both values sit within 0.005 of zero, i.e. noise around no-effect rather than a reversal. The
two conclusions that matter both reproduce:

| finding | off-recipe | canonical |
|---|--:|--:|
| frame-occupancy channels HARM `internal` (delta when removed, range over 12 arms) | +0.006 to +0.050 | +0.006 to +0.053 |
| start channel (ch3) is critical for `uORF` (attn / mamba4) | -0.075 / -0.090 | -0.069 / -0.073 |
| removing all three frame channels (f012) costs `ALL` (attn / mamba4) | -0.024 / -0.026 | -0.024 / -0.024 |

Removing a frame-occupancy channel makes the model call `internal` ORFs *better*. The channels
encode where the annotated reading frame is, which is exactly the prior that suppresses a
plausible ORF sitting inside a known CDS in a different frame. The `ch3` result is the mirror image:
the start-propensity channel is nearly free for annotated CDS (delta ~0.000) and costly only for
uORFs, whose defining evidence is a start site upstream of the annotated one.

Scoring for the ablation previously existed only as an inline shell heredoc, which
`docs/PIPELINE_POLICY.md` forbids. It is now `scripts/score_channel_ablation.py`, importing the same
`compare_dropin_calls.build_loader` filtering as every other drop-in number in the project.

## Regeneration divergence separates into a depth noise floor and a recipe term (2026-08-13)

The three mouse leukocyte packs were rebuilt on final-recipe alignments. They were expected to be
off-recipe; reading `scripts/heldout/align_and_filter.sh` showed they never were -- it already applied
EndToEnd, mm1 and `filter_tx_heldout`, and the filter received query-grouped input because STAR's
`Aligned.toTranscriptome.out.bam` is read-ordered regardless of `--outSAMtype BAM SortedByCoordinate`
(which sorts only the genome BAM). Their sole difference from canonical is two tightened STAR filters
that the 13-arm HUVEC sweep found immaterial.

Rebuilding them anyway turned an assumption into a measurement, and gave three near-same-recipe
regeneration numbers. All comparisons use the same metric as the HUVEC test,
`sum|new - old| / sum(old)` over per-nt P-sites on the pack's own universe and row order:

| pack | P-sites | net change | \|divergence\| | per-nt r | RNA coverage |
|---|--:|--:|--:|--:|---|
| gse155087_tcell | 77,961,556 | +0.13% | **0.13%** | 0.9986 | byte-identical |
| gse120762_lps | 38,481,730 | +0.40% | **0.47%** | 0.9901 | byte-identical |
| gse120762_nt | 14,253,922 | +3.08% | **3.16%** | 0.9668 | byte-identical |

Coverage came back byte-identical (`np.array_equal`) for all three, which independently confirms the
RNA input sets matched the originals -- the check that mattered, since a glob bug nearly fed Ribo
alignments in as RNA coverage on this exact path.

**Same-recipe divergence is a depth-dependent noise floor.** Across the three it scales inversely
with library depth (0.13% at 78M P-sites, 0.47% at 38M, 3.16% at 14M). The mechanism is RiboCode
`metaplots` re-deriving P-site offsets per sample: the offset choice is mildly stochastic, and its
effect on pooled totals shrinks as depth grows. All three changes are POSITIVE, consistent with
canonical being the looser recipe on the two swept mismatch filters.

**HUVEC's 3.04% is larger than any of these**, and HUVEC is the only one of the four compared against
a pack built by a genuinely different pipeline. That is consistent with a recipe term on top of the
regeneration floor.

> **SUPERSEDED in part, 2026-08-13.** An earlier version of this paragraph went further and claimed
> HUVEC's excess was "roughly 20x what the same-recipe trend predicts at that depth", extrapolating
> the n=3 depth trend above. The full 8-tissue Chothani table (next section) shows **no depth
> relationship at all**: the deepest library (Hepatocytes, 1.16B P-sites) and the shallowest (HUVEC,
> 96M) both diverge at ~3.1%, across a 12x depth range. The depth ordering in the three leukocyte
> packs does not generalise, and it should not be used to interpret any other number. What survives
> is the narrow statement: those three same-recipe rebuilds all came in under 3.2%, and they are the
> only same-recipe comparison available.

**Caveat.** n = 3 for the depth trend, spanning a 5.5x range of depth, and the three differ in tissue
and protocol as well as depth. It is a consistent ordering, not a fitted relationship, and it should
not be extrapolated. The load-bearing claim is narrower: at 96M P-sites the same-recipe floor is well
under 1%, so HUVEC's 3.04% is not regeneration noise.

## Task 79: the reproducibility test on all 8 Chothani tissues (2026-08-13, LANDED)

The training-data reproducibility claim previously rested on ONE tissue (HUVEC, 3.04%). All 74 Ribo
runs were re-aligned canonically (`scripts/riboseq_align.sbatch`, 71 new + 3 existing, 0 failures,
filter drop median 8.86%) and each tissue re-pooled and compared against the pack the model actually
trained on. Metric is `sum|new - old| / sum(old)` over per-nt P-sites on each pack's own universe --
the same metric as the original HUVEC test, which the generalised script reproduces EXACTLY (3.04%,
r=0.9567 vs the original 3.04%/0.957). That exact reproduction is the reason the other seven are
trustworthy.

| tissue | pack | n bam | training P-sites | canonical | net | \|divergence\| | per-nt r |
|---|---|--:|--:|--:|--:|--:|--:|
| Hepatocytes | `packed_union_Hepatocytes` | 5 | 1,159,703,745 | | | **3.16%** | 0.9679 |
| Fibroblast | `packed_union` | 32 | 1,093,425,533 | | | **3.71%** | 0.9374 |
| HCAEC | `packed_union_HCAEC` | 5 | 714,826,272 | | | **3.20%** | 0.9185 |
| Fat | `packed_union_Fat` | 6 | 581,557,882 | | | **11.66%** | 0.8442 |
| ES | `packed_union_ES` | 6 | 508,576,671 | | | **4.01%** | 0.9800 |
| VSMC | `packed_union_VSMC` | 11 | 486,734,365 | | | **5.85%** | 0.8373 |
| HA_EC | `packed_union_HA_EC` | 6 | 173,184,690 | | | **8.01%** | 0.4320 |
| HUVEC | `packed_union_HUVEC` | 3 | 95,959,146 | | | **3.04%** | 0.9567 |

Full values incl. net change: `results/chothani_regeneration/divergence_by_tissue.{tsv,md}`.

**Median 3.86%, range 3.04-11.66%; per-nt r 0.432-0.980.**

**The headline: one tissue was not enough, and it was the flattering one.** HUVEC's 3.04% is the
LOWEST divergence of all eight. Quoting it alone understates the typical case by about a third
(median 3.86%) and understates the worst case by nearly 4x (Fat, 11.66%). The eight-tissue table is
what belongs in the manuscript.

**Divergence and correlation are partly decoupled**, which means they are diagnosing different
things and both belong in the table. Fat moved the most signal (11.66%) but moved it coherently
(r = 0.844); HA_EC moved less (8.01%) but incoherently (r = 0.432). A single summary statistic would
hide that.

**No depth relationship.** Hepatocytes (1.16B P-sites) and HUVEC (96M) both diverge at ~3.1% across a
12x depth range, and the worst tissue (Fat, 582M) sits in the middle of the depth ordering. The
inverse-depth pattern seen in the three leukocyte rebuilds does not generalise; the note in the
previous section has been corrected accordingly.

### What was ruled out for the HA_EC outlier

HA_EC (r = 0.432, every other tissue >= 0.837) got its own diagnosis rather than a story. Three
explanations were tested and all three failed:

1. **A uniform P-site shift.** Best-lag cross-correlation peaks at lag 0 (r = 0.4320), not at any
   nonzero lag, so the signal is blurred rather than displaced. Both HA_EC and HUVEC show a local
   bump at lag +/-3, which is codon periodicity and confirms the test is measuring what it should.
2. **A defective old pack.** Start-codon metagenes rebuilt from all 8 stored `target_counts.npy`
   (`scripts/pack_start_metagene.py`) peak at exactly +0 nt with frame0 0.732-0.825. HA_EC's old pack
   is mid-range (0.744, better than HUVEC's 0.732). The training data was not broken.
3. **Changed P-site calling.** Old-pack frame0 vs new metaplots frame0 agree within tissue for every
   tissue; HA_EC moves +0.028, indistinguishable from HUVEC's +0.027 (which has r = 0.957). New is
   uniformly slightly higher, as expected once the ncRNA/cross-gene filter is applied.

A fourth hypothesis -- that RiboCode `metaplots` per-sample offset DISAGREEMENT drives low r -- was
proposed, recorded as a prediction, and **falsified**: ES has offset disagreement (len30: 6 vs 12) and
returned the HIGHEST r of any tissue (0.980), while HA_EC has the SMALLEST disagreement (3 nt) and the
worst r. Read-length-set heterogeneity and depth were also checked and rejected.

So HA_EC's cause is open. The reportable finding -- reproducibility is tissue-heterogeneous and the
single-tissue number understates it -- does not depend on explaining it, and no causal claim about
HA_EC should appear in the manuscript until it is settled.

## The off-recipe non-canonical F1 was INFLATED by alignment artifacts (2026-08-13, RESOLVED)

> **This section replaces an earlier version titled "The final-recipe alignment costs NON-CANONICAL F1
> on mouse", which had the causation backwards.** That version reported the same per-class numbers
> but framed them as the model losing ground, and proposed a train/eval recipe mismatch. The check
> below falsifies that reading. The numbers are unchanged; the interpretation is inverted.

The deployed checkpoint is frozen and its inputs -- sequence, ORF track, RNA-seq coverage -- did not
move (RNA was never re-aligned). So its predictions cannot have changed, and they did not:

| file, matched cell `attn_ribo-gse243134_rna-gse243134` | off-recipe md5 | canonical md5 | |
|---|---|---|---|
| `pred_preddepth_collapsed.txt` | `27d5084b` | `27d5084b` | **byte-identical**, 19,601 calls both |
| `pred_obsdepth_collapsed.txt` | `556409b1` | `15c6ab49` | differs (scaled by OBSERVED depth, which moved) |
| `real_collapsed.txt` | `488e747a` | `88c4e3cf` | differs, 15,771 -> 14,877 |

The standalone arm is bit-for-bit unchanged, so its F1 drop (0.8762 -> 0.8661) is **entirely the
reference moving**. Nothing about the model got worse.

**What the reference lost is the whole story.** Class composition of the observed calls:

| class | off-recipe | canonical | change |
|---|--:|--:|--:|
| annotated | 10,535 | 10,385 | **-1.4%** |
| uORF | 2,587 | 2,263 | **-12.5%** |
| novel | 1,228 | 954 | **-22.3%** |
| dORF | 414 | 362 | -12.6% |
| Overlap_uORF | 591 | 534 | -9.6% |
| internal | 355 | 324 | -8.7% |
| Overlap_dORF | 61 | 55 | -9.8% |
| TOTAL | 15,771 | 14,877 | -5.7% |

Canonical alignment strips 22% of novel and 12.5% of uORF reference calls while leaving annotated
CDS almost untouched (-1.4%). That is the signature of soft-clipped alignments and unfiltered
ncRNA/cross-gene reads manufacturing spurious non-canonical ORFs -- precisely what `EndToEnd` and
`filter_tx_heldout` exist to remove. Deep, unambiguous annotated CDS is unaffected; marginal
low-count non-canonical calls are where artifacts land.

**So the off-recipe non-canonical F1 was inflated: the model was being credited for matching calls
that were alignment artifacts.** The per-class table below is unchanged from the earlier version, but
now reads as deflation of an inflated baseline rather than degradation of a model.

| class | ceiling off | ceiling canon | d ceiling | model off | model canon | d model |
|---|--:|--:|--:|--:|--:|--:|
| annotated | 0.9925 | 0.9915 | -0.0010 | 0.9884 | 0.9862 | -0.0021 |
| uORF | 0.7418 | 0.7278 | -0.0140 | 0.6108 | 0.5627 | -0.0481 |
| novel | 0.7237 | 0.6957 | -0.0280 | 0.6240 | 0.5770 | -0.0470 |
| internal | 0.5459 | 0.5179 | -0.0280 | 0.0232 | 0.0253 | +0.0021 |
| dORF | 0.5236 | 0.4892 | -0.0344 | 0.1326 | 0.0945 | -0.0381 |
| Overlap_uORF | 0.6888 | 0.6719 | -0.0169 | 0.4867 | 0.4610 | -0.0256 |
| Overlap_dORF | 0.4500 | 0.4638 | +0.0138 | 0.2018 | 0.1651 | -0.0367 |

The model falls further than the ceiling because the two Ribo-seq experiments lose their artifacts
independently and largely disagreed on them anyway, so experiment-vs-experiment agreement on the
surviving real calls is comparatively stable. The model, by contrast, was systematically predicting
some of what the artifacts produced -- unsurprising, since soft-clipping displaces P-sites
systematically rather than randomly, and a sequence-driven model can partially track a systematic
displacement.

**Consequence for the manuscript: quote the canonical non-canonical numbers, and do not compare them
to any previously reported non-canonical F1.** The earlier values were measured against a reference
containing ~22% more novel calls than the corrected pipeline produces. This also means the canonical
recipe's value is larger than the P-site totals suggested: a 4.9% drop in P-sites removed 22% of
novel calls.

The train/eval recipe-mismatch hypothesis from the previous version is **withdrawn** -- with
`pred_preddepth` byte-identical there is no sense in which the predictions were ever in or out of
distribution. The final-recipe retrain is still worth scoring, but for the separate question of whether
training on cleaner targets improves non-canonical calling, not to explain this drop.


## The model's experiment-unsupported ORF calls are false positives, not deep discoveries (2026-08-14, LANDED)

The mouse-liver UpSet (`figures/P9_upset_mouse3x3`) shows **1,942 ORFs called by mamba4 and by no
single Ribo-seq experiment**. That plot cannot say whether they are model false positives or real
ORFs no single experiment was deep enough to see, because every observed arm in the 3x3 is one
dataset deep. This resolves it.

**The arbiter.** `results/merged_liver_ribocode` is a JOINT RiboCode call over all 28 canonically
aligned mouse-liver libraries (Janich 5 + GSE243134 21 + Wang 2), submitted as
`scripts/merged_liver_ribocode.sbatch` (job 36780045, 2h01m). Not a union of three call sets: 28
per-library P-site offsets feed one periodicity test on the pooled signal, so an ORF with weak but
consistent signal everywhere can clear significance here while failing in each dataset alone.

| | raw ORF calls | vs merged |
|---|--:|--:|
| **merged, 28 libraries** | **26,388** | -- |
| GSE243134 (21 libs) | 14,877 | +77.4% |
| Janich (5 libs) | 14,322 | +84.2% |
| Wang (2 libs) | 13,060 | +102.1% |

**The control that makes it interpretable.** A recovery rate alone means nothing: the merged
reference is deeper, so it recovers more of everything. The informative comparison is against
OBSERVED calls that are equally unsupported by the rest of the 3x3 -- ORFs called by exactly one
experiment and missed by the other two. Those are known-real and equally lonely, so they measure how
often "unsupported by the other arms" just means "too deep for them".

Scored by `scripts/score_merged_liver_overcall.py`, standalone arm (`pred_preddepth`), GSE243134 RNA,
genomic keying, `build_loader` filtering, 22,974-tx space (merged reference = 15,254 calls in it):

| stratum | n | recovered by merged | rate |
|---|--:|--:|--:|
| model calls that ARE experiment-supported (upper anchor) | 11,441 | 11,362 | **99.3%** |
| **experiment-unique observed calls (CONTROL)** | 1,054 | 785 | **74.5%** |
| **model-unique, mamba4 (the question)** | **1,942** | **149** | **7.7%** |
| model-unique, attn | 2,056 | 176 | 8.6% |
| model-unique, called by BOTH models | 1,504 | 131 | 8.7% |

**A ~10x gap. The extras are false positives.** The anchor at 99.3% proves the merged reference is
not the limitation: it finds essentially every model call that any experiment saw. The control at
74.5% shows a genuinely real, singly-observed ORF gets corroborated three times in four. The model's
extras are corroborated at 7.7%.

**Two models agreeing does not rescue them.** The 1,504 calls made by BOTH mamba4 and attn but no
experiment recover at 8.7% -- indistinguishable from either model alone. Architectural agreement here
is a shared systematic error, not independent evidence. Do not use cross-architecture consensus as a
confidence filter.

**Splitting on whether the locus has any signal at all** closes the "no data there" objection:

| | mamba4 | attn |
|---|--:|--:|
| extras on a gene the merged reference DOES call (locus demonstrably translated) | 881 -> **16.9%** | 1,067 -> 16.5% |
| extras on a gene with NO merged call anywhere | **1,061** | 989 |

Both halves say false positive. Roughly **55% of the extras sit on genes with no detectable
translation in 28 pooled libraries** -- and since the universe is already restricted to salmon
TPM >= 1, those are RNA-expressed transcripts with no ribosome signal, which is exactly the lncRNA
over-call mode. For the other 45%, where the locus is demonstrably deep enough, the specific ORF is
still rejected 83% of the time against a 74.5% control. (The 0% recovery in the no-merged-call
stratum is definitional, not a measurement -- recovery requires a merged call on the gene. Read its
n, never its rate.)

Per class (`overcall_recovery_by_class.tsv`), the gap holds everywhere, including canonical CDS:

| ORF class | experiment-unique (control) | model-unique mamba4 |
|---|--:|--:|
| annotated | 92-100% | **8.5%** |
| novel (lncRNA) | 75-89% | **7.7%** |
| uORF | 66-84% | 11.5% |
| dORF | 74-83% | 3.2% |
| internal | 33-54% | 10.0% |

**FIGURE (added 2026-08-14): `figures/B9_overcall_validation`.** P9 carries only the headline annotation; the four strata, the cross-architecture bar and the gene-coverage split now have their own rigor panel.

**Consequences.**
1. The P9 UpSet's "model only" bar must be captioned as ~92% false positives, not as candidate
   discoveries. The 873 lncRNA calls shared by both models are ~91% false positives (74/873).
2. This gives the FP-filter (track 1) and the Poisson calibration (track 2) a real target metric:
   raise model-unique recovery toward the 74.5% control, rather than optimising F1 against a
   one-dataset reference that cannot see the difference.
3. `internal` is the one class where the control itself is weak (33-54%), so it is the one class
   where a low model number is genuinely ambiguous.
4. The merged reference is a fourth, deeper arm for validation. It does NOT replace the three
   independent observed arms the 3x3 factorial needs.

## Absolute total unique peptides, reported for every proteomics panel (2026-08-14)

Standing instruction: all proteomics work reports total unique peptides, including where the number
disfavours the model. Previously only discovery density, novel counts and canonical cost were shown,
and those three can all favour a database that nonetheless hands back FEWER peptides than searching
GENCODE alone. Four panels now carry the absolute total. Derivation everywhere is
`d_gencode_peptides + novel_peptides - gencode_arm.novel_peptides`, with the baseline-novel term
subtracted explicitly rather than assumed zero (verified identical to a direct
`(canonical+novel)_arm - (canonical+novel)_baseline` recomputation).

### P11 -- 12 mouse macrophage populations: the model is COST-NEUTRAL, not positive

**The BMDM-only version of this panel overstated the result and is superseded.** It claimed the
model's database was "the only one that increases TOTAL unique peptides". BMDM (+31) is the best of
the 12 populations on that axis.

| | model | null: every AUG |
|---|--:|--:|
| beats GENCODE-only on total unique peptides | **5 / 12** | **0 / 12** |
| median delta total | **-8** (of ~55,000, 0.01%) | **-1,784** |
| DB size vs null | 121-429x smaller (median 156x) | -- |
| discovery density vs null | 33-171x higher (median 111x) | -- |

Per population: BMDM +31, SmallIntestinal +20, SpleenRecruited +11, LargeIntestinal +5, Microglia +1,
LungRecruited -8, Peritoneal -8, LungResident -15, Kupffer -25, SpleenResident -26, LiverRecruited
-37, RAW264 -50. The null is negative in all twelve (-1,604 to -2,292).

**The defensible claim is cost-neutrality**: the model adds 10-37 novel peptides per population at a
median cost of 8 total peptides out of ~55,000, while the naive null costs ~1,800 every time.

### F2b, D12, D13 -- 5 human datasets: the Poisson arm is positive 9 times in 10

| arm | above zero on total unique peptides |
|---|--:|
| model (Poisson) | **9 / 10** (dataset x model) |
| model (theta=1) | 6 / 10 |
| naive AUG null | **2 / 5** (datasets) |

The single Poisson failure is **A549 / attn = -35** on the tryptic proteome, where the matching
mamba4 arm is +11 -- the sign flips on model choice, so that point stays in every caption. The
near-cognate null loses **7,780** peptides on A549, the largest loss measured anywhere in this
project.

D12's panel (d) supports "the mamba4 Poisson arm is the only arm positive on all four datasets"
(+11 / +14 / +10 / +21). That is **mamba4-specific**, and D12 defaults to `--model mamba4`; its
docstring now carries the qualifier. It is a claim about never LOSING, not about winning: `null_atg`
takes SU-DHL-4 outright (+47) and `model_standard` takes DoHH2 (+95).

### These two results must not be merged

Five human cell lines (9/10 positive) and 12 mouse macrophage populations (5/12, median -8) are
different populations answering the same question differently. The statement that survives both is
that the model's database is roughly cost-neutral to slightly positive on total peptides while the
naive nulls are reliably and substantially negative. Each figure's `FIGURE_DATA_INPUTS.md` carries an
explicit "do NOT merge with" note.

## A QC panel silently stopped showing any training data (2026-08-14, FIXED)

`figures/S_riboseq_qc` reported Ribo-seq library quality for "all 16 dataset arms (9 training + 5
held-out + 2 cell lines)". It was actually showing **7**, and **not one of them was a training
tissue**.

**Cause.** The training entries read the source study's own
`experiments/biotype_probe/expression_context_human/data/ribocode_per_tissue/<T>/<T>_pre_config.txt`.
That directory no longer exists -- the final-recipe alignment redo replaced it. The generator collects
unreadable sources into a `missing` list, prints it, and renders anyway. So the panel drew 5
held-out + 2 pgx arms, wrote its values JSON, exported its TSV, and passed every check, while its
README row and its own `FIGURE_DATA_INPUTS.md` continued to claim 16 arms including 9 training
tissues. A panel whose entire purpose is to establish that the TRAINING libraries clear a
periodicity floor had stopped showing a single training library, and nothing failed.

**Fix.** Repointed to `results/chothani_regeneration/_canonical_<T>/canon_<T>_pre_config.txt`, the
configs produced by `scripts/chothani_canonical_diverge.sbatch`. The panel now draws **15** arms:

| role | n | f0 range |
|---|--:|---|
| training (Chothani, canonical) | **8** | 74.6% (HUVEC) - 85.8% (ES) |
| held-out validation | 5 | 77.5% - 87.7% |
| proteogenomics cell line | 2 | 74.8% - 78.7% |

**8 training arms is correct, not 9.** Brain was dropped from the source study's 9-tissue panel for
low periodicity (period_obs 0.044) before any pack was built, so no canonical config exists. The
per-tissue library counts now reconcile exactly with P3: Fibroblast 32, VSMC 11, ES/Fat/HA_EC 6,
HCAEC/Hepatocytes 5, HUVEC 3 = 69 training, plus Hepatocytes' 5 held out.

**Blast-radius scan** (by artifact, not by the panel I happened to notice -- 15 files reference the
dead path):

| consumer | behaviour on the missing path | verdict |
|---|---|---|
| `figures/S_riboseq_qc/make_riboseq_qc.py` | collects to `missing`, renders anyway | **SILENT -- fixed** |
| `scripts/replicate_concordance.py` | `assert len(samples) == 32` on an empty glob | hard-fails, safe |
| `scripts/eval_localization.py` | builds a concrete file path, `open()` raises | hard-fails, safe |
| `scripts/plot_prediction_examples.py` | reads a concrete `*_collapsed.txt` | hard-fails, safe |
| 11 others (`build_psite_target.py`, `repack_union_*.sbatch`, ...) | historical, outputs already on disk | not live |

The distinguishing pattern is **glob-or-collect versus open-a-named-file**. Every consumer that
names a specific file failed loudly; the one that gathered a set and rendered what it found did not.
Printing `missing` is not a gate. This is the same lesson as `feedback_dump_silent_skip_onehot_fasta`
in a different place, and the generator now documents it in place.

### The BMDM total-peptide number flips sign on search parameters (2026-08-14)

Found while adding absolute totals to the tutorial. Same databases, same spectra, same BMDM sample;
only the MSFragger search configuration differs:

| BMDM run | model arm, delta TOTAL unique peptides |
|---|--:|
| frozen search parameters (used by P11 and F2b) | **+31** |
| unfrozen (quoted in results.md Task 53 and the tutorial) | **-21** |

The null arms do not flip: -1,604 frozen vs -1,998 unfrozen for `null_atg`, both large and negative.

**Consequence.** The model's effect on total unique peptides at BMDM is within noise of zero, and
quoting either number without stating the search configuration is quoting noise as a result. This
independently supports the 12-population conclusion: the model's database is COST-NEUTRAL rather
than positive. It is also a second reason the superseded BMDM-only P11 headline was wrong -- it was
not merely the best of 12 populations, it was the best of 12 under the more favourable of two search
configurations.

The tutorial's proteogenomics table now carries a delta-TOTAL column and states both numbers
alongside the 12-population result. Ratios (154x DB, 28x canonical cost, 123x density) are NOT
affected by this and remain comparable across runs; database size is a property of the database.

## P12: a PREDICTED ORF database vs a REAL Ribo-seq one, 12 macrophage populations (2026-08-14)

> **SUPERSEDED 2026-08-15 -- numbers below are the 7-aa mamba build.** Rebuilt at the correct 30-aa
> tryptic floor on the attention arm; see "Macrophage proteogenomics with the TRANSFORMER". The
> CONCLUSION is unchanged (a fixed real-Ribo-seq database still wins), but every number moved: the
> Ribo-seq DB is 506 not 699 seqs, novel 26 vs 24, delta-total +18 vs -3 (8/11 vs 5/11), density
> 2.4x not 2.9x. Retained as the record of what was believed on 08-14.

P11 asks whether the model beats naive enumeration. The harder and more useful question -- does a
predicted database substitute for one built from actual ribosome profiling? -- had complete data in
`pgx_xsubtype/crosssubtype_table.md` and **no figure**. Built as `figures/P12_macrophage_xsubtype`.

Three databases, identical spectra, 12 populations, **18 mzML fractions each** (so the sweep is
depth-balanced and cross-population comparison is not confounded by acquisition depth):

| arm | database | size |
|---|---|--:|
| `model_predicted` | model ORFs, called PER POPULATION | ~1,587 (varies) |
| `model_ribocode_bmdm_nt` | REAL Ribo-seq, RiboCode on BMDM untreated | **699, fixed** |
| `model_ribocode_bmdm` | REAL Ribo-seq, RiboCode on BMDM NT + LPS | 2,173, fixed |

**The asymmetry that defines it.** The Ribo-seq databases are BMDM-derived and SHARED by every row;
the model's is rebuilt per population. So BMDM is the MATCHED case for the Ribo-seq arms and the
other 11 are a TRANSFER test. Every statistic below excludes BMDM -- including it would flatter the
Ribo-seq arms on their own home ground.

### Transfer populations (n=11)

| arm | novel peptides | delta TOTAL peptides | above zero | discovery density | DB size |
|---|--:|--:|--:|--:|--:|
| model (per population) | 26 | **-8** | **4/11** | **13.1** | ~1,587 |
| real Ribo-seq, BMDM NT | 27 | **+4** | **7/11** | **38.6** | **699** |
| real Ribo-seq, BMDM NT+LPS | 31 | -2 | 4/11 | 14.3 | 2,173 |

(medians; BMDM matched, excluded: model 37 novel / +31 total, real NT 24 / +6, real NT+LPS 35 / -1)

**A fixed 699-sequence database derived from BMDM Ribo-seq matches or beats the per-population
predicted database on every axis** -- comparable novel yield, better on total unique peptides, and
**2.9x the discovery density** at less than half the size. This is a negative result for the model in
this application and is reported as one.

**What it does and does not license.** Do NOT write "the predicted database beats Ribo-seq"; on the
transfer populations it does not. The defensible claim is that the model **does not REQUIRE
Ribo-seq**: the predicted arm needs only RNA-seq plus sequence, while every Ribo-seq arm here needed
a ribosome profiling experiment in some macrophage. That is a claim about cost and applicability,
not accuracy.

**Two caveats that must travel with the panel.**
1. The delta-total differences are NOT a ranking -- they are the same size as the frozen-vs-unfrozen
   search-parameter swing measured on BMDM (+31 vs -21, identical databases and spectra). Discovery
   density is the robust axis, being a ratio insensitive to that swing.
2. The real-Ribo-seq arms are **ATG-only by construction** (RiboCode reports N-terminal extensions
   without testing whether the upstream start is used), so their ncStart is structurally zero. That
   is the one axis where the predicted arm has a structural advantage, and it is a property of the
   caller rather than evidence about non-AUG initiation.

**Audit note.** Only `P11_bmdm_proteomics` consumed macrophage data before this, and only the
model-vs-null arm from `frozen_reports`. Several other 12-population macrophage sweeps remain
table-only with no figure: `macro_model_vs_null_{mamba4,attn}_union{,_cal}.md` (model and calibration
variants) and `table_mamba4_union_cal.md`. They are model/calibration variants of P11's comparison,
not new questions, so they are lower value than P12 was.

## P13: the two ORF databases find LARGELY DIFFERENT spectra (2026-08-14)

> **SUPERSEDED 2026-08-15 -- numbers below are the 7-aa mamba build.** Rebuilt at 30 aa on the
> attention arm: 866 / 1,065 / 386 shared = 25.0% of a 1,545 union (was 885 / 1,012 / 367 = 24.0% of
> 1,530). The finding is UNUSUALLY well replicated -- two ORF-length floors and two architectures
> agree to within one percentage point, Jaccard median 0.247 -> 0.253, overlap 11/11 in both -- so
> complementarity is a property of predicted-vs-measured databases, not of one build.

P12 shows a predicted database and a real BMDM-NT Ribo-seq database deliver comparable novel-peptide
yield across 12 macrophage populations. That leaves the question counts cannot answer: are they
finding the SAME spectra? Two arms can each report 27 novel PSMs and either agree completely or not
at all. Built as `figures/P13_macrophage_psm_venn` from
`proteogenomics/scripts/macro_novel_psm_venn.py`.

**Method.** Novel PSMs at each arm's own 1% class-specific FDR, matched on **`(spec_id, peptide)`** --
a scan is shared only if both arms assigned it the same peptide. FDR is the pipeline's own
`pgx.report.cut`, imported not reimplemented; the single-pass loader is asserted equal to
`pgx.score_compare.passing_psms` on BMDM (107 PSMs, cut 20.68), and **all 12 x 2 arm counts match
`novel_psms` in the frozen reports exactly**.

### Transfer populations (n=11, BMDM excluded as the Ribo-seq arm's home ground)

| | PSMs |
|---|--:|
| model-predicted database | 885 |
| real Ribo-seq database (BMDM NT, 699 seqs, fixed) | 1,012 |
| **shared** | **367** |
| model-only | **517** |
| Ribo-seq-only | **644** |
| union | 1,530 |

**Shared = 24.0% of the union.** Every population overlaps somewhat (11/11), but Jaccard never
exceeds 0.32 and reaches 0.057 at Microglia.

| population | model | Ribo | shared | J |
|---|--:|--:|--:|--:|
| SmallIntestinal | 156 | 106 | 63 | 0.317 |
| Peritoneal | 106 | 107 | 46 | 0.275 |
| LargeIntestinal | 62 | 112 | 37 | 0.270 |
| LungRecruited | 138 | 105 | 51 | 0.266 |
| SpleenRecruited | 88 | 105 | 39 | 0.253 |
| SpleenResident | 90 | 97 | 37 | 0.247 |
| RAW264 | 46 | 144 | 34 | 0.218 |
| LiverRecruited | 73 | 56 | 23 | 0.217 |
| LungResident | 58 | 102 | 24 | 0.176 |
| Kupffer | 54 | 36 | 10 | 0.125 |
| Microglia | 14 | 42 | 3 | 0.057 |
| BMDM (matched) | 107 | 61 | 28 | 0.200 |

### The two databases are COMPLEMENTARY, not redundant

The model finds **517 novel PSMs a real BMDM Ribo-seq database misses**; that database finds **644
the model misses**. Neither is a subset of the other in any population. This is the sharpest
available statement about what the predicted database adds: not more discoveries than Ribo-seq (P12
shows it does not), but DIFFERENT ones. The practical implication -- that the union beats either
alone -- is untested and would need its own search.

### `diff-pep` = 1 across all 12 populations, which is a QC result

Scans passing in both arms but assigned different peptides total **one** across the whole sweep. When
both databases pass a spectrum they essentially always agree on the peptide, so this figure is
measuring which spectra clear FDR, not rank-1 churn in peptide assignment. Had that number been
large, the Venn would have been uninterpretable -- the intersection would depend on how the
assignment was tie-broken.

**Caveat that must travel with it.** Each arm carries its own class-specific cut, so a spectrum can
be inside one circle and outside the other purely by threshold. That is a real difference in what a
database delivers at fixed error rate, not an artifact, but it does mean the non-overlap is partly a
threshold effect and not purely a sequence-space effect.

## Canonical retrain, mamba4: a wash on profile, and that does NOT settle the question (2026-08-14, PRELIMINARY)

The canonical mamba4 retrain (job 36769695) finished: 14h37m, 28 epochs, early-stopped from best at
epoch 21, val Pearson 0.5998. It exists to answer whether training on CLEANER targets -- the
final-recipe alignment recipe, which strips the artifacts that inflated the old non-canonical F1 --
improves non-canonical calling.

On held-out Hepatocytes profile metrics it is indistinguishable from the deployed model:

| | final-recipe retrain | deployed (union mm1) | delta |
|---|--:|--:|--:|
| Pearson, all | 0.6833 | 0.6851 | -0.002 |
| Pearson, protein-coding | 0.6978 | 0.6979 | -0.000 |
| Pearson, lncRNA | 0.4514 | 0.4793 | **-0.028** |
| period_pred | 0.2691 | 0.2632 | +0.006 |
| frame0_pred | 0.8229 | 0.8196 | +0.003 |

**This is reported as preliminary and must not be quoted as the answer, for three reasons.**

1. **The test sets are not matched.** Canonical n=70,564 vs deployed 70,883; lncRNA 1,005 vs 1,049,
   a 4.2% difference. The canonical ncRNA and cross-gene filters drop reads, so some transcripts fall
   below the scorable threshold. The lncRNA delta could be composition rather than learning.
2. **The reference moved as well as the model.** Canonical lncRNA `period_obs` is 0.0266 vs 0.0241
   for the deployed test target. Cleaner targets are the entire point of the final recipe, so the
   two numbers are not a like-for-like grade -- the same trap as the earlier "the model lost ground"
   framing, which turned out to be the ceiling moving.
3. **Profile Pearson is not how models are judged here** (`feedback_orf_call_eval_standing`). The
   evaluation that settles this is the RiboCode drop-in, both arms, scored per ORF class against a
   final-recipe reference with n_ref stated. It has not run.

The canonical runs also write no `pertx.tsv`, so even a matched-transcript profile comparison needs
an eval re-run rather than a re-read of what is on disk.

**Blocked on GPU contention**: job 36783778 (the CAR-T canonical dump, which P5 needs) is queued
ahead on the gpu partition, and the two attn canonical runs are still training.

## The 30-AA ORF floor: followed on ORF calls, NOT in proteogenomics -- and in HUMAN data it matters a lot (2026-08-14, RESOLVED)

**The rule.** Minimum ORF length 90 nt. RiboCode's `ORF_length` EXCLUDES the stop codon (verified:
1911 nt / 3 = 637 aa exactly for a 637-aa AAseq), so 90 nt is **exactly 30 amino acids**.

**ORF-call track: FOLLOWED, and correctly.** Two-stage by design -- call permissively
(`ribocode_dropin.py --min_aa 5`), then filter at SCORING, so the floor is applied through one code
path to the model's calls and the real Ribo-seq calls alike. `compare_dropin_calls.build_loader`
defaults to `min_len=90` and passes it regardless of `is_pred` (that flag gates only the enrichment
filter). `compare_dropin_ctg` and `compare_dropin_allalt` default to 90 too. **No scorer overrides
it.** Every ORF-call number here -- the 3x3, P5, P9, A4, B9, the merged-liver over-call work -- is
30-aa filtered symmetrically. Filtering at scoring rather than calling is deliberate: it guarantees
model and reference are cut at the same place, which a call-time threshold cannot.

**Proteogenomics track: NOT followed** -- `pgx/build_dbs.py` and `pgx/coding_potential.py` default to
`--min-aa 7` ("detectability"). The two tracks report on different ORF populations and nothing
declared it.

### The measurement, and why one dataset family was not enough

`proteogenomics/scripts/orf_length_audit.py` counts novel peptides that PASSED 1% class-specific FDR
and would be LOST at a 30-aa floor -- peptides whose supporting ORFs are **all** under 30 aa. A
peptide with any >= 30 aa supporter survives, so counting peptides that merely touch a short ORF
would overstate the cost. FDR from `pgx.report.novel_at_class_fdr`, imported not reimplemented.
`unmapped = 0` everywhere, so every peptide resolved to a supporting ORF.

| arm | MOUSE macrophage (12 pops) | HUMAN (5 datasets) |
|---|--:|--:|
| model (Poisson) | **0 / 295 (0.0%)** | **14 / 52 (26.9%)** |
| model (theta=1) | -- | **25 / 76 (32.9%)** |
| real Ribo-seq BMDM NT | 8 / 306 (2.6%) | -- |
| null_atg | 19 / 453 (4.2%) | 22 / 72 (30.6%) |
| null_nc | -- | 2 / 38 (5.3%) |
| **CPAT** | -- | **0 / 39 (0.0%)** |
| **CPC2** | -- | **0 / 55 (0.0%)** |

**The two families give opposite answers, and the mouse one alone would have been badly
misleading.** On mouse macrophages (tryptic) not one model discovery depends on a short ORF. On the
human datasets (mostly HLA-I) roughly a THIRD do.

**The comparison this actually threatens is model vs CPAT/CPC2, not model vs null.** Against the null
the effect is near-symmetric (26.9% vs 30.6%), so the 7-aa floor does not flatter the model there.
But CPAT and CPC2 have **zero** short-only discoveries -- their databases are 0-2% under 30 aa,
because composition-based coding-potential scoring selects long, ORF-like sequences by construction.
A 30-aa floor therefore cuts the model by ~30% and leaves the competitor untouched.

### D12's headline does not survive a 30-AA floor intact

D12 claims "the model wins all three HLA-I immunopeptidomes". Applying the floor to both sides:

| dataset | best model now | best CPAT/CPC2 | winner | best model @30aa | CPAT/CPC2 @30aa | winner |
|---|--:|--:|---|--:|--:|---|
| A549 (tryptic) | 11 | 25 | CPAT/CPC2 | 10 | 25 | CPAT/CPC2 |
| **HBL-1** | 14 | 10 | model | **8** | **10** | **CPAT/CPC2 -- FLIPS** |
| SU-DHL-4 | 22 | 14 | model | 16 | 14 | model |
| DoHH2 | 32 | 8 | model | 15 | 8 | model |

**Under a 30-aa rule the claim becomes "wins 2 of 3", not 3 of 3.**

### Decision: keep 7 aa, and state the dependence

Keeping it is still right -- HLA-I peptides are 8-11 aa, and microproteins under 30 aa are precisely
what a translation model finds that a composition-based selector rejects. That is a real capability,
not an artifact. But it is now a **declared, quantified dependence** rather than an undeclared one:
the model's HLA-I advantage over CPAT/CPC2 rests substantially on ORFs the ORF-call track excludes by
rule, and D12's three-of-three claim is contingent on the 7-aa setting.

Declared in `docs/PIPELINE_POLICY.md`; per-panel notes in the six proteomics
`FIGURE_DATA_INPUTS.md`. Per-dataset numbers in
`proteogenomics/data/orf_length_audit_human.json` and `.../pgx_xsubtype/orf_length_audit.json`.

**A CORRECTION TO AN EARLIER ENTRY IN THIS FILE.** The first version of this section concluded "it
costs zero discoveries" and recommended keeping 7 aa on that basis. That was measured on mouse
macrophages ONLY and does not hold for the human datasets, where the cost is ~30% and one D12
comparison flips. The conclusion (keep 7 aa) survives; the reason given for it did not.

**A reporting bug found and fixed inside the audit.** The first run reported `DB<30aa` as 2.0% where
the true novel-target fraction is 22.0%, because `fasta_lengths` returned every entry in the search
fasta -- 117,518 of them for 1,533 novel ORFs, the rest decoys and the full canonical proteome. The
short-only counts were never affected (they resolve specific `nuORF|` accessions), but the two
columns contradicted each other in one table. Fixed to novel targets only; the re-run reproduces
short-only exactly (0 / 8 / 19) and the corrected column matches an independent manual count.

## Macrophage proteogenomics with the TRANSFORMER, at the correct 30-AA floor (2026-08-15, LANDED)

Two changes at once, deliberately: the model (mamba4 -> attn) and the ORF floor (7 -> 30 aa, Rule 5,
tryptic). Running attn at 7 aa would have violated the standing rule; comparing attn@30 against the
existing mamba@7 would have confounded two variables. So BOTH models were searched at 30 aa, all
7 arms per population against **identical spectra** under frozen parameters, one pass per population.

Pipeline: attn predictions already existed (`pred_attn_union`, 12 populations) so no GPU was needed.
Calls (job 36789585) -> 30-aa databases (36789634) -> 7-arm search (36789647) -> reports.

### The transformer is the better database

| median across 12 populations | mamba4 | **attn** | real Ribo-seq | null |
|---|--:|--:|--:|--:|
| DB size | 1,252 | **1,062** | 506 | 84,226 |
| novel peptides | 27 | **26** | 26 | 43 |
| **discovery density (/1k)** | 17.9 | **21.5** | 51.4 | 0.5 |
| delta TOTAL peptides | -8 | **+3** | +16 | -1,130 |
| above baseline | 5/12 | **6/12** | 9/12 | 0/12 |

**attn finds the same peptides from a 15% smaller database**, so density rises 20% (17.9 -> 21.5).
Its ~15-20% fewer ORF calls (consistent across all 12 populations, ratio 0.77-0.91) were
disproportionately the ones that yield no confident peptide -- consistent with B9's finding that
~92% of experiment-unsupported calls are false positives. The conservatism costs nothing.

That settles the model choice on merit, not familiarity: attn is the smaller, denser, marginally
better database AND the CPU-runnable one (mamba_ssm requires CUDA).

**Do not rank the models on delta-total** (+3 vs -8). That gap is inside the +-50 frozen-vs-unfrozen
search-parameter swing. The DENSITY difference is the robust one -- a ratio, insensitive to that
swing, moving the same direction in 8 of 12 populations.

### CORRECTION: two P11 headline ratios were inflated at 7 aa

| | 7 aa (reported 2026-08-14) | 30 aa (correct) |
|---|--:|--:|
| DB size vs null | **156x smaller** | **78x smaller** |
| discovery density vs null | **111x higher** | **43x higher** |

At 7 aa the null carried 235,686 sequences, ~65% of them sub-30-aa ORFs **the ORF-call track excludes
by rule**. The correct floor shrinks the null to ~84,000 and the ratios halve. The earlier numbers
were inflated by counting sequences that should never have been in the comparison. Both P11's
docstring and its FIGURE_DATA_INPUTS.md now carry the corrected values and the reason.

The model-vs-null advantage is therefore about HALF what the 7-aa comparison implied -- still large,
still the right direction, but the honest figure is 78x/43x.

### P13 (novel-PSM overlap) holds for the transformer

Transfer totals (n=11, BMDM excluded): attn 866 PSMs, real Ribo-seq 1,065, **shared 386** = 25% of a
1,545 union. Still complementary rather than redundant: 479 attn-only vs 678 Ribo-only PSMs.
`same_scan_diff_peptide` remains 1 across all 12, so the disagreement is about which spectra clear
FDR, not about peptide assignment.

### P12 rebuilt: correcting the floor did NOT rescue the model here

Worth stating precisely because the 30-aa rebuild helped the model in P11. Against the real
Ribo-seq databases the transformer at 30 aa does slightly WORSE than mamba at 7 aa did:

| transfer medians (n=11) | model 7-aa mamba | **model 30-aa attn** | real Ribo-seq NT |
|---|--:|--:|--:|
| novel peptides | 26 | **24** | 27 -> **26** |
| delta TOTAL peptides | -8 | **-3** | +4 -> **+18** |
| above baseline | 4/11 | **5/11** | 7/11 -> **8/11** |
| discovery density (/1k) | 13.1 | **21.3** | 38.6 -> **51.4** |
| density gap | 2.9x | **2.4x** | -- |

The floor change lifted BOTH arms, and it lifted the Ribo-seq arm more on total peptides. The
density gap narrowed (2.9x -> 2.4x) while the total-peptide gap widened (+4/-8 -> +18/-3). P12
remains the honest counterweight to P11 and the conclusion is unchanged: **do not caption this as
the model beating Ribo-seq.**

### Two staleness bugs found and fixed while rebuilding these panels

1. **P13 was never actually regenerated on 08-15.** Its underlying `novel_psm_venn.json` was
   recomputed on the attention tree, but the FIGURE still read the 7-aa mamba JSON, and the status
   doc recorded it as done. Two panels were stale, not one. The lesson is the artifact's own
   recorded source (`*_values.json` -> `source`) is the only evidence a panel was rebuilt; a note
   saying it was rebuilt is not.
2. **`make_poster_layout.py` would have crashed**, not drifted: it indexed `transfer_summary` with
   the literal `"model_predicted"`, which no longer exists. It now reads the arm name from
   `P12_values.json["model_arm"]`. Its "MUST STATE" prose also hardcoded 699 seqs / 2.9x / "5 of 12"
   even though the file's own docstring claims it cannot drift, because only the headline numbers
   were pulled live and the guidance text was not. Both now derive from the JSON.

Databases sizes are now measured at plot time in P12 and P13 rather than written into labels, since
the same two Ribo-seq databases are 699/2,173 seqs at 7 aa and 506/1,265 at 30 aa.

## Brain and GSE39561 fail DIFFERENTLY, and low mapping does not predict which (2026-08-15)

Both are the project's "bad" Ribo-seq datasets and both have ~30% unique mapping. Under the SAME
final recipe they behave oppositely:

| | unique mapping | read lengths metaplots selected | verdict |
|---|--:|--:|---|
| GSE39561 (THP-1) | 29.5-32.7% | **0 of 3 libraries** | no periodicity -- UNUSABLE |
| Brain (Chothani) | 29.55% | **5 of 5 libraries** | periodic but weak -- usable, excluded on quality |

**Low unique mapping does not predict periodicity failure.** Brain selects sensible read lengths
(25-30 nt) at the canonical 12-nt offset in every library; GSE39561 selects none in any. So the
source study's exclusion of Brain for `period_obs 0.044` was a judgement about SIGNAL STRENGTH, not
about a broken library -- a materially different justification from the one assumed here previously.

**GSE39561 is now a VERIFIED negative control.** Its final-recipe rebuild also selected zero read
lengths, so `--alignEndsType EndToEnd` cannot rescue it: the recipe fixes P-site DISPLACEMENT but
cannot manufacture periodicity that is not in the library. Its final-recipe pack is coverage-only
(39,611 tx, ZERO P-sites, n_ribo=0 n_rna=4, RNA borrowed from GSE208041).

Brain was never fetched before this: 5 runs (SRR15513148-152) downloaded from GSE182371, md5-verified
against ENA (5/5), aligned with the final recipe (filter 9.08%, matching the 8.86% median across the other 71
Chothani runs). Metaplots at `results/chothani_regeneration/_canonical_Brain/`, which is the path
`figures/S_riboseq_qc` globs -- that panel becomes complete at 9 tissues on its next run.

## Provenance and identifier fixes (2026-08-15)

**Dataset registry** at `data/dataset_registry.tsv` + `scripts/dataset_labels.py`. Labels had grown as
a mix of surnames, GEO accessions, cell types and product names -- and the surname-labelled datasets
were surname-labelled because NO accession was recorded anywhere. Scheme is
`<species>_<tissue>_<accession>`: the accession is the unambiguous key, species+tissue is what a human
scans for, and the pair resolves the two-THP-1 collision (GSE208041 and GSE39561 are both THP-1) that
neither component resolves alone. Internal directory keys are deliberately NOT renamed. Accessions
supplied 2026-08-15: Janich **GSE67305** (was recorded nowhere), Wang **GSE94982** and Chothani
**GSE182371** (both already in methods.md, never used as labels).

**Chothani samplesheet was under-reporting the study.** It listed 71 runs / 7 tissues against 74
aligned BAMs. The 3 extras were **HUVEC** (SRR15513213/214/215) -- aligned, in the training set, in
P3 and S_riboseq_qc, but absent from the authoritative input record, so a rebuild from the sheet would
have silently dropped the shallowest training tissue. Now 79 runs / 9 tissues (HUVEC + Brain),
matching ENA's count for PRJNA756018 exactly.

**THP-1 GSE208041 is Vehicle-only, verified.** The study has 3 conditions x 5 reps (Vehicle, LPS,
LPSDex). Our 5 RPF runs are exactly the 5 Vehicle reps and our 4 RNA runs exactly the 4 Vehicle reps;
no treated sample is in the pack. The provenance doc recorded arm (RPF/RNA) but NEVER the treatment
condition, so nothing in the repo stated this -- the same class of gap as the missing accessions.

**CAR-T is now canonical and the ceiling moved, as predicted.** F1 0.904 -> 0.887 with n_ref
8,856 -> 8,253: the reference lost 603 calls when ~29% PCR duplicates came out. The model's
predictions barely changed; what moved is what it is graded against.

**P5 silently dropped an arm and was fixed.** Repointing to the canonical scored table made it a
3-arm panel -- THP-1 vanished, because the logic switched WHOLE FILES and the final-recipe rescore only
contained CAR-T. It now prefers canonical PER DATASET, falls back to the original, prints the
provenance of each, and records `source_per_dataset` in the values JSON. This is the third instance
this week of the same failure mode: a change that makes something more correct while silently
reducing coverage (see also the QC panel losing its training tissues).

## Task 86 SETTLED: final-recipe training does NOT change ORF calling (2026-08-15, LANDED)

Three final-recipe retrains had existed since 2026-08-13 with no ORF-call evaluation -- profile
metrics were a wash, and profile Pearson is explicitly not the criterion
(`feedback_orf_call_eval_standing`). Jobs 36790643_0 + 36790648_[1-4].

**Controlled by evaluating every arm on ONE substrate.** All five arms were dumped with
`RIBO_PACK_SUFFIX=canon`, so each scores against the same final-recipe Hepatocytes reference and
consumes the same RNA coverage; only the TRAINING substrate differs. Scoring the final-recipe model
against the final-recipe reference and the deployed model against the mm1 reference would have
compared two different reference sets, and any difference could then be reference quality rather
than model quality.

Control: `real` depends only on the pack, so all five arms must call the same ORFs. Verified --
ORF_ID sets and `ORF_type` identical across all five.

### Result: no material difference, both architectures

Standalone arm (`pred_preddepth`), shared reference, n_ref 3,028 non-canonical:

| run | overall F1 | annotated | non-canonical F1 | nc recall | nc precision |
|---|--:|--:|--:|--:|--:|
| final-recipe/attn | 0.9104 | 0.9894 | 0.5776 | 0.6050 | 0.5525 |
| final-recipe/attn_nof012 | 0.9110 | 0.9894 | 0.5813 | 0.6106 | 0.5546 |
| final-recipe/mamba4 | 0.9141 | 0.9892 | 0.5911 | 0.6073 | 0.5758 |
| deployed/attn | 0.9106 | 0.9893 | 0.5741 | 0.5928 | 0.5566 |
| deployed/mamba4 | 0.9142 | 0.9891 | 0.5914 | 0.6047 | 0.5787 |

**attn +0.0034, mamba4 -0.0003** on non-canonical F1 -- both inside the +-0.005 materiality band set
before the runs. **The deployed checkpoints stay deployed**; no `loto_canon` checkpoint goes into any
figure. The ORF-call criterion agrees with the profile metrics rather than overturning them.

Predicted by the reference sizes: the final-recipe Hepatocytes reference is 18,813 calls vs 18,964
on mm1, only **-0.8%**, so there was little headroom for cleaner targets to matter. The training
tissues were already on-recipe when their packs were built, so the rebuild changed almost nothing for
them; the held-out human datasets were not, which is why THP-1 (-7.5%) and CAR-T (-6.8%) moved and
Hepatocytes did not.

## The held-out F1 drop is the REFERENCE moving, not the predictions (2026-08-15)

THP-1 fell 0.851 -> 0.831 and CAR-T 0.904 -> 0.886 when their packs were rebuilt, which reads as
"standardising the pipeline made the model worse". It did not. A 2x2 cross-pairing of
{old pred, new pred} x {old reference, new reference} on THP-1 attn, all scored on one transcript
space, decomposes it:

| | ref = OLD (mm1) | ref = NEW (final recipe) |
|---|--:|--:|
| pred = old | 0.4925 | 0.4186 |
| pred = new | 0.4925 | 0.4185 |

**REFERENCE effect -0.0740. PREDICTION effect -0.0001.** The predictions are unchanged (15,247 vs
15,249 calls). The entire drop is the reference losing 1,079 calls (14,304 -> 13,225).

Mechanically consistent: rebuilding THP-1 moved the Ribo-seq target a lot (P-sites -8.6%) but the RNA
coverage that drives the predictions by **+0.005%**. Same inputs, same predictions, shorter ruler.

**WHAT REMOVED THOSE CALLS IS NOT THE SAME FOR THE TWO DATASETS.** An earlier version of this section
attributed both to PCR duplicates. That is wrong for THP-1 and the distinction matters:

| | UMIs? | what the rebuild removed |
|---|---|---|
| **THP-1 GSE208041** | **NO** -- the submitter stripped them before deposit | the ncRNA + cross-gene filter ALONE (~9.7% of records; measured 8.19% cross-gene + 1.45% ncRNA on this dataset) |
| **CAR-T GSE304796** | yes, 12 nt 5' UMI | ~29% PCR duplicates **plus** the ~9.7% filter. `umi_tools dedup` ran on the GENOME bam while the pack was built from the undeduplicated transcriptome bam |

Per `docs/PIPELINE_POLICY.md`: "GSE304796 has a 12 nt 5' UMI; nothing else does", and "GSE208041 and
GSE39561 have no UMIs at all ... for those two the filter is the only gap."

**The implication is not flattering either way.** The reference calls the final recipe removed were
ones the model HAD been matching, so removing them converts those predictions from TP to FP. For
THP-1 the model was being credited for calls supported by **cross-gene-ambiguous and ncRNA reads**;
for CAR-T, by those *and* by PCR duplicates. Both old numbers were inflated, for partly different
reasons.

## Recipe uniformity audit: two datasets are OFF-RECIPE (2026-08-15)

Prompted by the question "was the observed data all computed the same way?". Evidence is a recorded
`filter_tx_heldout` run per dataset, pack build time vs the 2026-08-12 adoption, and provenance
completeness -- not which scripts mention the filter.

| dataset | built | filter runs | verdict |
|---|---|--:|---|
| Chothani training (8 tissues) | 08-13 | 69 | on-recipe |
| Hepatocytes holdout | 08-13 | 5 | on-recipe |
| THP-1 GSE208041 | 08-15 | 5 | on-recipe |
| CAR-T GSE304796 | 08-14 | 3 | on-recipe |
| GSE39561 | 08-15 | 3 | on-recipe |
| **iPSC-CM Ruiz-Orera** | **07-14** | **0** | **OFF-RECIPE** |
| B721.221 | 08-06 | 0 in logs, but **7/7 BAMs carry filter provenance** | **ON-RECIPE** (audit was wrong) |

Cross-gene drop is a consistent **5.6-9.8%** across the 19 on-recipe dataset directories.

**CORRECTION, same day: B721 is ON-RECIPE and the audit was wrong about it.**
`proteogenomics/scripts/align_b721_ribo.sbatch` uses EndToEnd + mm1 and calls the filter, and all 7
BAMs carry its `@PG` provenance. The audit scraped `logs/*.out`, and B721's filter ran from a script
whose logs live elsewhere, so it scored zero. **Absence from a log scrape is not evidence of a
skipped step.** The converse also failed: janich and gse243134 show no filter evidence in `@PG`
because their filtered copies live in `data/liver3x3/pool_*/` while `heldout_bam/` keeps the raw
alignments. Neither method alone is authoritative -- the question is which artifact fed which pack.

That leaves **iPSC-CM alone**, and its status is UNKNOWN rather than confirmed off-recipe: its BAMs
and source hd5 are deleted and no align script survives, so there is nothing left to check. Its pack
predates the 2026-08-12 adoption, which is suggestive but not proof. It was re-downloaded from ENA
on 2026-08-15 (15 files, 34.4 GB, all md5-verified) so it can be rebuilt on the final recipe and the
question becomes moot.

The one 0.00% cross-gene entry is `_filtered_HUVEC` at 0.0022%, ~4,000x below the production range --
the coordinate-sorted-input failure `filter_tx_heldout` now hard-guards against, and the run that
produced the wrong "the filter is a no-op" conclusion on 2026-08-12. Nothing current reads it.

**The 30-aa floor passes everywhere**: 13.7% of raw calls are under 90 nt and all are removed at
scoring; 0 survive; no consumer overrides `min_len=90`.

## P5 is 4/4 ON-RECIPE, and the "off-recipe pair" was mostly a provenance gap (2026-08-16, LANDED)

Task #92. Started from a recipe audit that flagged iPSC-CM and B721 as off-recipe. **Both flags were
wrong or overstated, and only measurement settled it.**

### B721: never off-recipe. The audit method was.

`proteogenomics/scripts/align_b721_ribo.sbatch` uses `--alignEndsType EndToEnd`,
`--outFilterMultimapNmax 1` and calls `filter_tx_heldout.py`; all **7/7** BAMs carry that filter's
`@PG` provenance. The audit scraped `logs/*.out`, and B721's filter ran from a script whose logs live
elsewhere, so it scored zero. **Absence from a log scrape is not evidence a step was skipped.**

The converse also failed: janich and gse243134 show NO filter evidence in `@PG` because their
filtered copies live in `data/liver3x3/pool_*/` while `heldout_bam/` keeps the raw alignments.
Neither logs nor headers are authoritative alone -- the question is which artifact fed which pack.
A wrong off-recipe caveat had already been written onto A4, the project's only NON-null benchmark.

### iPSC-CM: rebuilt end-to-end, and it barely moved

Its BAMs, source hd5 and align script were all gone, so nothing survived to inspect -- status
UNKNOWN, not confirmed bad. Settled by rebuilding rather than inferring: FASTQs re-fetched from ENA
PRJEB65856 (15 files, 34.4 GB, **all md5-verified**, serial single-connection on the head node), 5
Ribo runs realigned on the final recipe (job 36791179, 5/5), pack rebuilt (36791252), dumped both
architectures (36791660).

Cross-gene drop on the new alignments is **7.8-8.4%**, squarely in the on-recipe band.

| | P-sites | n_ref | attn F1 | mamba4 F1 |
|---|--:|--:|--:|--:|
| OLD | 324,003,052 | 13,147 | 0.8763 | 0.8799 |
| NEW (final recipe) | 315,879,911 (**-2.5%**) | 13,087 (**-0.5%**) | 0.8751 | 0.8794 |

**-0.001 F1.** The old pack HAD been filtered after all; the -2.5% is the realignment alone.

### The correction is NOT uniform, and the blanket estimate was wrong for 3 of 4 arms

The poster session proposed printing "~0.02, likely overstated" for the off-recipe arms, derived from
THP-1 and CAR-T. Measured:

| arm | n_ref change | F1 change | why |
|---|--:|--:|---|
| THP-1 | -7.5% | **-0.020** | pack skipped the ncRNA/cross-gene filter entirely |
| CAR-T | -6.8% | **-0.018** | same, plus ~29% PCR duplicates |
| hepatocyte | -0.6% | **+0.002** | packed_union was ALREADY filtered; only realignment differs |
| iPSC-CM | -0.5% | **-0.001** | same |

**The ~0.02 transfers only to packs that skipped the filter.** Two did; two did not, and one of those
moved UP. Applying one blanket adjustment would have been wrong in both direction and magnitude for
half the panel.

### Final P5, attn, standalone, all four arms on one substrate

| arm | F1 | annotated | non-canonical | n_ref |
|---|--:|--:|--:|--:|
| hepatocyte (LOTO) | 0.911 | 0.997 | 0.594 | 16,236 |
| iPSC-CM (Ruiz-Orera) | 0.875 | 0.997 | 0.550 | 13,087 |
| THP-1 (GSE208041) | 0.831 | 0.997 | 0.525 | 13,225 |
| CAR-T (GSE304796) | 0.886 | 0.994 | 0.580 | 8,253 |

`mixed_substrate: false`. P5 no longer warns, and the four arms are mutually comparable for the first
time.

### Tooling added

`prepare_pack.py` had `--reuse-target` (new coverage, keep Ribo) but not its mirror. Added
**`--reuse-coverage`** (keep universe + coverage, write a new Ribo target) in both `prepare_pack` and
`prepare_from_bams` -- the "Ribo re-aligned, RNA untouched" path the final recipe requires, per
PIPELINE_POLICY's "RNA-seq stays Local" and "the universe is NOT redefined".

The rebuild asserts tx_order, lengths AND coverage byte-identical to the source before copying the
ORF track, because a misaligned track yields plausible numbers rather than a crash.

**Scope note:** 18.8 GB of iPSC RNA FASTQs were downloaded and NOT needed -- policy keeps RNA out of
the recipe change. Checking that first would have halved the download.

## Task 75 LANDED: depth explains about HALF the between-experiment disagreement (2026-08-16)

The matched-depth arm was deferred with an explicit criterion -- run it if the depth ordering
reappears at natural depth. **It reappeared unambiguously**, so this became necessary.

Natural depths: wang 78.3M < janich 110.0M < gse243134 150.5M. The observed-vs-observed matrix
tracked that ordering rather than lab or protocol -- wang scored **0.612/0.575** (novel recall) as the
COMPARED set but **0.805/0.834** as the REFERENCE. A shallow library misses a deep reference's calls
while a deep experiment recovers most of its few. That matters because **E2 uses this matrix as
"the between-experiment ceiling"**, so the ceiling confounded protocol with sequencing depth.

### Design

All three thinned to **75M**, INCLUDING wang (f = 0.4982 / 0.6818 / 0.9584). Thinning only the two
deeper sets would have left wang as the single arm carrying no thinning noise -- a subtler version of
the same asymmetry. **Three seeds**, because binomial thinning is stochastic; the seed spread is the
error bar. Jobs 36791789_0-8, 9/9 COMPLETED. `ribocode_dropin.py --subsample`, the depth_crossover
machinery, so no new code.

### Call counts converge, but not completely

| dataset | natural | matched 75M (3 seeds) |
|---|--:|--:|
| gse243134 | 14,877 | 13,688-13,717 |
| janich | 14,322 | 13,701-13,752 |
| wang | 13,060 | 13,035-13,050 |

Spread falls from **14%** to **5%**. gse243134 and janich become indistinguishable; wang stays ~700
calls (5%) below both. **That residual is the genuine protocol difference**, and it is much smaller
than the natural-depth gap implied.

### The asymmetry roughly halves

novel recall, wang as COMPARED vs wang as REFERENCE:

| | natural | matched |
|---|--:|--:|
| wang as compared (mean) | 0.594 | 0.656 |
| wang as reference (mean) | 0.820 | 0.773 |
| **asymmetry gap** | **0.226** | **0.117** |

uORF behaves the same way: gap 0.335 -> 0.205. Seed sd is 0.002-0.017, far below the effect.

### E2's ceiling band NARROWS, and that is the consequence for the figure

Min-max F1 over the 6 ordered pairs, which is exactly what E2 shades:

| class | natural | matched 75M | width |
|---|--:|--:|--:|
| annotated | 0.991-0.995 | 0.990-0.993 | 0.004 -> 0.003 |
| uORF | 0.682-0.832 | 0.720-0.812 | **0.150 -> 0.092** |
| novel | 0.680-0.785 | 0.702-0.788 | 0.104 -> 0.087 |
| Overlap_uORF | 0.665-0.771 | 0.660-0.757 | 0.105 -> 0.097 |
| dORF | 0.480-0.606 | 0.478-0.553 | 0.126 -> 0.075 |
| internal | 0.512-0.617 | 0.500-0.573 | 0.106 -> 0.073 |

**About a third of E2's ceiling WIDTH is depth, not biology** (uORF 0.150 -> 0.092, dORF 0.126 ->
0.075). The band's POSITION barely moves, so **no conclusion in E2 changes**: annotated still sits at
its ceiling, dORF and internal still fail badly, and the model's deficit is unchanged. What changes
is the honest description of the band -- it is narrower than the natural-depth version implies, and
part of its width was never protocol variability at all.

### What this does NOT license

Depth explains about half the disagreement, **not all of it**. The residual wang gap is real, so
"the experiments only disagree because of depth" would be as wrong as ignoring depth entirely. And
this is mouse liver at one target depth; it is not a general statement about Ribo-seq reproducibility.

Outputs: `results/liver3x3_matched_depth/{ribo_vs_ribo_matched_depth.tsv,.json,e2_ceiling_matched_depth.json}`.

---

## 2026-08-16 -- RNA provenance columns in the dataset census, and what the pack actually stores

Added `n_rna_libraries`, `rna_coverage_total`, `rna_source` and `accession_rna` to
`results/dataset_census/dataset_census.{tsv,json}` for all 16 rows.

### `rna_coverage_total` is a RAW total, so it may be labelled as depth

Verified, not assumed: it is exactly `coverage.npy.sum()` over the int32 per-nt counts. No
normalisation and no depth-correction. A pack's `global_mean_coverage` is exactly this value
divided by `sum_L`, i.e. the per-nt mean of the same quantity.

Caveat that affects the axis: it is **summed per-nt read depth over transcript positions**, while
`psites` is a **count of footprints**. Different units, separate axes. Both definitions are now
written into `dataset_census.json` under `quantities` so the labelling travels with the data.

### `n_rna_samples` reads 0 in every `_canon` pack despite real coverage

`coverage_norm.json` carries an `n_rna_samples` field that looks like the obvious source for a
library count. **It is 0 in every `_canon` pack**, including Fibroblast at 1,199,211,035,619
coverage. Those packs reused their union twin's RNA verbatim -- correct, because the final recipe
does not touch RNA -- but the rebuild never carried the count forward. Reading the field directly
would have reported "0 RNA libraries" across the whole training set.

Counts now resolve in three tiers: the provenance file list, then `n_rna_samples` when non-zero,
then the union twin accepted **only** when `coverage_total` matches exactly. Every row records
which tier applied in `rna_source`. Unresolved rows are left **empty, never 0**; HBL-1 is the one
such row (no pack in this project).

Two rows that look wrong and are not: **GSE39561** shows GSE208041's exact RNA totals because it
has no RNA of its own and borrows them; **B721.221** needed an explicit pack override because its
census row has `pack=None` (its depth comes from the external Ouspenskaia reference) while RNA is
its entire model input.

Fibroblast's 32 RNA libraries coincidentally equals its 32 Ribo libraries. Checked: provenance
lists 32 real RNA files with **zero overlap** against the Ribo samplesheet.

`accession_rna` is **GSE182372** for all nine Chothani tissues, distinct from Ribo-seq GSE182371.

### Mouse-liver model call runs on GSE243134 RNA -- the shallowest of the three arms

The fixed E3/P9 cell is `attn` / `rna_input=gse243134` / `pred_preddepth`. Because
`pred_preddepth` is standalone, that RNA is the *only* experimental input behind those recall
numbers. All three arms sit on the same 22,974-tx universe, so the totals are comparable:

| arm | RNA libraries | RNA coverage total | mean/nt | tx > 0 | top-100 share |
|---|--:|--:|--:|--:|--:|
| gse243134 | 19 | 3,242,315,370 | 57.54 | 22,503 | 31.1% |
| janich | 7 | 6,853,032,580 | 121.62 | 22,569 | 30.7% |
| wang | 2 | 6,884,416,652 | 122.17 | 22,605 | 39.6% |

**Library count and coverage run in opposite directions**: 19 libraries buys roughly half the
coverage of 2. `n_rna_libraries` is therefore NOT a depth proxy on this dataset. The upside is
that the winning model cell runs on the shallowest RNA, so the result is not bought with deeper
input.

**Correction to E3's stated mechanism.** Its rationale said Wang RNA collapses uORF recall because
"Wang is the shallowest RNA arm (2 samples)". The sample count is right; **the depth claim is
backwards**. Wang is the deepest arm by total and by mean/nt and covers the most transcripts. What
separates it is evenness: its top 100 transcripts hold 39.6% of coverage against ~31% for the
other two, from only 2 libraries. So the uORF collapse is a concentration/replicate effect, not a
depth one, and the precise mechanism is **not established** by these numbers -- recorded as
unexplained rather than replacing one guess with another. Comment fixed in
`make_shallow_experiment_win.py`; no plotted number changes.
