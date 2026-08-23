# Figures -- master index + regeneration commands

> Current vs superseded results: **`docs/STATUS_CURRENT_VS_ARCHIVED.md`**.

Every figure lives in its own subfolder with a reusable generator, a `FIGURE_DATA_INPUTS.md` (panel ->
data provenance), and PDF+PNG output. This index is the single place to see status + the exact command to
(re)produce each one. Env for all generators: `PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3`.

Figures are grouped A (validity / generalization), B (rigor / honesty), C (interpretability),
D (comparisons / application). Numbering matches the working figure discussion; the two *main* figures
(Fig 1 validity, Fig 2 utility) are composed from these.

Regenerate any built figure with: `cd figures/<folder> && $PY make_<name>.py` (exact command in each
folder's FIGURE_DATA_INPUTS.md "## Regenerate" block).

## A -- validity / generalization

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| A1 | `A1_localization_ceiling` | BUILT (mamba4, D1b) | predicted translation-localization BEATS the observed Ribo-seq periodicity ceiling on all-ORFs (Hep 0.886 vs 0.828; Ruiz-Orera 0.892 vs 0.825); on NON-CANONICAL it reaches 94-96% of ceiling but does not beat it -- the caption must say both | `cd A1_localization_ceiling && $PY make_localization_ceiling.py` |
| A2 | `A2_ribocode_dropin` | BUILT (mamba4, D1b) | predicted profile dropped into RiboCode reproduces its own ORF calls, F1 0.913 (attn supplemental 0.909), no Ribo-seq for the query | `cd A2_ribocode_dropin && $PY make_ribocode_dropin.py` |
| A3 | `A3_loto_spread` | BUILT (**pre-union recipe, labelled**) | even cross-tissue generalization; profile-Pearson spread (0.23-0.64) is held-out target quality (r=0.81 vs period_obs), count transfer even (CV 0.11) | `cd A3_loto_spread && $PY make_loto_spread.py` |
| **A4** | `A4_b721_ground_truth` | BUILT (**new 2026-08-14**) | **the only NON-NULL benchmark, and it had no figure.** B721.221: ORF calls predicted from Sarkizova RNA-seq vs **327 M measured** Ouspenskaia Ribo-seq footprints, 16,331 measured calls in scope. Against the assay's own SPLIT-HALF ceiling (F1 0.889, annot 0.973, **non-canon 0.498**) the model reaches **75% / 78% / 48-54%**. Poisson buys precision (0.78->0.89, 0.76->0.92) and pays recall -- the calibration dial confirmed against MEASURED translation. attn F1 0.666 vs mamba4 0.665, a third tie. Never quote the 0.24 non-canonical recall without the 0.498 denominator | `cd A4_b721_ground_truth && $PY make_b721_ground_truth.py` |
| examples | `prediction_examples` | BUILT | observed vs predicted per-nt profile exemplars (uORF / CDS / lncRNA), held-out Hepatocytes | `$PY scripts/plot_prediction_examples.py` |

## B -- rigor / honesty

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| B4 | `B4_replicate_ceiling` | BUILT (mamba4, D1b) | replicate-concordance ceiling vs released model: captures 73% (pc) / 54% (lncRNA) / 63% (uORF 5'UTR) / 36% (dORF 3'UTR) of the reproducible signal (attn supplemental: 70/51/61/34) | `cd B4_replicate_ceiling && $PY make_replicate_ceiling.py` |
| B5 | `B5_expression_independence` | BUILT (mamba4, D1b) | translation discrimination exceeds an expression-only baseline (periodicity +0.110 all / +0.067 non-canonical, length-controlled) | `cd B5_expression_independence && $PY make_expression_independence.py` |
| **B6** | `B6_input_ablation` | BUILT (released ckpt) | **input ablation on the deployed recipe** -- sequence-only keeps 98.5% of profile SHAPE but the count head loses 0.193 Pearson without RNA-seq. Cell-type specificity enters through the COUNT HEAD, not the shape. Closes task 61 / FIGURES_PLAN item 1 | `cd B6_input_ablation && $PY make_input_ablation.py` |
| B7 | `B7_multimap_posture` | BUILT | RNA-seq mm20 vs mm1 coverage identical on ~95% of tx (median Pearson 1.0000) -> mm1 retrain unwarranted | `cd B7_multimap_posture && $PY make_multimap_posture.py` |
| B8 | `S_fdr_rigor` | BUILT (supersedes `B8_fdr_calibration`) | a global cut does not control the non-canonical class in EITHER direction: 37.5x over-report on A549 (tryptic, 65k canonical peptides), down to 0.2x UNDER-report on DoHH2. The error tracks canonical-class size, not database size | `cd S_fdr_rigor && $PY make_fdr_rigor.py` |
| **B9** | `B9_overcall_validation` | BUILT (**new 2026-08-14**) | **the measurement behind P9's headline.** Are the model's 1,942 experiment-unsupported mouse-liver ORF calls false positives or depth-limited discoveries? Against a 28-library merged RiboCode reference: ANCHOR (experiment-supported model calls) **99.3%**, CONTROL (observed calls unique to one experiment) **74.5%**, the model's extras **7.7%** -> **~92% false positives**. Two architectures agreeing does NOT rescue them (8.7%). Gene-coverage split closes the 'no data there' objection: 1,061 of 1,942 sit on genes with zero translation in 28 pooled libraries | `$PY scripts/score_merged_liver_overcall.py` then `cd B9_overcall_validation && $PY make_overcall_validation.py` |
| S-QC | `S_riboseq_qc` | BUILT (**18 arms 2026-08-15**, was 15) | Ribo-seq library quality for all **18** dataset arms (8 training + 6 held-out + 2 pgx + **2 dropped**): read-length distribution, pooled frame-0 at annotated CDS starts, P-site offset vs read length. New: THP-1 GSE208041 **84.0%**, CAR-T **78.5%**, Brain **76.2%**. `dropped` is a new grey role for the project's EXCLUDED datasets. THP-1 GSE39561 CANNOT be added -- its pre_config has headers and zero data rows, so it lands in `missing` (use the 40.7% prose value, which is the same quantity). **Brain is above HUVEC on frame-0 (76.2 vs 74.6%) and 2x its depth**, so its exclusion rests on `period_obs` 0.044, not on library quality. `psites_at_cds` renamed **`psites_at_start_window`** -- it is a start-window sum at ~0.5% of depth, NOT usable as depth | `cd S_riboseq_qc && $PY make_riboseq_qc.py` |
| **S-REPRO** | `S_reproducibility` | BUILT (canonical) | **can the published pipeline regenerate the training data?** All 74 Chothani Ribo runs re-aligned and re-pooled per tissue. 3.04-11.66% of signal moves, median 3.86%, per-nt r 0.432-0.980. HUVEC -- the one tissue the claim used to rest on -- is the BEST of the eight, so the single-tissue number understated it. Panel B shows \|divergence\| and r are partly decoupled | `cd S_reproducibility && $PY make_reproducibility.py` |
| S-RS | `S_rescore_substrate` | BUILT | MS2Rescore is substrate-dependent: 0.16 confident novel per 1k novel PSMs on tryptic (rescoring collapses 10 -> 1) vs 21.9-100.7 on HLA-I (10 -> 26). Closes rigor point 4 | `cd S_rescore_substrate && $PY make_rescore_substrate.py` |

## C -- interpretability

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| C9 | `kozak` | BUILT | model V3 LEARNS the Kozak start context (learned kernel vs empirical PWM vs heuristic; r=0.807) | `$PY scripts/plot_kozak_weights.py` |
| C10 | `C10_saliency` | BUILT (released ckpt) | input saliency: the start codon is among the most salient nt in the whole transcript (rank 3/2036 uORF, 12/7914 lncRNA; annotated CDS is weaker at 319/2094), in-ORF saliency 3.4-8.8x mean -- the model reads the start context, not just coverage | `cd C10_saliency && $PY make_saliency.py` |
| depth | `depth_crossover` | BUILT | predicted profile beats measured Ribo-seq for ORF calling below ~29M P-sites (the "why not just do Ribo-seq" panel) | see `depth_crossover/FIGURE_DATA_INPUTS.md` (SLURM subsample grid + `plot_depth_crossover.py`) |
| **C11** | `C11_channel_ablation` | BUILT (canonical) | **ORF-track channel ablation, 6 arms x 2 models.** Double dissociation: the FRAME channels *suppress* internal ORFs (removing any one roughly triples internal F1), while the START channel is worth ~0.000 on annotated CDS and -0.069/-0.073 on uORFs. Also the reason the planned start-codon ISM sweep was cancelled | `cd C11_channel_ablation && $PY make_channel_ablation.py` |
| **C12** | `C12_codon_occupancy` | BUILT (canonical) | **codon occupancy the model was never given.** Against the right ceiling: mouse 0.263 = 98% of the across-species ceiling (0.269, never-seen transcripts), human 0.548 = 58% of within-species (0.944, but ~99.7% seen). Recipe-invariant (r=0.99952 across alignment recipes) | `cd C12_codon_occupancy && $PY make_codon_occupancy.py` |
| **C13** | `C13_memorisation` | BUILT | **how much is memorisation, at matched depth.** The naive all-seen-vs-never-seen gap (0.468/0.492) is ~4x inflated by depth; matched on log10 counts + length it is 0.100/0.108 with disjoint CIs. The honesty panel for the interpretability section | `cd C13_memorisation && $PY make_memorisation.py` |

## D -- comparisons / application

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| D11 | `D11_vs_seq2ribo` | BUILT | beats the purpose-built seq2ribo on within-transcript profile shape (0.526 vs published 0.05-0.19) | `cd D11_vs_seq2ribo && $PY make_vs_seq2ribo.py` |
| **F2b** | `F2b_discovery_forest` | BUILT (**2 panels**) | **the Fig 2 power fix.** (a) model/null discovery-density ratio, 5 datasets x 2 models x 2 arms = 20 points, ALL above 1.0 (median 93.6x, range 10.7-2744x). (b) the ABSOLUTE bottom line added 2026-08-14: change in TOTAL unique peptides vs GENCODE-only -- Poisson arm positive in **9/10**, theta=1 in 6/10, naive AUG null in **2/5** (A549 **-2,884**). The one Poisson failure, A549/attn **-35**, must stay in any caption | `cd F2b_discovery_forest && $PY make_f2b_forest.py` |
| D12 | `D12_cpat_cpc2` | BUILT (**4 panels**) | model vs CPAT/CPC2 (locked decision D3). Splits by assay: CPAT/CPC2 win tryptic (25/23 vs 11), the model wins all three HLA-I (14v10, 22v14, 32v8). On per-1k-sequence efficiency the Poisson arm leads all four (4.6-7.6 vs 0.03-0.11 for the null). Panel (d) added 2026-08-14: on TOTAL unique peptides the mamba4 Poisson arm is the only arm positive on all four (+11/+14/+10/+21) -- **mamba4-specific**, attn is -35 on A549 (see D13) | `cd D12_cpat_cpc2 && $PY make_cpat_cpc2.py` |
| D13 | `D13_discovery_errorbars` | BUILT (**4 panels**) | DB-design tradeoff: GAIN (novel peptides per 1k DB seqs), COST (dPSM vs the gencode baseline), DRIVER (DB size), and TOTAL (change in TOTAL unique peptides, added 2026-08-14 -- carries both released models, and the near-cognate null's **-7,780** on A549, the largest loss measured anywhere here). Current generator is `make_d13_pgx.py` on the released models + pgx + frozen params; `make_discovery_errorbars.py` is the f0-era version, retained and labelled SUPERSEDED | `cd D13_discovery_errorbars && $PY make_d13_pgx.py` |
| D14 | `D14_novel_antigen` | GATED (MS) | 2-3 highlighted novel PSMs: ms2pip mirror plots + ORF class + genomic context (Fig 2c) | build with Fig 2 MS pass |

## E -- the mouse-liver 3x3 factorial

3 Ribo-seq datasets x 3 RNA-seq inputs x 2 models on one shared universe (22,974 tx), so a
difference in these panels is a difference between EXPERIMENTS rather than between code paths.

**All three were rebuilt on canonical alignments 2026-08-13** and now read
`results/mouse_liver_3x3_canon/scored/`. They had been reading the off-recipe
`results/mouse_liver_3x3/scored/` and were missing from this index entirely, which is why nothing
flagged them as stale -- an unindexed figure is an unmaintained one, so new figures go in this table
when they are built, not later.

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| E1 | `E1_rna_provenance` | BUILT (canonical) | matching the RNA-seq to the Ribo-seq reference is worth only **+0.0035 mean F1** (range -0.0042 to +0.0070, negative in 2 of 12 cells) -- RNA provenance is a second-order effect | `cd E1_rna_provenance && $PY make_rna_provenance.py` |
| E2 | `E2_class_ceiling` | BUILT (canonical) | per-class between-experiment ceiling vs model. The model tracks the ceiling on uORF/Overlap_uORF but sits far below it on internal (0.022 vs 0.512-0.617) and dORF -- the ceiling, not 1.0, is the bar | `cd E2_class_ceiling && $PY make_class_ceiling.py` |
| E3 | `E3_shallow_experiment_win` | BUILT (canonical) | **the model beats a real (shallow) Ribo-seq experiment** on novel (+0.107) and uORF (+0.142, +0.103) recall vs Wang. Most Ribo-seq is shallow, so predicting better than a shallow experiment for free is the applied argument | `cd E3_shallow_experiment_win && $PY make_shallow_experiment_win.py` |

## P -- poster panels (48 x 36 in landscape, 4 columns)

Built for the poster; each is also a standalone figure. The other ten poster panels reuse existing
figures (`arch_attn`, `S_riboseq_qc`, `A1`, `A2`, `E1`, `E2`, `E3`, `F2b`, `D12`). Poster reports the
DEPLOYED checkpoints, not the canonical retrains.

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| **P3** | `P3_training_data` | BUILT | training composition: 7 cell types / **69** Ribo-seq libraries, Hepatocytes (5) held out entirely, 84,472-tx universe, **11x depth spread** (Fibroblast 32 libs vs HUVEC 3). Panel B states the inference inputs, because "predicts from sequence" is ambiguous about whether Ribo-seq is used | `cd P3_training_data && $PY make_training_data.py` |
| **P5** | `P5_heldout_human` | BUILT (**CAR-T provisional -- check `cart_source_is_canonical`**) | **the central validation claim.** 4 independent human held-outs (hepatocyte LOTO / iPSC-CM cross-study / THP-1 / CAR-T), both arms. Standalone F1 0.849-0.913. Panel B is the honest half: ~0.99 recall on annotated CDS vs **0.55-0.60 non-canonical** | `cd P5_heldout_human && $PY make_heldout_human_panel.py` |
| **P9** | `P9_upset_mouse3x3` | BUILT (**8 outputs**: 2 model families x 4 ORF-class views) | where the model's disagreement with observed Ribo-seq lives. Three mouse-liver experiments + the standalone arm, one universe (22,974 tx), genomic keying. The "model only" bar is a MEASUREMENT, not a guess: a 28-library merged RiboCode reference corroborates **7.7%** of mamba4's 1,942 experiment-unsupported calls, against a **74.5%** control (observed calls equally unsupported by the other arms) and a **99.3%** anchor -- so **~92% are false positives**. Both models agreeing does NOT rescue them (8.7%): shared architecture error, never a confidence filter | `$PY scripts/score_merged_liver_overcall.py` then `cd P9_upset_mouse3x3 && $PY make_upset_mouse3x3.py` |
| **P11** | `P11_bmdm_proteomics` | BUILT (frozen search params, **12 populations**, **30 AA + ATTENTION** 2026-08-15) | mouse macrophage proteogenomic lift across all 12 populations: DB median **78x smaller**, discovery density median **43x higher** than a naive AUG null. On TOTAL unique peptides the model is **cost-neutral** (above the GENCODE-only baseline in 6/12, median **+3** of ~55,000) while the AUG null is **0/12**, median **-1,100**. SUPERSEDES the BMDM-only version (+31 was the best of 12 and read as general) AND the 7-aa mamba version, whose 156x/111x ratios were inflated by sub-30-aa null sequences the ORF-call track excludes by rule. Reproduce the old panel with `--source .../pgx_xsubtype/frozen_reports --arm model_predicted` | `cd P11_bmdm_proteomics && $PY make_bmdm_proteomics.py` |
| **P12** | `P12_macrophage_xsubtype` | BUILT (new 2026-08-14, **30 AA + ATTENTION** 2026-08-15) | **the full macrophage sweep: a PREDICTED database vs a REAL Ribo-seq one**, 12 populations x 3 arms, 18 mzML each (depth-balanced). BMDM is the MATCHED case (the Ribo-seq DBs came from it) and is excluded from every median. On the 11 TRANSFER populations a fixed **506-seq** BMDM Ribo-seq database matches the per-population model on novel yield (26 vs 24), beats it on total unique peptides (+18 vs -3; 8/11 vs 5/11) and is **2.4x denser** (51.4 vs 21.3 per 1k). The model's case is NOT needing Ribo-seq, not beating it. Correcting the floor narrowed the density gap (2.9x -> 2.4x) but widened the total-peptide gap. DB sizes are measured at plot time, never hardcoded | `cd P12_macrophage_xsubtype && $PY make_macrophage_xsubtype.py` |
| **P13** | `P13_macrophage_psm_venn` | BUILT (new 2026-08-14, **30 AA + ATTENTION** 2026-08-15) | **per-population Venn of NOVEL PSMs**, model-predicted DB vs the real BMDM-NT Ribo-seq DB (506 seqs), 12 populations. Unit is `(spec_id, peptide)`, so a scan is shared only if both arms assigned the SAME peptide; same-scan-different-peptide is annotated separately and **totals 1 across all 12** (one spectrum, LungResident -- assignments are stable, so the disagreement is about which spectra clear FDR). Transfer totals (n=11): model 866 PSMs, Ribo-seq 1,065, shared 386 = **25.0%** of the union, overlap in 11/11 but Jaccard never above 0.34 -- the two databases are COMPLEMENTARY. Held to within 1 point of the 7-aa mamba build (24.0%), so this is a property of predicted-vs-measured, not of one build | `$PY proteogenomics/scripts/macro_novel_psm_venn.py && cd P13_macrophage_psm_venn && $PY make_macrophage_psm_venn.py` |

## Poster layout

`figures/make_poster_layout.py` -> `docs/POSTER_LAYOUT.md`: the 4-column, 13-panel layout with every
headline number read LIVE from each panel's `*_values.json`, so the layout cannot drift from the
figures. Each panel carries a MUST-STATE line recording the narrower claim that survived
measurement, plus a "numbers deliberately NOT on the poster" list. Re-run after regenerating any
panel; it warns if P5's CAR-T arm is still the pre-canonical pack.

## Plotted data as TSV (every figure)

`figures/export_figure_tsvs.py` converts each generator's `*_values.json` into a TSV in the same
folder, so the numbers behind a figure can be read without parsing nested JSON or re-running a model.
It reads the JSON the generator wrote AFTER drawing, so a TSV cannot drift from its figure, and new
figures are picked up automatically. Run `$PY export_figure_tsvs.py` after regenerating anything.

Output is WIDE (one row per record) where the JSON holds a record table, LONG (`path`/`value`)
otherwise; the shape is noted in each TSV's header comment. Four generators write their own richer
TSVs directly because their plotted data is per-nucleotide or per-region rather than a summary:
`prediction_examples`, `profile_exemplars`, `o2_venn`, `benchmarks`.

## Profile exemplars and the encoding panel (added 2026-08-09)

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| exemplars | `profile_exemplars` | BUILT (mamba4) | observed-vs-predicted profiles chosen so agreement SURVIVES deleting the top observed peaks; human + mouse, per ORF class. Replaces `prediction_examples`, two of whose three exemplars were spike-dominated (RPL32 top1_frac **0.538**) | `$PY scripts/select_profile_exemplars.py` then `$PY scripts/plot_profile_exemplars.py --tsv ... --per-class 1` |
| S-ENC | `S_encoding` | BUILT (pre-union recipe) | one-hot vs RiNALMo / Orthrus / HydraRNA, identical recipe and n=33,918. one-hot best on pc profile (0.6250); RiNALMo wins lncRNA and periodicity. No retraining needed -- Task 41 data | `cd S_encoding && $PY make_encoding_panel.py` |

## Main figures (assembled)

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| **Fig 1** | `main/` | BUILT (mamba4, D1b) | validity, 5 panels: exemplars / drop-in F1 0.913 / localization vs ceiling / depth crossover / LOTO spread. 7.20 x 10.19 in | `cd main && $PY make_main_figures.py --only Figure1_validity` |
| **Fig 2** | `main/` | BUILT (5 datasets) | utility, 3 panels: discovery forest / vs CPAT-CPC2 / DB tradeoff. 7.20 x 8.67 in. No D14 PSM panel and no cross-species MS arm -- task 62 superseded at 5 of 6 datasets | `cd main && $PY make_main_figures.py --only Figure2_utility` |

Panels are composited as vector PDF pages (via `pypdf`), not re-rastered, and each figure writes a
`*_manifest.json` naming every source panel and its mtime -- so a stale panel shows up as a stale
timestamp. Layout is a data structure in the generator; re-specifying a figure means editing one list.
**Read `main/FIGURE_DATA_INPUTS.md` before using these**: it lists five known limitations, including
that panel 1e is still on the pre-union recipe and that Fig 1 is taller than one printed page.

## Schematics (not data figures)

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| arch | `arch_attn` | BUILT | model architecture panels for BOTH candidate models -- `arch_rinalmo` (attn, 2 transformer blocks, 5,071,106 params) and `arch_rinalmo_mamba` (mamba4, 4 Bi-Mamba blocks, 7,521,026 params), identical except the mixer region so they sit side by side. RiNALMo Fig-1 grammar for the network, Orthrus Fig-1 grammar for the inputs; the ORF-track example is computed by `build_orf_track.py` at render time | `cd arch_attn && $PY make_arch_panel.py` |

## Checkpoint provenance (audited 2026-08-08)

Run directories are named for the DATASET, not the CHECKPOINT, so reading a figure's source path tells
you nothing about which model produced it. A scan of every generator -- not only the ones under
suspicion -- found that **five of the six Fig 1 panels were built from
`orf_v2_attn_onehot_holdout_Hepatocytes`**, the pre-nokozak / pre-mm1 / pre-union checkpoint, while the
project ships the union models. Only `arch_attn` was on a released model.

| generator | before | now |
|---|---|---|
| A1 localization ceiling | stale | rebuilt from `eval_localization_released.sbatch` output |
| A2 drop-in | stale (F1 0.923) | released (F1 **0.909**); asserts the npz `meta` tag before plotting |
| A3 LOTO spread | `results/loto_9fold/onehot_orf_v2_attn_holdout_*` | **unchanged and labelled** -- no union 9-fold exists; re-running means a 9-fold retrain |
| B4 replicate ceiling | hardcoded pre-union dict | reads `extra_metrics.json` of the released run |
| B5 expression independence | stale | rebuilt from the released `localization_metrics.json` |
| C10 saliency | stale ckpt + Kozak-heuristic track | released ckpt + `orf_track_v2_nokozak.npy` (both had to move together) |
| prediction_examples | stale | released |
| arch_attn | released | unchanged |

Rule going forward: a generator that quotes a model number either asserts the checkpoint tag (A2's
`check_checkpoint()`) or reads it from a file and prints it (B4, and the two tutorial held-out pages).
Hardcoded model numbers are how this happened.

### A second, independent problem found the same day: the WRONG MODEL

Locked decision **D1b (2026-07-31) makes `orf_v2_mamba4` the main Figure 1 model**, with `orf_v2_attn`
moving to supplemental while staying shipped as the CPU inference path. Every Fig 1 panel was still
built on attn: the decision was recorded in the plan and never applied to the generators. The older
"lead with the transformer" rationale (attn 0.639 vs mamba 0.599 whole-tx Pearson) is pre-union and was
what D1b superseded, after the union 3-seed comparison reversed it (mamba4 0.6799 vs attn 0.6595, with
non-overlapping seed ranges).

A1, A2, B4 and B5 now take `FIG_MODEL`, default **mamba4**, and write the model into their values JSON.
`FIG_MODEL=attn` builds the supplemental variant to `*_attn.pdf` / `*_values_attn.json` -- a separate
path, because the first attn build silently overwrote the mamba4 PDF at the shared filename.

Two panels stay on attn and must say so in their captions:

- **C10 saliency** -- documented CPU-only, and mamba_ssm needs CUDA. D1b itself names attn as the CPU
  inference path, so this is coherent rather than a compromise.
- **B6 input ablation** -- both ablation arms were trained with attn. A mamba4 version means two more
  12-17 h trainings. The double dissociation is unlikely to be architecture-specific, but that has not
  been measured and should not be asserted.

## Status summary (2026-08-07)
- BUILT (16): A1, A2, A3, B4, B5, B7, C10, D11, D13 + prediction_examples, kozak, depth_crossover,
  arch_attn, and four added on 2026-08-07: **D12** (CPAT/CPC2, closes D3), **S_fdr_rigor** (replaces
  the gated B8), **S_riboseq_qc** (closes D5), **S_rescore_substrate** (closes rigor point 4).
- STILL GATED on MS: **D14** only (novel-PSM mirror plots). Everything else the MS pipeline was
  blocking is now built on the released mamba4 + pgx + frozen search params.
- `B8_fdr_calibration/` is retained as the empty placeholder folder; the built figure is
  `S_fdr_rigor/`, whose measurement corrected the placeholder's "global inflates ~10x" one-liner --
  the error runs in both directions and is driven by canonical-class size (see its
  FIGURE_DATA_INPUTS.md).
- Immunopeptidome panel is 4 of 6 as of this date (A549 tryptic + HBL-1 / SU-DHL-4 / DoHH2 HLA-I,
  with THP-1 in progress). B721.221 MS is BLOCKED (MassIVE FTP unreachable from prism); B721 still
  contributes the drop-in ground-truth check.

## Notes
- BUILT figures (prediction_examples, kozak, depth_crossover) predate the make_<name>.py + PDF
  convention; their generators live in `scripts/` and their FIGURE_DATA_INPUTS.md already carry full regen
  commands. PDF versions can be produced by re-running those generators with a `.pdf` output path.
- All numbers are read from result JSON/TSV files, so every BUILT generator regenerates automatically if
  the model is retrained (e.g. the now-closed mm1 RNA-coverage ablation).
