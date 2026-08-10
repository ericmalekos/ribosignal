# Figures -- master index + regeneration commands

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
| examples | `prediction_examples` | BUILT | observed vs predicted per-nt profile exemplars (uORF / CDS / lncRNA), held-out Hepatocytes | `$PY scripts/plot_prediction_examples.py` |

## B -- rigor / honesty

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| B4 | `B4_replicate_ceiling` | BUILT (mamba4, D1b) | replicate-concordance ceiling vs released model: captures 73% (pc) / 54% (lncRNA) / 63% (uORF 5'UTR) / 36% (dORF 3'UTR) of the reproducible signal (attn supplemental: 70/51/61/34) | `cd B4_replicate_ceiling && $PY make_replicate_ceiling.py` |
| B5 | `B5_expression_independence` | BUILT (mamba4, D1b) | translation discrimination exceeds an expression-only baseline (periodicity +0.110 all / +0.067 non-canonical, length-controlled) | `cd B5_expression_independence && $PY make_expression_independence.py` |
| **B6** | `B6_input_ablation` | BUILT (released ckpt) | **input ablation on the deployed recipe** -- sequence-only keeps 98.5% of profile SHAPE but the count head loses 0.193 Pearson without RNA-seq. Cell-type specificity enters through the COUNT HEAD, not the shape. Closes task 61 / FIGURES_PLAN item 1 | `cd B6_input_ablation && $PY make_input_ablation.py` |
| B7 | `B7_multimap_posture` | BUILT | RNA-seq mm20 vs mm1 coverage identical on ~95% of tx (median Pearson 1.0000) -> mm1 retrain unwarranted | `cd B7_multimap_posture && $PY make_multimap_posture.py` |
| B8 | `S_fdr_rigor` | BUILT (supersedes `B8_fdr_calibration`) | a global cut does not control the non-canonical class in EITHER direction: 37.5x over-report on A549 (tryptic, 65k canonical peptides), down to 0.2x UNDER-report on DoHH2. The error tracks canonical-class size, not database size | `cd S_fdr_rigor && $PY make_fdr_rigor.py` |
| S-QC | `S_riboseq_qc` | BUILT | Ribo-seq library quality for all 16 dataset arms (9 training + 5 held-out + 2 cell lines): f0 74.8-87.9% vs a 33% floor, footprints 28-31 nt, offset 12 nt. Closes D5 | `cd S_riboseq_qc && $PY make_riboseq_qc.py` |
| S-RS | `S_rescore_substrate` | BUILT | MS2Rescore is substrate-dependent: 0.16 confident novel per 1k novel PSMs on tryptic (rescoring collapses 10 -> 1) vs 21.9-100.7 on HLA-I (10 -> 26). Closes rigor point 4 | `cd S_rescore_substrate && $PY make_rescore_substrate.py` |

## C -- interpretability

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| C9 | `kozak` | BUILT | model V3 LEARNS the Kozak start context (learned kernel vs empirical PWM vs heuristic; r=0.807) | `$PY scripts/plot_kozak_weights.py` |
| C10 | `C10_saliency` | BUILT (released ckpt) | input saliency: the start codon is among the most salient nt in the whole transcript (rank 3/2036 uORF, 12/7914 lncRNA; annotated CDS is weaker at 319/2094), in-ORF saliency 3.4-8.8x mean -- the model reads the start context, not just coverage | `cd C10_saliency && $PY make_saliency.py` |
| depth | `depth_crossover` | BUILT | predicted profile beats measured Ribo-seq for ORF calling below ~29M P-sites (the "why not just do Ribo-seq" panel) | see `depth_crossover/FIGURE_DATA_INPUTS.md` (SLURM subsample grid + `plot_depth_crossover.py`) |

## D -- comparisons / application

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| D11 | `D11_vs_seq2ribo` | BUILT | beats the purpose-built seq2ribo on within-transcript profile shape (0.526 vs published 0.05-0.19) | `cd D11_vs_seq2ribo && $PY make_vs_seq2ribo.py` |
| **F2b** | `F2b_discovery_forest` | BUILT | **the Fig 2 power fix.** Model/null discovery-density ratio, 5 datasets x 2 models x 2 arms = 20 points, ALL above 1.0 (median 93.6x, range 10.7-2744x). Turns the one-off into a reproducible pattern | `cd F2b_discovery_forest && $PY make_f2b_forest.py` |
| D12 | `D12_cpat_cpc2` | BUILT | model vs CPAT/CPC2 (locked decision D3). Splits by assay: CPAT/CPC2 win tryptic (25/23 vs 11), the model wins all three HLA-I (14v10, 22v14, 32v8). On per-1k-sequence efficiency the Poisson arm leads all four (4.6-7.6 vs 0.03-0.11 for the null) | `cd D12_cpat_cpc2 && $PY make_cpat_cpc2.py` |
| D13 | `D13_discovery_errorbars` | BUILT | DB-design tradeoff, 3 panels: GAIN (novel peptides per 1k DB seqs), COST (dPSM vs the gencode baseline), DRIVER (DB size). Current generator is `make_d13_pgx.py` on the released models + pgx + frozen params; `make_discovery_errorbars.py` is the f0-era version, retained and labelled SUPERSEDED | `cd D13_discovery_errorbars && $PY make_d13_pgx.py` |
| D14 | `D14_novel_antigen` | GATED (MS) | 2-3 highlighted novel PSMs: ms2pip mirror plots + ORF class + genomic context (Fig 2c) | build with Fig 2 MS pass |

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
