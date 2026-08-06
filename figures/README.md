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
| A1 | `A1_localization_ceiling` | BUILT | predicted translation-localization matches/beats the observed Ribo-seq periodicity ceiling (Hep 0.891 vs 0.828; Ruiz-Orera 0.894 vs 0.825) | `cd A1_localization_ceiling && $PY make_localization_ceiling.py` |
| A2 | `A2_ribocode_dropin` | BUILT | predicted profile dropped into RiboCode reproduces its own ORF calls, F1 0.923, no Ribo-seq for the query | `cd A2_ribocode_dropin && $PY make_ribocode_dropin.py` |
| A3 | `A3_loto_spread` | BUILT | even cross-tissue generalization; profile-Pearson spread (0.23-0.64) is held-out target quality (r=0.81 vs period_obs), count transfer even (CV 0.11) | `cd A3_loto_spread && $PY make_loto_spread.py` |
| examples | `prediction_examples` | BUILT | observed vs predicted per-nt profile exemplars (uORF / CDS / lncRNA), held-out Hepatocytes | `$PY scripts/plot_prediction_examples.py` |

## B -- rigor / honesty

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| B4 | `B4_replicate_ceiling` | BUILT | replicate-concordance ceiling vs model: captures 67% (pc) / 41% (lncRNA) / 55% (per-codon) of the reproducible signal | `cd B4_replicate_ceiling && $PY make_replicate_ceiling.py` |
| B5 | `B5_expression_independence` | BUILT | translation discrimination exceeds an expression-only baseline (periodicity +0.120 all / +0.087 non-canonical, length-controlled) | `cd B5_expression_independence && $PY make_expression_independence.py` |
| B7 | `B7_multimap_posture` | BUILT | RNA-seq mm20 vs mm1 coverage identical on ~95% of tx (median Pearson 1.0000) -> mm1 retrain unwarranted | `cd B7_multimap_posture && $PY make_multimap_posture.py` |
| B8 | `B8_fdr_calibration` | GATED (MS) | class-specific vs global FDR + decoy calibration; global inflates novel discovery ~10x. Needs the immunopeptidome search consolidation (Fig 2). | build with Fig 2 MS pass |

## C -- interpretability

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| C9 | `kozak` | BUILT | model V3 LEARNS the Kozak start context (learned kernel vs empirical PWM vs heuristic; r=0.807) | `$PY scripts/plot_kozak_weights.py` |
| C10 | `C10_saliency` | BUILT | input saliency: the start codon is among the most salient nt in the whole transcript (rank 8/2036 uORF, 15/7914 lncRNA), in-ORF saliency 3.6-13.6x mean -- the model reads the start context, not just coverage | `cd C10_saliency && $PY make_saliency.py` |
| depth | `depth_crossover` | BUILT | predicted profile beats measured Ribo-seq for ORF calling below ~29M P-sites (the "why not just do Ribo-seq" panel) | see `depth_crossover/FIGURE_DATA_INPUTS.md` (SLURM subsample grid + `plot_depth_crossover.py`) |

## D -- comparisons / application

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| D11 | `D11_vs_seq2ribo` | BUILT | beats the purpose-built seq2ribo on within-transcript profile shape (0.526 vs published 0.05-0.19) | `cd D11_vs_seq2ribo && $PY make_vs_seq2ribo.py` |
| D12 | `D12_cpat_cpc2` | GATED (MS) | model-selected ORFs vs CPAT/CPC2 coding-potential baseline in the immunopeptidome DB (locked decision D3) | build with Fig 2 MS pass |
| D13 | `D13_discovery_errorbars` | BUILT | deterministic DB-design tradeoff (raw MSFragger 1% FDR): model 2-3.4x more novel per ORF, null displaces 1.5-3.1x more CDS, null 3x decoy load = FDR-recalibration driver | `cd D13_discovery_errorbars && $PY make_discovery_errorbars.py` |
| D14 | `D14_novel_antigen` | GATED (MS) | 2-3 highlighted novel PSMs: ms2pip mirror plots + ORF class + genomic context (Fig 2c) | build with Fig 2 MS pass |

## Schematics (not data figures)

| fig | folder | status | one-line | regenerate |
|-----|--------|--------|----------|------------|
| arch | `arch_attn` | BUILT | model architecture panels for BOTH candidate models -- `arch_rinalmo` (attn, 2 transformer blocks, 5,071,106 params) and `arch_rinalmo_mamba` (mamba4, 4 Bi-Mamba blocks, 7,521,026 params), identical except the mixer region so they sit side by side. RiNALMo Fig-1 grammar for the network, Orthrus Fig-1 grammar for the inputs; the ORF-track example is computed by `build_orf_track.py` at render time | `cd arch_attn && $PY make_arch_panel.py` |

## Status summary (2026-07-21)
- BUILT (11): A1, A2, A3, B4, B5, B7, C10, D11 (this session) + prediction_examples, kozak, depth_crossover (prior).
- GATED on the immunopeptidome MS pipeline (Fig 2): B8, D12, D13, D14. The pipeline is consolidating --
  DoHH2 done, SU-DHL-4 re-running (predict PENDING on GPU), B721.221/THP-1/mouse not yet started.
- All non-MS-gated exploratory figures are now built.

## Notes
- BUILT figures (prediction_examples, kozak, depth_crossover) predate the make_<name>.py + PDF
  convention; their generators live in `scripts/` and their FIGURE_DATA_INPUTS.md already carry full regen
  commands. PDF versions can be produced by re-running those generators with a `.pdf` output path.
- All numbers are read from result JSON/TSV files, so every BUILT generator regenerates automatically if
  the model is retrained (e.g. the now-closed mm1 RNA-coverage ablation).
