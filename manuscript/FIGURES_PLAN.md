# Brief Communication -- figures plan (working draft, 2026-07-20)

Target: brief communication, **2 main figures** + supplementary. NOT building figures yet -- this is the
plan + open decisions. Numbers below are pulled from results.md (grounded, not aspirational).

## Core message (one sentence)

A lightweight, CPU-runnable model predicts ribosome P-site translation from sequence + RNA-seq alone --
matching or beating experimental Ribo-seq at low-to-moderate depth and generalizing across tissues and
species -- and this enables cell-type-specific cryptic-ORF discovery in MHC-I immunopeptidomes with no
Ribo-seq for the query sample.

- **Fig 1 = the model predicts translation and generalizes** (validity).
- **Fig 2 = the immunopeptidome application** (utility).

## Model recap (what is deployed)

- Deployed model = one-hot `orf_v2_attn`: dilated CNN + 2 transformer layers + non-AUG ORF channel, ~5.0M
  params, CPU-runnable. Inputs: one-hot sequence + per-nt RNA-seq coverage. Dual heads (per-nt profile +
  count). FM token embeddings give NO lift over one-hot (Task 15) -- one-hot is the recommendation.
- Mixer options (Task 17C): transformer (deployed) vs bidirectional Mamba. Mamba is O(L) vs O(L^2) and
  trades whole-tx Pearson DOWN (0.599 vs 0.639) for periodicity UP; ~5.77M params. Not equivalent -- a
  tradeoff, and its niche is very long transcripts.

---

## PROPOSED STRUCTURE (user proposal + refinements flagged inline)

### Main Fig 1 -- the model predicts translation and generalizes
User proposed: single architecture (they suggested Mamba) + one-hot vs token-emb + exemplar uORF/CDS/lncRNA
from LOTO Hepatocyte predictions (predicted vs observed concordance).

**REFINEMENT (important): lead with the deployed TRANSFORMER `orf_v2_attn`, not Mamba.** Rationale:
1. It is the deployed/recommended model -- consistency with every downstream result.
2. It is BETTER on whole-tx Pearson (0.639 vs 0.599), which is exactly the quantity the exemplar
   predicted-vs-observed concordance panels display; Mamba would make the exemplars look visibly worse.
3. Mamba's advantages (periodicity, O(L)) are niche; params are nearly identical (5.0M vs 5.77M), so
   "least resource intensive" is not a strong differentiator here.
   -> Put the transformer-vs-Mamba mixer comparison in Supp 1-2, framed as "Mamba is a viable O(L)
   alternative for very long transcripts, trading whole-tx Pearson for periodicity."

**REFINEMENT: Fig 1 is overloaded** (architecture + encoding + 3 exemplars). Move architecture + encoding
grids fully to Supp. Proposed Main Fig 1 panels:
- (a) LOTO cross-tissue generalization: 9-fold spread of the key metric(s) (drop-in F1 ~0.92, localization
  AUROC), showing even generalization across held-out tissues.
- (b) **Exemplars (the visual centerpiece):** uORF, CDS, lncRNA ORF from held-out Hepatocyte -- predicted
  vs observed per-nt P-site profile, with the 3-nt periodicity visible. Pick clean, high-period_OBS cases.
- (c) Depth crossover (Task 18): predicted beats measured Ribo-seq below ~29M P-sites -- the "why not just
  do Ribo-seq" punchline. Strong practical argument; belongs in a main figure.
- (d) small panel: one-hot ~= token-emb (drop-in F1 onehot 0.923 vs rinalmo 0.915 vs orthrus 0.916) --
  framed as "the simplest, dependency-free encoding is as good" (accessibility, not a horse race).

### Main Fig 2 -- immunopeptidome application (utility)
User proposed: scoring metrics vs other held-out datasets + cross-species; peptidomics benefit in the
MS2Rescore case, model vs expressed-ORF-universe null (ignore PRICE); highlight 2-3 high-confidence novel PSMs.

Proposed panels:
- (a) Cross-dataset / cross-species generalization: drop-in F1 and localization AUROC on held-out human
  (Ruiz-Orera, F1 0.931) and mouse (Wang, F1 0.929) -- the model transfers.
- (b) **Peptidomics benefit (the headline):** model-selected DB vs 3-frame null DB (both on the
  cell-type-specific EXPRESSED universe), MS2Rescore-rescored novel peptides at 1% class-specific FDR --
  **across MULTIPLE immunopeptidome datasets** (see expansion below), not just HBL-1. Drop PRICE per user.
  Show discovery-rate ratio model/null per dataset (forest-plot style) to demonstrate a reproducible pattern.
- (c) **2-3 highlighted novel PSMs:** annotated spectra with ms2pip mirror plots (predicted vs observed
  fragments) + the ORF class and genomic context. These are the concrete "new antigen candidates."

### Supp Fig 1-2 -- architecture x encoding grid + resources
- Metric grid: {transformer, Mamba} x {one-hot, token-emb (RiNALMo, Orthrus)} on pc/uORF/dORF/lncRNA
  profile Pearson + periodicity + localization AUROC.
- Resource comparison: trainable params, peak train memory, train time/epoch, inference time/1k tx, CPU vs
  GPU. Message: all within a narrow score band; one-hot transformer is the accessible default; Mamba is the
  O(L) option for very long transcripts.

### Supp Fig 3 -- more held-out HUMAN exemplars (Ruiz-Orera)
Additional clean predicted-vs-observed ORF concordance examples.

### Supp Fig 4 -- more held-out MOUSE exemplars (Wang)
Same, cross-species.

---

## KEY WEAKNESS TO FIX: Fig 2 statistical power

Current immunopeptidome evidence = 1 dataset (HBL-1), small counts (model 26 / null 4 rescored, 2 reps).
A brief communication headline needs a REPRODUCIBLE PATTERN. Two moves:

1. **Expand MHC-I to 3-5 datasets** (below) -- turn the one-off into a pattern (Fig 2b forest plot).
2. **Add a coding-potential baseline** (CPAT / CPC2) to the ORF-selection comparison -- "model beats
   CPAT-selected ORFs" is a stronger claim than "model beats naive 3-frame null." CPAT/CPC2 outputs already
   exist in the sibling project (results/cpat_pc, cpc2_pc) and are wired into the 11-scorer matrix.

## MHC-I dataset expansion -- candidate datasets (VERIFY open-access + matched RNA-seq before download; log in DATA_PROVENANCE.md)

Requirement for each: open MHC-I immunopeptidome raw MS + open RNA-seq for the SAME cell line (RNA-seq drives
the model's expressed universe + coverage input). Ideally spread across cancers/labs.

| candidate | source | why | notes |
|---|---|---|---|
| Ruiz-Cuevas DoHH2 + SU-DHL-4 | PXD020620 + PRJNA647736 (same study as HBL-1) | +2 DLBCL lines, immediate | within-study reproducibility; same protocol; RNA-seq present |
| Ouspenskaia 2022 (B721.221) | MSV000084787 + matched RNA/Ribo | cross-study lymphoblastoid; has RNA + Ribo + mono-allelic immunopeptidome | Ribo-seq could give ORTHOGONAL validation IF quality is good (unlike HBL-1) |
| THP-1 AML immunopeptidogenomics | published HLA-A*02:01 immunopeptidome + triplicate RNA-seq | cross-cancer (AML), different lab | verify raw MS + RNA-seq open |
| melanoma / other solid tumor line | e.g. via CrypticProteinDB's 26-dataset list | cross-cancer breadth | many are EGA/controlled (Chong 2020 = EGA, reject); pick an open one |

Priority: (1) DoHH2 + SU-DHL-4 (fast, locks reproducibility), (2) Ouspenskaia B721.221 (cross-study +
possible orthogonal Ribo-seq), (3) one cross-cancer open line (AML/melanoma). 3-4 total datasets makes a
credible pattern for a brief communication.

## What else would strengthen the approach (my suggestions)

1. **No-RNA-seq ablation.** Model with vs without the RNA-seq coverage channel -> shows the model uses
   expression, which is WHY it is cell-type-specific. Directly supports the "cell-type-specific" claim.
2. **CPAT/CPC2 baseline** in ORF selection (above) -- coding-potential is the obvious comparator a reviewer
   will ask for.
3. **FDR rigor panel (supp):** class-specific vs global FDR, and an entrapment/decoy calibration curve, to
   pre-empt the "are these real?" question -- this project already learned the global-vs-class-specific FDR
   lesson; show it.
4. **MS2Rescore substrate-dependence as a rigor point:** tryptic whole-proteome (no lift, class too sparse)
   vs immunopeptidome (large lift) -- demonstrates the team understands when rescoring helps. Small supp
   panel; strengthens credibility.
5. **Accessibility framing throughout:** CPU-runnable, ~5M params, no Ribo-seq, no FM dependency (one-hot).
   This is the practical selling point of a brief communication and should be explicit in Fig 1 + abstract.
6. **Orthogonal support for highlighted PSMs:** ms2pip mirror plots are the minimum bar; if any highlighted
   novel ORF has cross-species conservation or independent evidence, note it. Synthetic-peptide validation
   is the gold standard but likely out of scope for a brief communication.
7. **Cross-species immunopeptidome (stretch):** if an open mouse MHC-I dataset with RNA-seq exists, one mouse
   immunopeptidome point would extend the pattern beyond human -- big credibility boost, aligns with the
   cross-species theme of Fig 2.

## Decisions (LOCKED 2026-07-20)

- ~~D1: **Main Fig 1 = deployed TRANSFORMER `orf_v2_attn`, one-hot.** Mamba mixer -> Supp 1-2.~~
  **SUPERSEDED 2026-07-31 by D1b** (original text kept above for provenance).
- D1b (user decision 2026-07-31, replaces D1): **Main Fig 1 = `orf_v2_mamba4`, one-hot.** It wins the
  3-seed union comparison on held-out profile Pearson with non-overlapping seed ranges (mean 0.6799 vs
  0.6595; mamba4's worst seed 0.6753 beats attn's best 0.6603). `orf_v2_attn` moves to supplemental
  but stays maintained and shipped, because it is the only one of the two that runs without a GPU
  (mamba-ssm needs CUDA kernels) and it is therefore the released CPU inference path.

  **State the scope of mamba4's win honestly in the caption.** The +0.0204 Pearson edge is real on
  per-nucleotide profile shape and does NOT carry into any downstream task measured so far:
  B721.221 drop-in ORF calling is a tie (F1 0.665 mamba4 vs 0.666 attn, and the same tie in
  precision/recall trade), the 4-dataset MS panel is a tie (49 vs 48 novel peptides, with the sign of
  the difference flipping by dataset), and the 12-population macrophage run is a tie. Where mamba4
  does lead on discovery (BMDM, 40 vs 33 novel peptides) it is by calling MORE ORFs (1,533 vs 1,244
  sequences) at indistinguishable per-sequence density (26.1 vs 26.5 per 1,000), not by better
  per-ORF discrimination. Also note mamba4 is ~5x noisier across seeds (spread 0.0098 vs 0.0018), so
  a single-seed mamba4 number must never be quoted without the spread.
- D2: **MHC-I expansion = HBL-1 (have) + DoHH2 + SU-DHL-4 + Ouspenskaia B721.221 + THP-1 AML + one mouse**
  (6 total). Start order: DoHH2/SU-DHL-4, then B721.221, THP-1, mouse. 30-min data-access cap per dataset.
- D3: **Include CPAT/CPC2 coding-potential baseline** in the ORF-selection comparison (A549 + immunopeptidomes).
  **DONE 2026-08-06, all 4 datasets, figure `figures/D12_cpat_cpc2/`.** Both arms enumerate from the
  identical candidate pool as `null_atg` (regression-tested), so only the selection rule differs. The
  result splits by assay: CPAT/CPC2 win the tryptic whole proteome (25 / 23 novel peptides vs the
  model's 11), the model wins all three HLA-I immunopeptidomes (14 vs 10, 22 vs 14, 32 vs 8). On
  discovery efficiency the model's Poisson arm leads all four. Fig 2 must present the split, not just
  the immunopeptidome half of it.
- D4: **Depth-crossover in Main Fig 1.** PRICE dropped from Fig 2 (may keep a one-line caveat in supp).
- D5 (new): **Add Ribo-seq quality metagene plots (ribotish) as supplemental figures** for every study used
  that has Ribo-seq (HBL-1, Ouspenskaia B721.221 if acquired, Chothani training tissues for reference).
  **DONE 2026-08-07 as `figures/S_riboseq_qc/`, built from RiboCode `metaplots` output rather than
  ribotish.** ribotish needs a GENOME BAM plus a GTF; this project keeps transcriptome BAMs for
  everything except the training set, so honouring the letter of D5 would have meant re-aligning a
  dozen studies to regenerate quantities already on disk. `metaplots` had already been run on every
  dataset and writes the same three axes (per-read-length P-site offset, frame sums at annotated
  start codons, periodicity test), so reading those files gave a LARGER panel (all 16 library groups)
  at zero compute. See that folder's FIGURE_DATA_INPUTS.md.

## Execution status (task list #17-24)

- #17 update this doc (DONE)
- #18 DoHH2 + SU-DHL-4 pipeline
- #19 Ouspenskaia B721.221 pipeline
- #20 THP-1 AML pipeline
- #21 mouse immunopeptidome pipeline
- #22 CPAT/CPC2 baseline
- #23 ribotish metagene QC supp figs
- #24 consolidate multi-dataset model/null (Fig 2b) + pick highlight PSMs (Fig 2c)

================================================================================
AMENDMENT 2026-08-08 -- numbers above are pre-union; panel structure is unchanged
================================================================================

APPENDED, NOT EDITED. Nothing above this line was altered; the pre-amendment file is preserved at
`FIGURES_PLAN.md.bak.2026-08-08`. Read this section before quoting any figure from the plan above.

### Why

Every held-out drop-in number in this plan was measured on `orf_v2_attn_onehot_holdout_Hepatocytes`,
the pre-nokozak / pre-mm1 / pre-union checkpoint, not on the shipping union models. Run directories are
named for the DATASET rather than the CHECKPOINT, so the staleness was invisible from the path -- it
lives only in the npz `meta` field. It was found by scanning every figure generator, which showed five
of the six Fig 1 panels in the same state.

### Corrected numbers (released models, attn / mamba4)

| where the plan says | actual, on the released models |
|---|---|
| "drop-in F1 ~0.92" (Fig 1a) | Hepatocytes standalone F1 **0.909** (attn) |
| "Ruiz-Orera, F1 0.931" (Fig 2a) | shape (`pred_obsdepth`) 0.929 / 0.934 holds; **standalone `pred_preddepth` 0.876 / 0.880** |
| "Wang, F1 0.929" (Fig 2a) | shape 0.923 / 0.927 holds; **standalone 0.867 / 0.867** |
| "onehot 0.923 vs rinalmo 0.915 vs orthrus 0.916" (Fig 1d) | pre-union values; the one-hot >= FM conclusion is a within-checkpoint comparison and is unaffected, but the three absolute numbers must not appear beside released-model numbers |

The pattern: the SHAPE claim survived intact everywhere; the STANDALONE claim (predicted shape AND
predicted depth) dropped 0.05 because the count head over-calls on the 84,472-tx union universe. Any
standalone F1 in a caption must now be stated with its calling arm (theta=1 or Poisson theta*).

### Panel structure: unchanged, with three additions

The proposed Fig 1 (a-d) and Fig 2 (a-c) structure still stands. Additions from work since the plan:

- **New candidate panel, B6 `figures/B6_input_ablation/`** -- the input ablation on the deployed
  recipe. It is the direct evidence for the cell-type-specificity claim that Fig 2 rests on, and it
  NARROWS that claim: sequence-only recovers 98.5% of profile shape, so the specificity lives in the
  count head (-0.193 count Pearson without RNA-seq), not in the shape. Wherever the manuscript says
  the model is cell-type-specific, it must mean "which ORFs clear the depth threshold", and cite
  -0.193. Do not write "cell-type-specific profile shape" -- it is measurably close to false.
- **Both calling arms, always** (feedback_two_arm_orf_calling). `results/released_two_arm_orf_calls.json`
  has all 16 rows. Poisson raises precision on all 8 dumps but improves F1 on only ONE of four datasets
  (Wang) and halves non-canonical F1 every time, so a figure that shows only the Poisson arm overstates
  the method for discovery use and only the theta=1 arm overstates it for precision use.
- **Fig 1a caveat.** A3 (LOTO spread) is the one panel still on the pre-union recipe -- there is no
  union 9-fold and building one is nine retrains. Its CLAIM (that the fold spread tracks held-out
  target periodicity at r=0.81, rather than being a generalization gradient) is a property of the
  design and survives; its absolute Pearson values are checkpoint-dependent. Label the panel.

### Status of the two gates named in the plan

- Fig 2b power fix: DONE. `figures/F2b_discovery_forest/`, 5 datasets x 2 models x 2 arms = 20 points,
  all above 1.0.
- Fig 2c (D14 highlighted PSMs): still the only MS-gated item. B721.221 MS is blocked (MassIVE FTP
  unreachable from prism); the immunopeptidome panel is at 4 of the 6 datasets decision D2 asked for.
