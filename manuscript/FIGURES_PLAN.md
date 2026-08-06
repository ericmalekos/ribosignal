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

- D1: **Main Fig 1 = deployed TRANSFORMER `orf_v2_attn`, one-hot.** Mamba mixer -> Supp 1-2.
- D2: **MHC-I expansion = HBL-1 (have) + DoHH2 + SU-DHL-4 + Ouspenskaia B721.221 + THP-1 AML + one mouse**
  (6 total). Start order: DoHH2/SU-DHL-4, then B721.221, THP-1, mouse. 30-min data-access cap per dataset.
- D3: **Include CPAT/CPC2 coding-potential baseline** in the ORF-selection comparison (A549 + immunopeptidomes).
- D4: **Depth-crossover in Main Fig 1.** PRICE dropped from Fig 2 (may keep a one-line caveat in supp).
- D5 (new): **Add Ribo-seq quality metagene plots (ribotish) as supplemental figures** for every study used
  that has Ribo-seq (HBL-1, Ouspenskaia B721.221 if acquired, Chothani training tissues for reference).

## Execution status (task list #17-24)

- #17 update this doc (DONE)
- #18 DoHH2 + SU-DHL-4 pipeline
- #19 Ouspenskaia B721.221 pipeline
- #20 THP-1 AML pipeline
- #21 mouse immunopeptidome pipeline
- #22 CPAT/CPC2 baseline
- #23 ribotish metagene QC supp figs
- #24 consolidate multi-dataset model/null (Fig 2b) + pick highlight PSMs (Fig 2c)
