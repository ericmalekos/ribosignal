# What is CURRENT, what is SUPERSEDED, what is IN PROGRESS

Last reviewed 2026-08-15. One page, so that "which numbers are live?" never needs archaeology.

The hazard this exists to prevent is real and has happened here more than once: a superseded number
surviving in a README for a week after the data contradicted it (P11's BMDM-only headline), and a QC
panel silently rendering with no training data after its source directory was replaced.

---

## 1. Models

| | path | status |
|---|---|---|
| **DEPLOYED** | `results/loto/orf_v2_{attn,mamba4}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` | **CURRENT.** Every figure, the poster and the manuscript use these. |
| final-recipe retrains | `results/loto_canon/` (3 runs) | **TRAINED, NOT DEPLOYED, NOT ORF-CALL-EVALUATED.** Task #86. See that dir's README. |
| 9-fold LOTO, rinalmo folds | `results/loto_9fold/`, `results/fold*_rinalmo*` | historical, superseded by the union recipe |

**Do not swap a `loto_canon` checkpoint into a figure.** SETTLED 2026-08-15: the ORF-call
evaluation the standing rule demands was run (5 arms on one substrate) and agrees with the
profile metrics -- non-canonical F1 moves +0.0034 (attn) and -0.0003 (mamba4), both inside the
+-0.005 materiality band. Final-recipe training changes nothing; the deployed checkpoints are
correct. See `results/loto_canon/orfcall_comparison/`.

---

## 2. Ribo-seq / ORF-call results

| topic | CURRENT | SUPERSEDED |
|---|---|---|
| mouse-liver 3x3 | `results/mouse_liver_3x3_canon/` | `results/mouse_liver_3x3/` (emptied, README) + `results/_archive_offrecipe_2026_08_13/` |
| saliency | `results/aggregate_saliency_canon/` | `results/aggregate_saliency/` (README) |
| codon occupancy | `results/codon_occupancy_canon/` | `results/codon_occupancy/` (README) |
| ORF-track channel ablation | `results/orf_channel_ablation_canon/` | `results/orf_channel_ablation/` (README, **3.2 GB**, first deletion candidate) |
| human ORF calls (THP-1, CAR-T) | `results/human_orf_calls/` | -- but **CAR-T is provisional**, see section 4 |
| merged 28-library liver reference | `results/merged_liver_ribocode/` + `overcall/` | new 2026-08-14, no predecessor |
| reproducibility, 8 tissues | `results/chothani_regeneration/_canonical_<T>/` | the source study's `ribocode_per_tissue/` no longer exists |

**The final-recipe alignment recipe is the dividing line.** Anything computed before it was aligned
without `--alignEndsType EndToEnd` and pooled without the ncRNA + cross-gene filter. Off-recipe
non-canonical F1 was INFLATED (~22% more novel calls than the corrected pipeline produces), so
off-recipe and final-recipe numbers must never be compared. See `docs/PIPELINE_POLICY.md`.

---

## 3. Proteogenomics

| | path | status |
|---|---|---|
| **CURRENT, all figures** | `.../pgx_xsubtype_aa30/` | **30 AA + ATTENTION model** (2026-08-15). P11 and P13 now read this; P12 regeneration pending. 7 arms searched per population against identical spectra, including BOTH `model_predicted` (mamba) and `model_predicted_attn`. |
| SUPERSEDED | `.../pgx_xsubtype/` | 7 AA, mamba only. Violates Rule 5 for tryptic data. Its DB-size and density ratios (156x/111x) were INFLATED by sub-30-aa null sequences; correct values are 78x/43x. Reproducible via P11's `--source`/`--arm`. |
| A549 30-AA rebuild | **ARCHIVED 2026-08-15, will NOT happen** | A549 is frozen at 7 AA and is the project's ONLY Rule 5 exception. See `proteogenomics/data/_archive_a549_2026_08_15/README.md`. Its 5 consuming figures now carry the exception note. **It stays in every panel**: it is the only tryptic human dataset and carries the results least flattering to the model (CPAT/CPC2 win the tryptic proteome; the one Poisson failure, -35; the largest global-vs-class FDR gap, 37.5x). |
| human datasets | `proteogenomics/data/{A549,HBL1,SUDHL4,DoHH2,THP1,B721}_pilot/pgx_*/` | **CURRENT.** HLA-I ones (HBL1, SUDHL4, DoHH2, THP1, B721 -- all `nonspecific`/`termini=0`) stay at 7 AA permanently and correctly. A549 is tryptic and archived at 7 AA. |
| unfiltered macrophage | `.../macrophage_tissue/superseded_unfiltered/` | superseded (pre-FDR-filter era, Task 52) |
| frozen vs unfrozen reports | `pgx_xsubtype/frozen_reports/` vs `pgx_xsubtype/report_*.json` | **BOTH are current, and they differ.** Figures use FROZEN. results.md Task 53 and the tutorial quote UNFROZEN. Never mix: BMDM's delta-total is +31 frozen and -21 unfrozen for identical databases and spectra. |

**Assay determines the ORF floor** (Rule 5): tryptic whole-cell lysate 30 AA, MHC/HLA 7 AA. The 12
macrophage populations and A549 are tryptic; HBL-1, SU-DHL-4, DoHH2, THP-1, B721.221 are MHC.
Classify from `fragger.params`, never from a name.

---

## 4. In progress / known-provisional

| item | state | task |
|---|---|---|
| **P5's four held-out arms** | **2 of 4 are ON-RECIPE.** THP-1 dumped + rescored 2026-08-15 (F1 0.851 -> 0.831, n_ref 14,304 -> 13,225); CAR-T on-recipe since 08-14 (0.904 -> 0.886, 8,856 -> 8,253). **hepatocyte and iPSC-CM are OFF-RECIPE** and likely OVERSTATED -- P5 now prints every arm's substrate and warns. The drops are the REFERENCE shrinking, not the model: 2x2 decomposition gives reference -0.0740 / prediction -0.0001. THP-1 has NO UMIs (filter only, ~9.7%); CAR-T had ~29% PCR duplicates on top. | #92 |
| final-recipe retrain ORF-call eval | **SETTLED 2026-08-15**: no material difference (attn +0.0034, mamba4 -0.0003 non-canonical F1, both inside +-0.005). Deployed checkpoints stay deployed. | #86 done |
| matched-depth 3x3 arm | not started | #75 |
| tutorial push | blocked on explicit instruction | #83 |

**Poster package REFRESHED 2026-08-15** (`poster_package/`, `poster_package.tar.gz`, 21.3 MB,
md5 `44e8bef8e736f430c87ba4f4a6a89791`). Three structural fixes, not just new numbers:

- Its README is now `docs/POSTER_PACKAGE_README.md` and is COPIED in. The builder opens with
  `shutil.rmtree(out)`, so the previous hand-written 18 KB README inside `poster_package/` was one
  rebuild away from being deleted with no trace.
- Ships `DATASET_LABELS.md` (generated from `data/dataset_registry.tsv`) + the registry itself, so
  internal keys in filenames map to display labels for captions.
- **P5 now renders both architectures.** `--model attn` writes `P5_heldout_human_attn.*`; the mamba4
  default (locked decision D1b) is untouched. The old README printed P5's mamba4 numbers under an
  "attn" heading -- every panel's model is now carried in its own `*_values.json`.

**Complete as of 2026-08-15:** macrophage attn sweep at 30 AA (#91), **P11 + P12 + P13 all rebuilt on
it**, Brain acquired + aligned + metaplots (9th Chothani tissue), GSE39561 verified negative control,
THP-1 final-recipe pack, dataset registry, Chothani samplesheet reconciled to 79 runs / 9 tissues.

> **Two entries in this table were WRONG when written, and the correction is the point of the page.**
> It claimed "P11 and P13 done, P12 pending". P13's *underlying JSON* had been recomputed on the
> attention tree but the *figure* still read the 7-aa mamba JSON, so two panels were stale, not one.
> A stale duplicate of this section also survived below it asserting CAR-T was PROVISIONAL with a
> queued job, four lines under the entry saying it is CANONICAL. Both are fixed. Check the artifact's
> own recorded source (`*_values.json` -> `source`), never the note that says it was regenerated.

---

## 5. Figures

39 folders under `figures/`, **all with a `FIGURE_DATA_INPUTS.md`** (verified 2026-08-14).
`figures/README.md` is the index; `docs/POSTER_LAYOUT.md` is generated by
`figures/make_poster_layout.py` and reads every headline number live from each panel's
`*_values.json`, so it cannot drift from the figures.

**Superseded figure content, corrected in place:**

- **P11** was BMDM-only and claimed the model is "the only one that increases TOTAL unique peptides".
  BMDM (+31) is the best of 12; the median is -8 and the model is above baseline in 5/12. Rebuilt
  across all 12 populations.
- **S_riboseq_qc** silently rendered 7 arms with **no training tissue** after its source directory
  was replaced by the recipe standardisation, while claiming 16. Repointed; now 15 arms, 8 training.
- **P9's** "model only" bar was labelled "candidate FPs"; it is now the measured ~92%.
- **D13/D12/F2b/P11** gained absolute total-unique-peptide panels (standing rule: report totals even
  when they disfavour the model).

**Claims that were measured and found FALSE** -- these are in the poster layout's do-not-print list:

- "two architectures agreeing is evidence" (8.7%, indistinguishable from one model)
- "the predicted database beats Ribo-seq" (it does not, on 11 transfer populations)
- P11's BMDM-only headline
- any non-canonical F1 predating the recipe standardisation

---

## 6. Housekeeping

Disk, largest superseded item first (group ceph quota is shared and near-full,
memory `feedback_group_ceph_15tb_quota`):

- `results/orf_channel_ablation/` -- 3.2 GB, nothing depends on it, safe to delete
- `results/_archive_offrecipe_2026_08_13/` -- retained deliberately; the canonical-vs-off-recipe
  delta is itself a documented finding
