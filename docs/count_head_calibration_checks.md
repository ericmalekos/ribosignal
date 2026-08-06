# Count-head calibration -- task tracker

Living document: each check carries a **Status** (PLANNED / IN-PROGRESS / DONE) and a **Result**
block filled in as we go. Progress log at the bottom.

## CURRENT STATE (snapshot 2026-07-24)

The final model (`results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes`) is trained and
localization-equivalent to the deployed model; the open work is the de novo ORF-calling calibration/validation
(O1/O2/O3), most of which is CPU-done, with the honest cross-dataset numbers gated on a GPU queue.

**Landed (CPU, in hand):**
- **O1 dial** -- CDS-anchored stringency dial + two-arm protocol (standard theta=1 + downsample/Poisson).
  Poisson injection is the winning lever (+0.10 novel precision). Checks 2/3/4 DONE on Hepatocytes.
- **O2 ceiling** -- the de novo bar computed: between-experiment F1(Wang_obs, Janich_obs) = 0.675/0.735 overall,
  **0.505 for non-canonical** (`o2_between_dataset.py`). Model-vs-ceiling still pending (needs the model dumps).
- **O3 track 1 (FP-filter)** -- DONE in-sample: shape features + per-class logreg (`orf_features.py`,
  `fp_filter.py`, `fp_filter_hep.joblib`). +0.10 novel / +0.12 uORF precision at FIXED CDS recall 0.956;
  AUC 0.669; ceiling-capped by model smoothness. Phase-1 refuted the uniformity premise.

**Running / queued (GPU, order shown; GPUs fully booked):**
```
dump_gse_nk 35889920    DONE + dropin 35889921 DONE  ->  Check 5 no-Kozak GSE landed (see below)
dump_wang_nk 35889941   DONE + dropin 35889942 DONE  ->  O2 model-vs-ceiling landed: CDS 0.84 of ceiling de novo
dump_janich_nk 35890103 PENDING (Priority)    ->  O2 symmetric check (Janich model vs ceiling)
loto_peaky 35890525     PENDING (Priority)    ->  O3 track 2 A/B (anti-smoothing model, PK=0.2)
dropin_janich_nk        PENDING (Dependency, fires after its dump)
```

**Landed: O2 model-vs-ceiling, BOTH arms (Wang) -- Check 7.** Standard arm: CDS de novo 0.621 (0.84 of the 0.735
ceiling). **Poisson arm: CDS de novo 0.744 at theta=0.02 -- PASSES the ceiling** (calibrated sequence+RNA-seq
predicts liver CDS as well as a 2nd real Ribo-seq experiment). The standard-arm gap was pure CDS over-calling
(21,315 vs ~14,700 real). Non-canonical: 0.067 -> 0.195 under Poisson but still << its 0.505 ceiling (that gap =
FP-filter + peaky model). See Check 7. Janich symmetric check + peaky training still queued.

**Also in flight (CPU, off the GPU queue): macrophage test data download** (Nat Commun 2022 s41467-022-35095-7):
IPX0001245001 proteomics (216 .raw, ~173 GB, iProX, ~12h) + PRJNA482293 RNA-seq (36 single-end fastq, 175 GB,
ENA md5-verified, ~4h) -> `proteogenomics/data/macrophage_tissue/`. Matched proteome+transcriptome for 12
tissue-macrophage populations, for the proteogenomics pipeline (model-informed DB search). See DATA_PROVENANCE.md.

**Landed since snapshot: Check 5 no-Kozak GSE (2026-07-24).** The final no-Kozak mm1 model reproduces the
over-calling pattern and is SLIGHTLY WORSE on non-canonical than the deployed heuristic model (NT novel 73x vs
64x, dORF 76x vs 46x over-call; predF1 0.655 vs 0.683) -- dropping the Kozak track (mandated for profile
generalization) removes its suppression of weak-context non-canonical starts. Deep-reference validation: extra
NT calls 11.8% real overall, novel 2.7% / dORF 2.4% (spurious) vs the 87.3% ceiling. Fix = O3, not re-adding
Kozak. Standard arm only; Poisson arm + FP-filter on GSE = follow-up.
`results/gse120762/{compare_deployed_vs_nokozak_nt_lps,overcall_check_nokozak}.txt`.

**Pending actions (fire when the matching GPU job lands):**
1. ~~dropin_gse_nk~~ DONE (above). Remaining GSE follow-up: Poisson arm + FP-filter (need the GSE no-Kozak
   `pred_profiles.npz`).
2. **dump/dropin_wang_nk (running) + _janich_nk** -> (a) O2 model-vs-ceiling: is F1(model, datasetB_obs) ~ the
   0.505 non-canonical ceiling? (b) the HONEST cross-dataset FP-filter: train on Hepatocytes, apply to
   Wang/Janich (in-sample +0.10 is an optimistic upper bound).
3. **loto_peaky 35890525** -> the O3 A/B: drop-in -> `orf_features.py` -> `fp_filter.py train` -> CDS-anchored
   table. Wins if feature-AUC clears 0.669 AND novel over-calling drops at matched CDS recall AND val_pearson
   holds. If val tanks, dial PK to 0.05-0.1 and re-launch.

**Standing rules in force:** two-arm protocol (standard + Poisson) on every ORF-calling test; always evaluate
by how ORF CALLS change (RiboCode drop-in real/pred_obsdepth/pred_preddepth), never a single metric; validate
precision against a DEEP independent reference, never a shallow same-dataset one; KOZAK never set (`--kozak
none`) on all retrains. Nothing committed to git (awaiting explicit instruction).

## Goal and strategic framing

The model predicts per-nt Ribo-seq from sequence + RNA-seq. Two heads: a **profile** (shape) and a
**count head** (per-tx total). The drop-in split shows the shape is good (`pred_obsdepth` precision
~0.85) but the standalone count (`pred_preddepth`) **over-calls**, concentrated in short non-canonical
ORFs, because RiboCode gates on absolute P-site thresholds and the count head's absolute scale is
tied to training depth, not the target's.

**The tool's value is calling ORFs when Ribo-seq is absent or low quality.** So we prioritize
**high-depth observed cases** (where truth is measurable), NOT shallow matching, and expose a
**CDS-anchored stringency dial**: CDS calls are trusted, so the user picks an operating point as a
fraction of CDS recovery ("the stringency that recovers 90% of true CDS") and reads off which
non-canonical ORFs come along at that confidence. 90%-CDS is stricter than a 110%-of-CDS-count point.

## Objectives

- **O1 -- CDS-anchored stringency dial.** A calibrated knob letting a user call ORFs at a chosen
  CDS-confidence level (sensitivity vs specificity), with a confidence estimate per non-canonical ORF.
- **O2 (headline validation) -- match the between-experiment reproducibility ceiling.** Two
  independent Ribo-seq datasets of the SAME tissue+species agree only up to biological + technical
  reproducibility. If the model (sequence + RNA-seq, NO Ribo) predicts dataset B's ORF calls about as
  well as dataset A's real Ribo-seq does, the model is "as informative as a second real experiment" ->
  strong evidence for de novo ORF calling. Benchmark on mouse liver (Wang vs Janich), at CDS-anchored
  operating points, per ORF class.
- **O3 (forward goal) -- a more robust novel-ORF caller, beyond RiboCode's frame test.** RiboCode
  calls ORFs with a SINGLE univariate periodicity test (Wilcoxon / Mann-Whitney of frame-0 vs
  frame-1/2 P-sites per codon). That one feature is exactly what the model's smooth predicted density
  fools (Check 4's over-calling: even spurious short ORFs look cleanly periodic). Poisson injection
  patches it, but a MULTIVARIATE / ML caller that combines periodicity with other discriminative
  features (coverage uniformity/entropy, ORF length, start-codon context, RNA-seq support, 3'-stop
  drop-off) should filter false positives far better. Two tracks: (A) co-opt an existing ML-based
  caller (RibORF = SVM on periodicity + %-max-entropy + uniformity; RP-BP = Bayesian periodicity;
  DeepRibo = CNN on signal+sequence; also Ribo-TISH / ORFquant / PRICE); (B) build a lightweight
  per-ORF false-positive-filter classifier. Key enabler: O2's cross-dataset reproducibility +
  the deep-observed validation give a PRINCIPLED real/false label with no extra experiments, so
  FP-filtering becomes a supervised problem. Evaluated on the SAME CDS-anchored dial + two-arm +
  over-call framework; success = higher novel precision at matched CDS recall than RiboCode+Poisson.
  Ideally BAM-free (drives off the model's predicted density, like the RiboCode drop-in).

Every check is scored by the standing rule: RiboCode drop-in + over-call-vs-deeper-reference
validation ([[feedback_orf_call_eval_standing]]).

## Evaluation protocol (STANDING, from 2026-07-24): always run BOTH arms

Every ORF-calling test / dataset going forward is run through TWO calling arms and reported
side by side, so the calibration benefit is always visible:
- **Standard arm** -- `pred_preddepth` at `--pred_scale 1.0`, deterministic (raw model rate ->
  round -> RiboCode). The naive standalone prediction.
- **Downsample+Poisson arm** -- `--pred_poisson` at the CDS-anchored `--pred_scale theta` (the
  recommended de novo recipe: simulate a Ribo-seq experiment of effective depth theta by drawing a
  Poisson sample from the predicted rate). theta is the dataset's CDS-anchored operating point
  (e.g. ~90% CDS recall).
Check 4 established Poisson as the winning lever (+0.10 novel precision at matched CDS recall on
Hepatocytes); this protocol makes the comparison mandatory on every subsequent dataset.

## High-depth calibration / validation set

| dataset | species | status | role |
|---|---|---|---|
| **Hepatocytes** (Chothani holdout) | human | READY (real 15,441 calls + mm1-nokozak preds + drop-in) | primary calibration anchor |
| **Wang liver** | mouse | pack ready; need nokozak track + dump + drop-in | O2 dataset A (model-run) |
| **Janich pooled liver** (5-6 samp) | mouse | observed calls exist; pack TBD via `prepare/` | O2 dataset B (+ deeper anchor) |

Shallow sets (GSE120762 NT) are a stress test only, never a calibration target.

---

## Check 0 -- assemble the high-depth set    **Status: IN-PROGRESS**
Method: Hepatocytes done. Wang: nokozak track -> `dump_pred_profiles` -> `ribocode_dropin` (mouse annot,
3 variants). Janich: `prepare_pack.py` from psites+coverage -> track -> dump -> drop-in. Confirm Wang and
Janich share the mouse annotation so ORF_IDs are comparable (required for O2).
**Result:** _(pending)_

## Check 1 -- high-depth over-call baseline    **Status: PARTIAL**
Method: `compare_dropin_models.py` per dataset -- precision/recall/F1 + over-call ratio by ORF_type,
deployed vs mm1-nokozak.
**Result:** Hepatocytes: deployed `pred_preddepth` P=0.741, nokozak P=0.713; over-call uORF 1.6->2.0x,
novel 2.1->2.4x, dORF 2.8->3.4x (no-kozak slightly worse). Wang/Janich pending.

## Check 2 -- CDS-anchored operating curve (core)    **Status: IN-PROGRESS**
Method: stringency theta = scale on the PREDICTED effective depth (`ribocode_dropin.py --pred_scale`,
theta<1 -> fewer calls). Sweep theta; at each, vs observed calls compute CDS recall, CDS precision,
n_calls/n_real_CDS. Deliverable: `cds_recall(theta)` (scripts: dropin_sweep_hepatocytes.sbatch +
cds_anchored_curve.py).
**Result:** BUG CAUGHT + FIXED -- the first sweep used `--subsample`, which only thins the REAL variant
(no-op for pred_preddepth; all 10 thetas gave identical 19,110 calls). Added `--pred_scale` (scales the
predicted total for pred_preddepth) + `--pred_poisson` (noise injection, Check 4b). **DONE (Hepatocytes,
mm1-nokozak):** the dial is monotonic and CDS is a stable anchor (CDS precision ~0.91-0.94 across the
whole range). Curve: results/loto/.../dropin_sweep/{cds_anchored_curve.tsv,dial_summary.txt}.

## Check 3 -- non-canonical yield per CDS operating point    **Status: DONE (Hepatocytes)**
**Result (mm1-nokozak, Hepatocytes; real CDS=10,947; precision vs observed deep calls):**
```
CDS_recall  theta  CDS_prec  calls/CDS    uORF(prec)    novel(prec)   dORF(prec)
   0.322    0.02    0.906      0.37       93(0.82)      24(0.71)      4(0.00)
   0.733    0.10    0.929      0.91      922(0.72)     175(0.59)     64(0.28)
   0.863    0.20    0.934      1.18     1913(0.65)     400(0.54)    137(0.25)
   0.930    0.35    0.936      1.39     2895(0.58)     691(0.45)    249(0.18)
   0.954    1.00    0.936      1.75     5091(0.43)    1271(0.34)    733(0.11)
   0.955    2.50    0.936      2.07     6744(0.35)    1870(0.25)   1775(0.06)
```
Tightening (lower theta) trades novel yield for novel precision monotonically: uORF precision 0.35->0.82,
novel 0.25->0.71. This IS the dial. Wang/Janich cross-dataset validation of these operating points = Check 5.

## Check 4 -- which lever wins    **Status: DONE (Hepatocytes)**
Method: novel/uORF precision at matched CDS recall (compare_levers.py). (a) depth-scale (--pred_scale),
(b) Poisson injection (--pred_poisson), (d) pval threshold; (c) isotonic ceiling = pred_obsdepth.
**Result: POISSON INJECTION WINS.** novel_prec / uORF_prec at matched CDS recall:
```
CDS_recall   depth        poisson       pval
   0.85     0.54/0.65    0.59/0.78    0.45/0.63
   0.90     0.49/0.61    0.59/0.75    0.45/0.69
```
- **(b) Poisson injection: best** -- +0.10 novel precision over depth-scale at 0.90 CDS recall. Sampling
  the smooth predicted density as Poisson(rate) gives it real-Ribo-like noise, so spurious weak short
  ORFs don't survive the frame test.
- (a) depth-scale: decent simple baseline.
- (d) pval threshold: worst -- removes calls uniformly, not preferentially the spurious ones.
- (c) isotonic/pred_obsdepth: NOT useful for novel precision (it sits at full real depth: CDS_recall 0.97
  but novel_prec only 0.31). Per-tx count-scale calibration does not reduce over-calling; the stringency
  reduction (depth/poisson) does.
**-> the dial should use --pred_poisson (+ --pred_scale for the operating point).**

## Check 5 -- calibration transfer    **Status: IN-PROGRESS (dial mechanism transfers)**
Method: apply the Poisson dial to unseen datasets; check theta->operating-point map + precision.
**Result (deployed-model GSE, CPU while GPUs jammed):** the KNOB TRANSFERS -- theta=0.05 -> ~90% CDS
recall on human Hepatocytes (0.925), GSE-NT mouse BMDM shallow (0.907) and GSE-LPS deep (0.866); CDS
precision stable ~0.87-0.93. BUT non-canonical PRECISION at an operating point is NOT a portable
number: novel prec at theta=0.05 = 0.58 (Hep, deep ref) / 0.19 (GSE-LPS, deep) / 0.07 (GSE-NT, SHALLOW
ref). This is a REFERENCE-DEPTH artifact -- precision is vs each dataset's own observed calls, and a
shallow reference misses real ORFs so the model's recovered ORFs count as false. => must evaluate vs a
DEEP INDEPENDENT reference (= the O2 test).
**Deep-reference follow-up (GSE-NT theta=0.05 vs deep GSE-LPS obs) -- the artifact is CLASS-SPECIFIC:**
CDS 0.89->0.83 (solid); uORF 0.23->0.46 (reference-depth artifact REAL: half the shallow "misses" are
corroborated by deep data); novel 0.07->0.11 vs a 0.54 baseline (real NT-novels in deep-LPS) = model
novels validate at only ~20% of the real-novel rate -> NOVEL is GENUINELY over-called, NOT an artifact.
So: knob transfers, CDS trustworthy, uORF was unfairly shallow-penalized, but NOVEL over-calling is a
real defect (deployed model) = exactly the O3 (robust caller) target.
**no-Kozak GSE redo -- DONE (2026-07-24, dump 35889920 / dropin 35889921 COMPLETED; standard arm)
`results/gse120762/{compare_deployed_vs_nokozak_nt_lps,overcall_check_nokozak}.txt`:**
- The over-calling PATTERN reproduces on the final no-Kozak mm1 model, and is SLIGHTLY WORSE on non-canonical
  than the deployed heuristic model (NT over-call ratios: novel 73x vs 64x, dORF 76x vs 46x, uORF 14x vs 12x;
  LPS similar). predF1 dips accordingly (NT 0.655 vs 0.683, LPS 0.716 vs 0.743) via lower precision + higher
  recall. Interpretation: the Kozak track had been suppressing weak-context non-canonical starts; dropping it
  (mandated for PROFILE generalization, Task 20) removes that suppression -> a bit MORE non-canonical over-
  calling. The fix is O3 (FP-filter + anti-smoothing), NOT re-adding Kozak (deployed stays heuristic-pinned).
- **Deep-reference validation (no-Kozak NT extras vs deep LPS obs) confirms the over-calling is real + class-
  specific:** of 9,774 extra NT calls, 11.8% are in deep LPS-real (vs the 87.3% reproducibility ceiling);
  by class novel 2.7%, dORF 2.4% (genuinely spurious), uORF 15.6%, Overlap_uORF 20.3%, annotated 22.4%. Same
  verdict as the deployed model: CDS + uORF defensible, novel/dORF genuinely over-called.
- **Two-arm gap:** this is the STANDARD arm only (pred_preddepth theta=1). The Poisson arm + the FP-filter
  applied to these GSE calls are the pending follow-ups (need the GSE no-Kozak `pred_profiles.npz`).

## Check 6 -- expose the dial    **Status: PLANNED**
Method: wire calibration into the calling path as `--cds-recall 0.90` (emit calibrated theta + per-ORF
confidence).
**Result:** _(pending)_

## Check 7 -- between-dataset reproducibility ceiling (O2)    **Status: model-vs-ceiling DONE (standard arm)**
Method: coord-key-matched F1(Wang_obs, Janich_obs) = the CEILING (two real liver Ribo experiments);
then model F1(Wang_model, Janich_obs). o2_between_dataset.py (model-optional).
**Result -- CEILING: F1(Wang_obs, Janich_obs) = ALL 0.675, CDS 0.735, non-canonical 0.505.** Two real
liver experiments agree only ~0.68 overall; non-canonical ORFs are just ~0.50 reproducible even between
real data. So the model's non-canonical precision must be judged against ~0.505, NOT 1.0.
**Result -- MODEL vs CEILING (2026-07-24, Wang dump 35889941 + dropin 35889942, `results/o2_liver/o2_summary.txt`):**

| class | ceiling F1(Wang,Janich) | model de-novo F1 | model/ceiling | with-Ribo F1 |
|---|---|---|---|---|
| CDS | 0.735 | 0.621 | **0.84** | 0.759 |
| non-canonical | 0.505 | 0.067 | 0.13 | 0.321 |
| ALL | 0.675 | 0.204 | 0.30 | 0.632 |

**Headline: for CDS the model reaches 84% of the between-experiment ceiling DE NOVO** (sequence + RNA-seq, NO
Ribo) -- nearly as informative as a second real Ribo-seq experiment. **Non-canonical de novo is far short
(13%).** BUT this is the UN-CALIBRATED standard arm: the model made 122,636 calls vs 17,881 real (~7x
over-called at theta=1), which tanks non-canonical precision. The two-arm calibrated result (Poisson/CDS-anchored
dial + the O3 FP-filter, and the track-2 peaky model) is the pending follow-up that directly targets this gap --
run the Poisson arm + FP-filter on the Wang de-novo calls next. Net read: CDS is solid de novo (0.84 of ceiling);
non-canonical is the open gap that calibration + the peaky model must close. Janich model-vs-ceiling (its dump
35890103) still queued for the symmetric check.
**Result -- POISSON ARM (2026-07-24, `dropin_sweep_wang_poisson.sbatch` + `o2_poisson_eval.py`,
`results/o2_liver/o2_poisson_summary.txt`): calibration CLOSES the CDS gap and PASSES the ceiling.** The
standard-arm CDS shortfall was pure over-calling (21,315 CDS calls vs ~14,700 real, precision 0.51). Poisson
at the CDS-anchored operating point fixes it -- CDS F1 vs the Janich ceiling (0.735):

| arm | CDS calls | CDS prec | CDS rec | CDS F1 | non-canon F1 (ceiling 0.505) |
|---|---|---|---|---|---|
| standard theta=1 | 21,315 | 0.514 | 0.786 | 0.621 | 0.067 |
| Poisson theta=0.05 | 16,121 | 0.665 | 0.769 | 0.713 | **0.195** (best nc) |
| **Poisson theta=0.02** | 13,251 | **0.763** | 0.726 | **0.744** | 0.174 |

At theta=0.02 the de-novo model's CDS F1 (0.744) **edges past the between-experiment ceiling (0.735)** -- calibrated
sequence+RNA-seq predicts liver CDS calls as well as a second real Ribo-seq experiment (matching the with-Ribo arm
0.759, now with NO Ribo). Non-canonical improves ~3x with Poisson (0.067 -> 0.195) but stays far below its 0.505
ceiling -- that gap is the FP-filter + peaky-model's job, not the depth dial. **Two-arm O2 verdict: CDS de novo
reaches the reproducibility ceiling once calibrated; non-canonical is the remaining open problem.**

## Check 8 -- robust novel-ORF caller (O3)    **Status: track 1 DONE; track 2 DONE (NEGATIVE)**
Two tracks, both scored on the SAME CDS-anchored dial + two-arm + over-call framework:
- **Adopt an ML-based caller.** Candidates: RibORF (SVM on 3-nt periodicity + %-max-entropy + read
  uniformity), RP-BP (Bayesian periodicity + offset), DeepRibo (CNN on signal+sequence), Ribo-TISH,
  ORFquant, PRICE. First filter = which can run on a SUBSTITUTED per-nt density (like
  ribocode_dropin.py) so the model's predicted profile drives them WITHOUT a BAM. RibORF's SVM
  feature set is the closest fit to a drop-in FP-filter.
- **Lightweight FP-filter classifier.** Per-candidate-ORF features from the predicted profile:
  frame-0 fraction, periodicity strength, coverage entropy/uniformity, ORF length, start-codon
  context (Kozak-free), RNA-seq coverage, 3'-of-stop drop-off. Model: logistic / gradient-boosting,
  applied as a post-filter after candidate enumeration.
- **Label source (the enabler):** O2 cross-dataset reproducibility (Wang<->Janich) + the
  over-call-vs-deep-reference check already produce real/false ORF labels -> the FP-filter is a
  supervised problem needing NO new experiments.
**Deliverable:** a caller/filter beating RiboCode+Poisson novel precision at matched CDS recall,
ideally BAM-free. Current Poisson-calibrated RiboCode is the interim baseline it must beat.
**Mini-plan (2026-07-24):** `docs/o3_orf_caller_miniplan.md` -- survey (verdict: no wholesale co-opt; every
periodicity caller inherits RiboCode's bug on a smooth density) + two tracks + phased execution + full results.
**Result (2026-07-24):**
- **Phase 1 refuted the uniformity premise:** the model over-smooths REAL ORFs too, so PME/uniformity alone
  gives AUC ~0.50 for novel; the discriminative shape info is largely ABSENT from the model's output. Best
  multivariate CV = logreg = RF = 0.669 (novel) -- linear + modest.
- **Track 1 (FP-filter, `fp_filter.py`) DONE:** per-class balanced logreg on the shape features. Novel base
  0.235 -> AUC 0.668; at fixed CDS recall 0.956 it lifts novel precision 0.294 -> 0.398 (+0.10) and uORF
  0.425 -> 0.540 (+0.12), WITHOUT the CDS-recall cost that lowering theta incurs (theta=0.05 hits 0.590 novel
  but drops CDS recall to 0.925). Complementary levers. In-sample (Hepatocytes-trained/applied) -> optimistic
  upper bound; cross-dataset (Wang/Janich) pending GPU. `fp_filter.py apply` IS the `--fp-filter` step.
- **Track 2 (anti-smoothing model) TRAINING:** `model.py::profile_entropy_gap` penalizes `relu(H(pred)-H(obs))`
  so the model stops hedging with smooth profiles; `train_loto_noBrain_peaky.sbatch` (PK=0.2) queued as
  35890525, an exact A/B against the current no-Kozak mm1 model. If it makes shape features discriminative +
  reduces over-calling, the smoothness was the bottleneck and track 1's ceiling lifts.

---

## Notes
- **Denoising is already solved** (`pred_obsdepth` clean when Ribo exists); calibration serves the
  **de novo / no-Ribo** case.
- **Why CDS as anchor:** the one ORF class we trust unconditionally, spanning the full signal range;
  its recall is a depth- and dataset-agnostic yardstick.
- **O2 annotation caveat:** Wang and Janich observed calls must be on the same mouse GENCODE annotation
  for ORF_ID matching; verify in Check 0 before comparing.

## Supplemental figure data (PRESERVE)
The CDS-anchored sweep is a planned supplemental figure. Source files (persistent, group fs), under
`results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes/`:
- `dropin_sweep/cds_anchored_curve.tsv`          lever (a) depth-scale
- `dropin_sweep_poisson/cds_anchored_curve.tsv`  lever (b) Poisson injection  <- RECOMMENDED dial
- `dropin_sweep_pval/cds_anchored_curve.tsv`     lever (d) pval threshold
- `dropin_sweep/lever_comparison.txt`            Check-4 head-to-head
Each row: theta, cds_recall, cds_precision, calls_over_cds, {uORF,novel,dORF}_{n,prec}. The per-theta
`theta_*/pred_preddepth_collapsed.txt` (raw RiboCode calls) are kept too. Do NOT delete.

**Recommended (Poisson) dial, HUMAN Hepatocytes (Chothani LOTO holdout):** at ~90% CDS recall
(theta=0.05): ~226 novel @0.58, 660 uORF @0.73; stricter (theta=0.02, CDS recall 0.835): ~86 novel
@0.59. Poisson preserves strong CDS at low theta (CDS recall floors ~0.84) while noise-rejecting weak
novels -- the desired behaviour.
**Seed-stability (3 seeds 0/1/2, `dropin_sweep_poisson/seed_stability.txt`): DIAL IS STABLE** -- per
theta the min..max ranges are tight (CDS recall +/-0.004, novel precision +/-0.01, novel count +/-~15 of
~226). Figure-ready with error bars; the single-seed dial is representative. Naming: "human Hepatocytes"
= Chothani LOTO holdout (this dial); "mouse liver (Wang/Janich)" = the separate cross-species O2 sets.

## Progress log
- 2026-07-24: plan created; O2 (between-dataset ceiling) added as headline validation.
- 2026-07-24: Checks 2/3/4 DONE on Hepatocytes. Sweep bug caught (--subsample was a no-op on
  pred_preddepth) -> added --pred_scale/--pred_poisson. CDS-anchored dial built; **Poisson injection
  is the winning lever** (+0.10 novel precision at matched CDS recall); over-calling is a smoothness,
  not count-scale, problem. Two-arm protocol (standard + downsample+Poisson) made standing.
- 2026-07-24: Janich RNA fastq found corrupt (multi-connection download); re-downloading clean via
  prefetch (task #35). Core O2 unaffected (uses Janich observed calls). Wang + GSE dumps GPU-queued.
- 2026-07-24: added **O3 / Check 8** (robust novel-ORF caller beyond RiboCode's frame test; ML FP-filter
  trained on O2 reproducibility labels) as a forward goal.
- 2026-07-24: **O3 survey done** -> verdict: no wholesale co-opt (every periodicity caller inherits the bug
  on a smooth density); keep RiboCode for candidate generation + build a shape-based FP-filter
  (`docs/o3_orf_caller_miniplan.md`).
- 2026-07-24: **O3 Phase 1** (`orf_features.py`) -> uniformity premise REFUTED (PME AUC ~0.50 novel; model
  smooths real ORFs too); best multivariate CV = logreg = RF = 0.669. Strategic pivot: the fix is upstream
  (make the model predict structured profiles), not a richer post-hoc filter.
- 2026-07-24: **O3 track 1 (FP-filter) DONE in-sample** (`fp_filter.py`). At fixed CDS recall 0.956: novel
  precision 0.294->0.398 (+0.10), uORF 0.425->0.540 (+0.12); complementary to theta (which costs CDS recall).
  `apply` is the `--fp-filter` step. Cross-dataset (Wang/Janich) validation pending GPU.
- 2026-07-24: **O3 track 2 (anti-smoothing model) built + queued** -- `model.py::profile_entropy_gap` +
  `train_loto.py --peakiness_weight`; `train_loto_noBrain_peaky.sbatch` PK=0.2 launched as 35890525 (exact A/B
  vs the current nokozak mm1 model), behind the O2 dumps in the GPU queue.
- 2026-07-24: added the **CURRENT STATE snapshot** at the top of this doc (queue + landed + pending actions).

---

## Check 9: CDS-anchoring the macrophage MODEL-DB threshold (2026-07-31, LANDED)

**Context.** Checks 2-6 anchor the Poisson/theta depth dial routed through RiboCode's frame test. That
stack was never wired into the macrophage proteogenomics path: a grep of `proteogenomics/` for
`theta|pred_scale|pred_poisson|cds_recall` returns zero hits. The macrophage MODEL database is instead
selected by `build_a549_dbs.py --thresh 0.5`, i.e. every candidate ORF whose model-predicted in-frame
fraction `pred_frame0` clears a FIXED 0.5. That 0.5 was never calibrated against anything.

Two problems with a fixed cut:
1. It is not commensurate across checkpoints. Two models can place the same ORF on opposite sides of 0.5
   purely because their score distributions differ, so an attn-vs-mamba4 DB comparison at a shared fixed
   cut confounds "which model ranks ORFs better" with "which model emits larger numbers".
2. On the deployed checkpoint it is far too permissive (below).

**Method.** `proteogenomics/scripts/calibrate_f0_threshold.py` applies the same CDS-anchoring principle
to the caller the macrophage path actually uses: sweep `pred_frame0`, measure recall of ANNOTATED-CANONICAL
ORFs (`orf_class == "canonical"` in `candidates.tsv`), and take the most stringent cut still retaining the
target CDS recall (0.90, matching the theta dial). Like `calibrate_dial.py` this needs **no observed
Ribo-seq** -- the anchor is the GENCODE CDS annotation -- which is essential here because the macrophage
populations have no ribosome profiling at all. It deliberately does NOT route through RiboCode: Task 19b
showed that caller collapses to CTG on predicted profiles, which is why `enumerate_score_orfs.py` scores
ORFs directly from the predicted profile in the first place.

**Result -- the legacy 0.5 cut barely filters.** BMDM, deployed checkpoint:

| arm | thresh | CDS recall | novel ORFs | dORF share of novel |
|---|--:|--:|--:|--:|
| legacy fixed (uncalibrated) | 0.50 | **0.997** | 55,479 | 62% |
| CDS-anchored @ 0.90 recall | **0.74** | 0.936 | 7,037 | 28% |

At 0.5 the model retains 99.7% of annotated CDS, i.e. it is essentially not filtering on the canonical
axis, and admits **7.9x more novel ORFs** than an anchored cut. Most of the excess is 3'UTR dORFs.

Per-checkpoint anchored thresholds (median over 12 populations, target CDS recall 0.90):

| checkpoint | thresh* | spread | novel @ thresh* (BMDM) | novel @ 0.5 (BMDM) | ratio |
|---|--:|---|--:|--:|--:|
| deployed attn (pre-union) | 0.74 | 0.74 - 0.74 | 7,037 | 55,479 | 0.13 |
| `orf_v2_mamba4` union | 0.74 | 0.72 - 0.74 | 18,784 | 73,614 | 0.26 |

Two things worth noting. The anchored threshold is **stable across all 12 populations** (0.72-0.74),
so it is tracking the model's score distribution rather than population noise. And both checkpoints
anchor to nearly the SAME threshold while retaining very different novel counts there (mamba4 keeps
~2.7x more) -- which is the point of anchoring: it equalises the canonical operating point so the
remaining difference is genuinely about how each model ranks non-canonical ORFs relative to canonical.

**Consequence for the published churn table.** The model DB in `macro_model_vs_null.md` was built at 0.5,
so it was much closer to the null DB than intended (45,117 vs 187,043 unique novel sequences, only 4.1x
apart). That compresses the model-vs-null contrast, so the reported 440 vs 374 novel peptides is likely
an UNDERSTATEMENT of the model's advantage.

**Two-arm protocol for the rebuild.** Changing checkpoint and threshold in one step would confound them,
so every rebuilt checkpoint gets both: `db<TAG>` at the legacy 0.5 (isolates the checkpoint effect) and
`db<TAG>_cal` at the CDS-anchored cut (isolates the threshold effect). Driver:
`proteogenomics/scripts/macro_rebuild_chain.sh <TAG> <dump_jobid>`.

**Infrastructure note.** `orf_v2_mamba4` is GPU-only: `mamba_ssm` ships CUDA kernels with no CPU
fallback and is present only in the `rnazoo-orthrus` image (1.2.0.post1, torch 2.2.0+cu118). A
`--device cpu` dump dies at import. The attn checkpoint runs on either, and on GPU a full population
dump takes ~5 min against an 8 h CPU budget.
