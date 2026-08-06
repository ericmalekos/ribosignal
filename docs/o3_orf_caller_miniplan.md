# O3 mini-plan: alternates to RiboCode for novel-ORF prediction

Companion to `docs/count_head_calibration_checks.md` (Check 8). Forward goal; scoped as a mini-effort,
CPU-feasible, reusing the existing evaluation framework + O2 labels. Not a from-scratch deep caller.

## Motivation (diagnostic, not speculative)
RiboCode calls ORFs with a SINGLE univariate periodicity test (Wilcoxon/Mann-Whitney of frame-0 vs
frame-1/2 P-sites per codon). On our SMOOTH model-predicted density this over-calls short non-canonical
ORFs: even spurious weak ORFs look cleanly periodic and pass. Quantified (Check 5, GSE-NT deployed): novel
calls validate at only 0.11 vs a 0.54 real-novel baseline against a deep reference -> genuine over-calling.
Poisson injection patches it (+0.10 novel precision) but is a workaround. Goal: a caller / FP-filter that
raises novel-ORF precision at matched CDS recall.

## The deciding requirement
The failure mode is SMOOTH spurious ORFs. The signal RiboCode ignores is **read UNIFORMITY / distribution
shape**: a spurious model ORF has artificially uniform predicted coverage, whereas a real ORF has a
characteristic ramped/peaked profile (start enrichment, pauses, 3'-of-stop drop-off). So the caller MUST
use uniformity / read-distribution features, not periodicity alone -- and should ideally run on a
SUBSTITUTED per-nt density {tx_id: np.array} (BAM-free), the way we already drive RiboCode's detectORF.

## The approach (post-survey): RiboCode candidate generation + a custom FP-filter
Track A (co-opt a caller wholesale, BAM-free) is DROPPED -- no surveyed tool is both right-featured AND cheap
to run on our injected density (survey below). RiboCode stays the CANDIDATE GENERATOR (+ the Poisson dial for
the operating point); a lightweight FP-filter re-scores its calls to kill the smooth-spurious novel ORFs.

**FP-filter feature vector** (all computed in numpy on the predicted `{tx_id: np.array}` -- no BAM, no genePred):
- **RibORF-derived core:** f1/f2/f3 in-frame read fractions; **PME** (percentage-of-max-entropy = codon-level
  uniformity); **f1max** (fraction of codons whose in-frame position dominates). The multivariate periodicity
  +uniformity signal RiboCode lacks.
- **Over-smoothing detectors:** coverage **Gini / coefficient-of-variation / max-to-mean / fraction-zero**
  (a real short ORF stays spiky+ramped in the model output; a spurious one flattens); **5' start-ramp**; **3'-of-
  stop drop** (ribosome release; ORFscore / RRS-style).
- **Context:** ORF length, start-codon identity/context (Kozak-free), RNA-seq coverage support.

**Classifier:** logistic regression / gradient-boosting on real-vs-spurious labels from O2 reproducibility (an
ORF is real if reproduced in a deep INDEPENDENT experiment). DeepRibo's framing (a learned classifier over
coverage shape + signal-to-noise) is the template; discard its prokaryote Shine-Dalgarno features.

## The training signal (the enabler)
O2 cross-dataset reproducibility (Wang<->Janich) + the over-call-vs-deep-reference check already produce
real/false ORF labels with NO new experiments: an ORF is "real" if reproduced in a deep INDEPENDENT
experiment. So (B) is a supervised problem, and (A)'s co-opted model can be re-scored/re-thresholded on
these labels.

## Evaluation (reuse the existing framework, no new metrics)
Score EVERY candidate on the SAME CDS-anchored dial + two-arm + over-call framework:
- non-canonical precision at matched CDS recall, vs the RiboCode+Poisson baseline (Hepatocytes novel 0.58);
- the BAR = the O2 non-canonical reproducibility ceiling (F1 0.505) -- reaching it = as good as a 2nd real
  experiment;
- validated against a DEEP independent reference (never a shallow same-dataset one).

## Candidate survey (2026-07-24 research pass) -- VERDICT: no wholesale co-opt; re-implement RibORF's features
Deciding fact: EVERY periodicity-based caller inherits the bug. A smooth predicted density is trivially
periodic wherever the model puts in-frame P-sites, so RiboCode's Wilcoxon frame test, Ribo-TISH's rank-sum,
Ribotricer's phase score, RP-BP's Bayesian periodicity, and ORFquant's multitaper spectrum ALL pass spurious
short ORFs the same way. The discriminating signal is read-distribution SHAPE / UNIFORMITY, not periodicity.

| tool | features | uniformity? | BAM-free | license | fit |
|---|---|---|---|---|---|
| **RibORF** | f1/f2/f3 in-frame + **PME** (codon uniformity) + **f1max** | **YES** | features re-implementable in numpy (~40 lines); pretrained SVM won't transfer to predicted density | no LICENSE file (docs say GPL-3.0) | **the right feature set** |
| Ribotricer | phase score only | no | easiest host (Python, array profile) | GPL-3.0 | periodicity-only -> same bug |
| RP-BP | Bayesian periodicity | no | hard (self-selects offsets from BAM) | **MIT** | same bug |
| Ribo-TISH | Wilcoxon rank-sum (= RiboCode family) | no | BAM hard-wired | GPL-3.0 | same bug |
| ORFquant | multitaper spectral (most precise in benchmarks) | no | very hard (R/RiboseQC S4) | GPL-3.0 | precise via periodicity+coverage, NOT uniformity |
| PRICE | EM read-distribution model | partial | very hard (Java/Gedi) | GPL-3.0 | conceptually apt, impractical to drive from a density |
| DeepRibo | CNN(seq)+RNN(coverage) ML | learned | prokaryote-only (SD features) | GPL-3.0 | wrong features; borrow the FRAMING |

Newer ML (RiboTIE/TRISTAN transformer, smORFer; 2022+) don't fit: transformer callers just learn the
smoothness; smORFer is prokaryote/Fourier. No eukaryotic ML FP-filter for "smooth predicted density" exists;
DeepRibo (2019) is the closest template. Benchmark (BiB 2024 bbae268): RiboCode + RibORF are the sORF
over-callers, ORFquant most precise; and the benchmark's own recommended FP filter "fraction in-frame" will
NOT help us (our model makes in-frame high everywhere) -- confirming the fix is a SHAPE feature.

**VERDICT:** drop track A (no tool is both right-featured AND cheap to run BAM-free on an injected density).
Keep RiboCode for CANDIDATE GENERATION. The path is track B: a custom FP-filter that RE-IMPLEMENTS RibORF's
features (license-clean -- PME/entropy are standard statistics, not copyrightable code) plus shape features,
trained on O2 reproducibility labels.

## Phased execution (CPU-feasible)
- **Phase 0 -- survey**: DONE (2026-07-24). Verdict: no wholesale co-opt; re-implement RibORF's features +
  shape features in a custom FP-filter. Track A dropped.
- **Phase 1 -- features + labels**: DONE (2026-07-24). `orf_features.py` extracts the feature vector per
  candidate ORF from the predicted `{tx_id: np.array}` (drop-in candidate enumeration, no BAM). Labels: real
  if the ORF's coord key is in the deep observed calls. Ran on Hepatocytes theta=2.5 pool -> 24,318 ORFs
  (`o3_features_hep.tsv`). **KEY FINDING (uniformity premise REFUTED):** PME (codon uniformity) alone gives
  AUC ~0.50 for novel -- the model smooths REAL ORFs too, so the shape signal RiboCode lacks is largely absent
  from the model's own output. Best single feature frac_f0 (in-frame fraction) AUC ~0.62; multivariate 5-fold
  CV logreg = RF = AUC 0.669 (novel), i.e. the signal is LINEAR and modest, not a rich separable manifold.
- **Phase 2 -- baseline**: RiboCode+Poisson at the CDS-anchored operating point (in hand).
- **Phase 3 -- train the FP-filter**: DONE (2026-07-24). `fp_filter.py train` -- per-ORF-class balanced
  logistic regression (+ a pooled non-canonical model), 5-fold CV, saves scaler+coef+per-recall thresholds.
  Honest within-dataset CV (Hepatocytes, `fp_filter_hep.joblib`):

  | class | base rate | AUC | P@recall .9 | P@recall .5 | lift@.5 |
  |---|---|---|---|---|---|
  | novel | 0.235 | 0.668 | 0.268 | 0.356 | 1.52x |
  | uORF | 0.361 | 0.647 | 0.391 | 0.500 | 1.39x |
  | dORF | 0.035 | 0.644 | 0.038 | 0.057 | 1.63x (still ~0.05 abs -- dORF is 96.5% spurious) |
  | annotated (CDS) | 0.936 | 0.591 | 0.949 | 0.950 | 1.01x (untouched, as intended) |
  | POOLED noncanon | 0.261 | 0.674 | 0.299 | 0.410 | 1.57x |

- **Phase 4 -- compare**: DONE (2026-07-24). `fp_filter.py apply` on the CDS-anchored / two-arm framework.
  **The FP-filter's niche is precision at FIXED CDS recall** (theta is a global dial that trades CDS recall for
  precision; the filter targets ONLY non-canonical calls, leaving CDS alone). Hepatocytes, in-sample:

  | operating point | CDS recall | novel prec | uORF prec | dORF prec |
  |---|---|---|---|---|
  | theta=1.0 raw (standard arm) | 0.956 | 0.294 | 0.425 | 0.062 |
  | theta=1.0 + FP-filter(recall .5) | **0.956** (held) | 0.398 (+0.10) | 0.540 (+0.12) | 0.098 |
  | theta=0.05 raw (low-theta) | 0.925 (-.03) | 0.590 | 0.733 | 0.274 |

  Read: the FP-filter buys +0.10 novel / +0.12 uORF precision WITHOUT costing CDS recall; lowering theta buys
  MORE precision but at a CDS-recall cost. Complementary levers -- combine for the highest precision.
  **Redundancy caveat:** at the already-stringent theta=0.05 point the filter at recall .9 barely moves
  anything (68/12,010 dropped) -- stringent theta and the shape-filter remove the SAME smooth-spurious ORFs,
  so the filter's marginal value is largest when you must hold CDS recall high (can't lower theta).
  **In-sample caveat:** the filter is Hepatocytes-trained and applied to Hepatocytes calls (a superset it was
  fit on), so +0.10 is an optimistic upper bound. The honest CROSS-dataset number (train Hepatocytes, apply
  Wang/Janich) lands when the GPU dumps clear.
- **Phase 5 -- integrate**: DONE (2026-07-24). `fp_filter.py apply --model <joblib> --profiles
  pred_profiles.npz --candidates <collapsed.txt> --out <filtered.txt> [--target-recall R | --threshold T]`
  IS the post-RiboCode `--fp-filter` step: drop-in emits `pred_profiles.npz` + `*_collapsed.txt`; the filter
  re-scores and rewrites a filtered `*_collapsed.txt`. Deployment flow = train once on a labeled dataset,
  apply to new data. Left as a standalone step (the filter model is dataset-specific; not baked into the
  drop-in).

**VERDICT (track 1):** the FP-filter is a REAL but MODEST lever (~0.67 AUC, +0.10 non-canonical precision at
fixed CDS recall). Its ceiling is set by the model's smoothness (Phase 1: real ORFs are smoothed too, so the
discriminative shape info is largely absent from the output). The higher-ceiling fix is track 2 (anti-smoothing
model) -- if that model predicts structured profiles, re-running `orf_features.py` + `fp_filter.py train` on ITS
output should show higher feature-separation AUC and a bigger precision lift. That is the direct A/B test.

Phases 1-5 done CPU-only from the Hepatocytes predictions in hand; Wang/Janich add the honest cross-dataset
labels once their GPU dumps land.

## Track 2 -- anti-smoothing model (the higher-ceiling fix; trained IN PARALLEL with the current model)
Phase 1 refuted the premise that shape features can rescue precision post-hoc: the model over-smooths REAL
ORFs, so the discriminative structure is not in its output. The fix is upstream -- make the model predict
structured/peaked profiles. Built 2026-07-24, training queued (35890525) behind the O2 dumps so they keep
priority.

**The loss term** (`model.py::profile_entropy_gap`, wired via `train_loto.py --peakiness_weight`): the
profile head's multinomial NLL is minimized by matching the empirical profile, but a smooth prediction still
scores well (cross-entropy tolerates a hedged distribution). The new term penalizes the Shannon-entropy excess
`relu(H(pred) - H(obs))` per transcript -- how much SMOOTHER the prediction is than the (deep, ~6.6 P-sites/nt,
so real-periodicity-bearing) pooled target. relu keeps only over-smoothing; a prediction peakier than the noisy
target is never pushed back. Smoke-tested: penalizes smooth >> peaky (0.96 vs 0.003), gradient nonzero and
descent reduces the gap at realistic (non-uniform) points; Adam accumulates the shallow near-uniform gradient.

**The run** (`train_loto_noBrain_peaky.sbatch`, PK=0.2): identical model / config / holdout (Hepatocytes) /
no-Kozak mm1 track as the current model (`..._nokozak_mm1_holdout_Hepatocytes`, already trained) PLUS the
entropy-gap term -> `..._nokozak_mm1_peaky0p2_holdout_Hepatocytes`. Directly comparable.

**The A/B test** when it finishes: (1) does val_pearson hold (peakiness must not wreck profile accuracy)?
(2) re-run `orf_features.py` on the peaky model's drop-in -> does feature-separation AUC rise above 0.669?
(3) re-run `fp_filter.py` + the CDS-anchored table -> does novel precision at fixed CDS recall beat 0.398?
(4) does the peaky model OVER-CALL less at matched CDS recall (fewer raw novel calls)? If yes on (2)-(4), the
smoothness was the bottleneck and the fix generalizes; if val_pearson tanks, dial PK down (0.05-0.1) and rerun.

**A/B RESULT (2026-07-25, `peaky_ab.sbatch`, `peaky_ab_35906078.out`) -- NEGATIVE. Track 2 does not pay off.**
- (1) val_pearson dropped modestly: peaked at 0.50 (epoch 6) vs the standard nokozak-mm1 0.527 -- a ~5% loss
  (CORRECTION: the earlier "22% vs 0.64" was a wrong baseline; 0.64 was the older deployed mm20/heuristic model,
  not the nokozak-mm1 standard this is a variant of). The verdict below is unchanged -- it rests on (2)/(4).
- (2) feature-AUC did NOT improve: novel 0.656 (peaky) vs 0.669 (standard); pooled non-canon 0.671 vs 0.674.
  The anti-smoothing forces peakiness on REAL and SPURIOUS ORFs alike, so a peakier density separates them no
  better than a smooth one. The FP-filter's ~0.67 ceiling is fundamental, NOT a smoothness artifact.
- (4) over-calling reduced only trivially: at theta=1, novel 1271->1150 (-9.5%), dORF 733->443 (-40%), uORF
  5091->4982 (-2%), annotated preserved (11160->11157). Real but negligible next to the Poisson dial, which
  cuts novel over-calling ~85% FOR FREE (122K->17K on Wang) without hurting val_pearson.
- **VERDICT: bad trade, and a gentler PK cannot rescue it** -- the strongest anti-smoothing setting already
  failed to move the AUC, so a weaker push can only do less. The FP-filter + Poisson calibration remain the O3
  answer. Value of the experiment: it RULES OUT "the model is just too smooth" and pins the FP-filter ceiling
  to a real limit -- the model cannot distinguish real from spurious non-canonical ORFs from sequence + RNA
  alone, however the output is shaped. The `--peakiness_weight` term stays in the codebase (default 0, off).

## Decision / kill criteria
- **Adopt** if a candidate beats RiboCode+Poisson novel precision at matched CDS recall by a meaningful margin
  (target +0.05) without costing CDS recall, AND is BAM-free-able + license-clean for redistribution.
- **Defer / keep RiboCode+Poisson** if nothing beats the baseline or adaptation/licensing cost is prohibitive.

## Scope guard
Mini-effort only. A from-scratch end-to-end learned ORF caller on the predicted density is the bigger future
option -- noted, out of scope here. Do NOT lose sight of the pending GPU dumps (Wang/Janich/GSE) that feed the
O2 model-vs-ceiling; those come first.
