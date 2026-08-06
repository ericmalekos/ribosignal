# Design reference: magnitude, reference depth, and multi-dataset training

Scope: how the model predicts translation MAGNITUDE (the count head), why the profile
(shape) already transfers, why the current count head does not transfer cleanly, and how
including more datasets changes things. This is a design-rationale reference, kept separate
from the working `methods.md`, the `results.md` log, and `manuscript/methods.md`. Its central
decision was RESOLVED 2026-07-15 (see Section 9): the count head transfers, so no refactor is
needed for deployment; the refactor only matters for a future multi-dataset training expansion.

## 0. TL;DR

- The model has two heads. The **profile** head (shape: where ribosomes sit) is the
  deliverable and it transfers across datasets and species, because shape is set by sequence
  and the genetic code, which are universal, and because the ORF caller's decision is
  shape-driven.
- The **count** head (magnitude: how many footprints total) is the only component with a
  transferability question. As built it predicts an ABSOLUTE count on the training data's
  scale, i.e. "expected footprints if you sequenced at the training reference depth." That is
  a fine units choice for the use case (a dataset with RNA-seq but no Ribo-seq has no native
  Ribo-seq depth to match anyway).
- BUT the count target is currently ABSOLUTE, and that makes multi-dataset training
  INCOHERENT: the same normalized input maps to different absolute targets across datasets of
  different depth, so the head averages and eats noise. So, as currently defined, adding
  datasets makes the count head worse, not better.
- The fix that makes "more datasets" help the count head is to put its target on a
  dataset-invariant scale (translation efficiency = footprints per unit RNA, or a per-million
  share of the translatome). Then multi-dataset training is coherent, each dataset adds real
  signal, and reference depth becomes an explicit inference-time knob.
- Whether the count head is worth refactoring was settled by the held-out drop-in (Section 9):
  the `pred_obsdepth` vs `pred_preddepth` F1 gap is ~0 (human, all three backends), so the count
  head transfers and the refactor is NOT needed for deployment -- only for multi-dataset training.

## 1. Two heads, two fates

`RiboSignalModel` (see `manuscript/methods.md` Section 7) emits:

- **Profile head** (`profile_multinomial_nll`): a softmax distribution over transcript
  positions. The loss divides out the per-transcript total (`n = counts.sum()`), so it is
  trained ONLY on the shape.
- **Count head**: a single scalar regressed to `log1p(total P-sites)` -- the absolute pooled
  footprint count on the transcript.

The two behave differently under dataset change because they depend on different things.

## 2. Why the profile (shape) transfers

Where ribosomes sit is set by start/stop codons, reading frame, the ribosome's 3-nucleotide
step, ORF structure, and Kozak context. Those are properties of the genetic code and the
cytoplasmic translation machinery, which are identical across human tissues, across a human
cell line the model has never seen, and across the human/mouse boundary. A model that learned
"given this sequence and its relative RNA abundance, here is the profile" is learning
something a new dataset also obeys.

Crucially for downstream ORF calling: RiboCode's `detectORF` decision is shape-driven. The
p-value comes from a Wilcoxon test of frame0 vs frame1/frame2 per-codon counts (periodicity);
the absolute count only gates "is there enough signal to test" (per-transcript sum >= 5,
per-ORF frame0 sum >= 5, >= 5 nonzero frame0 codons). So ORF identity and reading frame come
from the profile; the count only has to be "deep enough to clear the gate."

## 3. The count head as built = a reference-depth prediction

The count head's INPUT is the mean-normalized coverage channel
(`log1p(coverage / dataset_mean_depth)`), so a transcript with a given RELATIVE RNA abundance
produces the same predicted total regardless of how deep the new experiment was. Its OUTPUT is
trained against the absolute total on the training data. Net effect: it predicts "expected
footprints at the training reference depth," independent of the target dataset's own depth.

For the intended use (predict Ribo-seq for a dataset that has RNA-seq but no Ribo-seq), this is
exactly right: there is no target Ribo-seq depth to match, so a fixed reference is the only
well-defined option, and a deep, high-quality reference gives a consistent, sensitive calling
threshold across datasets. Because calling is shape-driven (Section 2), the count does not need
to be ACCURATE, only in a sane range; the reference-depth prediction supplies that.

Consequence for interpreting a shallow held-out: predicting at a DEEP reference while the
target's real Ribo-seq is shallow means the standalone prediction can legitimately OUT-call the
real data (recover ORFs the shallow experiment lacked power to detect). A predicted call that
misses a real call should be checked against depth before being counted a false positive.

## 4. The problem: an absolute count target makes multi-dataset training incoherent

The deployed LOTO models already train on 8 Chothani tissues (not one), so the count head's
reference is already an 8-tissue blend. That exposes the real issue.

Take a housekeeping transcript T translated in both Fibroblast (32 samples pooled, deep;
mean per-nt depth ~7,733... note the target is footprints not coverage, but the depth scales
similarly) and Brain (few samples, shallow; mean depth ~106). T's RELATIVE coverage shape is
similar in both, so the normalized coverage input is similar, and the sequence/ORF inputs are
IDENTICAL (same transcript). But the absolute total P-sites differ by the depth x
sample-count ratio, easily 10 to 50x. The count head has no depth or dataset covariate, so it
sees near-identical inputs with very different targets and can only regress to the average,
taking the rest as irreducible error.

This is why the count head (Pearson ~0.78) is softer than the profile, and why its "reference
depth" is an ill-defined blend rather than a clean unit. The profile head does NOT have this
problem precisely because its target is normalized (`counts/N`): T's profile is the same shape
in both tissues, so the target is consistent. The normalization that makes the profile
transferable is exactly what the count target lacks.

Therefore, as currently defined, adding datasets makes the count head WORSE (more scale
ambiguity, more averaging). The count head is the one component allergic to more data.

## 5. The fix: a dataset-invariant count target

Put the count target on a scale every dataset shares:

- **Translation efficiency (TE)**: footprints per unit RNA (a rate). Most biologically
  meaningful and most transferable ("how efficiently is this mRNA translated"). Both footprints
  and RNA in normalized units (per-million / TPM) makes TE dataset-invariant.
- **Per-million share**: footprints on transcript / total footprints x 1e6 ("share of the
  translatome"). Also dataset-invariant.

With a normalized target, every dataset speaks the same units. Multi-dataset training becomes
coherent (no depth-ambiguity noise), each new dataset adds real signal (sequence and condition
diversity), and the quantity is independently useful ("predicted TE"). Reference depth becomes
an EXPLICIT inference knob: predict TE, then scale by (target RNA x a chosen reference total)
to produce callable integer counts. You set the reference; it is not a training accident.

The profile head is unchanged; only the count target + loss change. This is a target change,
not an architecture change.

## 6. Including other datasets: what changes

- **Profile (shape):** more datasets straightforwardly help. More diverse translation examples
  improve generalization on the harder strata (lncRNA, uORF) and, if a second species is added
  to TRAINING, cross-species. No units issue. LOTO already shows near-lossless protein-coding
  transfer across tissues; more datasets firm up the rest.
- **Count (magnitude):** more datasets help ONLY after the target is normalized (Section 5).
  Before that, they hurt (Section 4).
- **Evaluation tension:** a dataset used for training can no longer measure transfer. Keep
  Ruiz-Orera and Wang (and ideally one more independent species) as STRICT held-outs; add OTHER
  public Ribo-seq to training. Otherwise the generalization claim evaporates.
- **Protocol / batch effects:** different labs differ in periodicity quality, read-length
  distributions, rRNA contamination, and 3' bias. Pool naively and the model learns protocol
  artifacts. Mitigations: per-dataset P-site calling (already done), QC-gating shallow/poor
  datasets, a common count normalization, and possibly a dataset covariate the model can use at
  train time and that is zeroed at inference.

## 7. Coverage-shape caveat (separate from depth)

Depth normalization (Section 3) handles the SCALE of the coverage input. It does not handle the
SHAPE of coverage (3' bias, coverage non-uniformity) differing across protocols. That is a
distinct, secondary transfer concern for the count head's input and is bounded empirically by
running the held-out with and without the coverage channel (`input_mode` both vs emb).

## 8. Drop-in variant semantics (how to report)

- `pred_preddepth` = predicted shape x count-head total = the fully standalone, no-Ribo-seq,
  reference-depth prediction. This is the DEPLOYMENT-REALISTIC variant and should be the
  headline for any transfer claim.
- `pred_obsdepth` = predicted shape x the REAL per-transcript total = a diagnostic that borrows
  the real Ribo-seq total to isolate whether the error is in the shape or the magnitude. It is
  NOT available at deployment (a no-Ribo-seq dataset has no real total).

Earlier framing led with `pred_obsdepth`; the correct emphasis is `pred_preddepth` as the real
result, `pred_obsdepth` as the ablation.

## 9. Decision -- RESOLVED 2026-07-15 (count head transfers; no refactor needed for deployment)

**Measurement landed.** On the human Ruiz-Orera cross-study drop-in, `pred_obsdepth` F1 and
`pred_preddepth` F1 are essentially identical across all three backends: one-hot 0.931 vs 0.929
(+0.002), Orthrus 0.925 vs 0.926 (-0.001), RiNALMo 0.924 vs 0.926 (-0.002). The gap is ~0, i.e.
`pred_preddepth ~= pred_obsdepth` (results.md Task 16).

**Decision:** the current absolute-count head TRANSFERS fine for the ORF-calling use case; the
TE / per-million refactor (Section 5) is NOT needed for single-dataset, no-Ribo-seq deployment.
The predicted depth (regressed to the Chothani reference) gates RiboCode calling as well as the
real observed depth would, so the fully standalone `pred_preddepth` path is validated on
independent human data. The earlier worry that the count head would not transfer was wrong for
the calling use case (calling is shape-driven, and the count only has to clear the signal gate).

The refactor's remaining value is narrower than first thought: it matters ONLY for coherent
MULTI-DATASET TRAINING (Section 4), where one normalized input mapping to different absolute
targets across datasets of different depth/sample-count is genuinely incoherent. It is no longer
on the critical path for deployment; prototype it if/when a multi-dataset training expansion is
taken, otherwise defer.

For the record, the decision rule this measurement settled:
- `pred_preddepth` ~= `pred_obsdepth` (OBSERVED): count head transfers fine; no refactor needed.
- `pred_preddepth` << `pred_obsdepth` (not observed): count head would be the weak link and the
  refactor would be justified.

The mouse cross-SPECIES drop-in (Wang liver, vM38) now confirms the same conclusion a second time:
`pred_preddepth` F1 0.919 vs `pred_obsdepth` 0.929 (gap -0.010), so the standalone predicted-depth
path transfers across the human->mouse boundary as well as it does cross-study in human. (This landed
2026-07-15 after fixing a hardcoded human tx->gene map in `compare_dropin_calls.py` -- see results.md
Task 16.)

## 10. Implementation sketch (only if the refactor is taken)

- Change the count target from `log1p(total P-sites)` to a normalized quantity: e.g.
  `log(footprints_CPM / RNA_TPM)` for TE, or `log(footprints_CPM)` for a per-million share.
  Requires per-dataset totals (footprint library size, RNA library size) to normalize.
- Keep the profile head and its multinomial NLL unchanged.
- At inference, choose an explicit reference total R; predicted callable counts =
  `profile x (predicted_TE x target_RNA x R)` or `profile x (predicted_share x R)`, then round
  to integers for RiboCode (the existing drop-in rounding logic).
- State R explicitly in the manuscript ("predicted at a reference depth of R footprints per
  transcriptome") instead of the implicit "8-tissue-blend x posture-A" scale.
- Re-run the held-out drop-in and confirm `pred_preddepth` closes the gap to `pred_obsdepth`.
