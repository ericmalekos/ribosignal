# Autonomous work plan, 2026-08-11 (~2-3 days)

Agreed mode: work through this list with minimal input. After EACH step, update `methods.md`,
`results.md`, and the new tutorial section (below). Commit at each milestone; push at the end of each
phase.

New tutorial section: **"What the model learned"** (`tutorial/source/learned/`), covering the
memorisation control, the attribution results, and the codon analysis. Registered in the toctree
under its own caption.

## Guiding constraints (carried from established project rules)

- Every ORF-call number: `build_loader` filtering, genomic key, expression-restricted, state `n_ref`.
- Both prediction arms, always.
- `--kozak none` everywhere.
- No claim about a class with n<100 reference calls without an interval.
- Attribution primary substrate is MOUSE (transcripts never seen); human Hepatocytes is the
  deliberately CONTAMINATED arm, not the trusted one. See `attribution_experiments_plan.md`.

---

## Phase 1 — close the Chothani thread (in flight)

**1.1** Wait out `align_chothani_refetch` (74 runs). Verify: 74 BAMs, all `SO:coordinate`, indexed,
mean unique% recorded, zero ABORT in logs.

**1.2** Pool P-sites per tissue (8 tissues) from the new BAMs via `prepare_from_bams` on each
tissue's existing universe. Per-tissue metaplots, as standard.

**1.3** **Measure the divergence** against the existing training packs' `target_counts.npy`, per
tissue. This is the deliverable, not a match: predict a Janich-style ~1-2% gap. Report per-tissue
total P-sites old vs new, nt differing, and % of signal.

**1.4** Write `docs/chothani_regeneration.md` + results.md. If divergence is small, the training data
is now reproducible and that is worth stating plainly; if large, state what that implies for the
released models (nothing retroactively, but future retrains will not reproduce them).

## Phase 2 — the memorisation control (cheap, settles an open question)

**2.1** Score the **98 Hepatocytes transcripts never seen in any held-in tissue** against a
size- and signal-matched sample of seen transcripts, both models, both arms. Report per-transcript
profile correlation and ORF-call agreement, with intervals (n=98 is small).

**2.2** If the unseen set performs materially worse, the 99.7% overlap IS inflating held-out numbers
and every Hepatocytes figure needs a caveat. If it does not, the mouse evidence stands and the
overlap is cosmetic. Either result is publishable and must be written up as found.

**2.3** Document in the new tutorial section + `attribution_experiments_plan.md`.

## Phase 3 — attribution (tests 1 and 2)

**3.1** Test 1, aggregate saliency: mouse GSE243134 (primary) and Hepatocytes (contaminated arm),
both models. Aggregate |d target / d one-hot| by position relative to start/stop, by frame, by ORF
class. Settles whether the C10 n=3 asymmetry (start-codon rank 3 for uORF vs 292 for CDS) is real.

**3.2** Test 2, ORF-track channel ablation: zero each of the 5 channels, re-dump, re-call, record
per-class F1 delta. Decides whether `internal` (0.022 vs a 0.539-0.620 ceiling) fails for want of
SEQUENCE signal or COVERAGE signal.

**3.3** Figures + tutorial page per test. Report the contaminated-vs-clean comparison explicitly:
features present in both substrates are learned rules, features only in human are memorisation
candidates.

## Phase 4 — codon-signal relationship (the new experiment)

The 99 existing `pred_profiles.npz` each carry BOTH `obs_flat` and `pred_flat` on the same
transcripts, so the empirical and model codon statistics are the SAME computation on two arrays. No
new inference is needed for the core result.

**4.1** `scripts/codon_occupancy.py`: for every transcript with an annotated CDS and sufficient
signal, normalise per-nt counts by the transcript's CDS mean (removes expression), assign each codon
position to its codon identity, and average genome-wide. Emit a 61-vector (sense codons) per
(dataset, source) where source in {observed, predicted}. Report BOTH P-site and A-site (P+3)
assignment; the A-site is the conventional dwell proxy and the choice must not be silent.

**4.2** **Empirical reproducibility ceiling for codon dwell.** Correlate the observed 61-vectors
BETWEEN datasets (mouse liver x3, Chothani tissues, THP-1, CAR-T). Codon dwell is known to be
protocol-sensitive (cycloheximide artefacts), so this ceiling is required before any model
comparison means anything -- the same logic as the ORF-call ceiling.

**4.3** **Model vs empirical.** Correlate `pred_flat`-derived codon dwell against `obs_flat`-derived,
per dataset. The model never receives codon identity or tRNA abundance as input, so any recovered
structure was learned from sequence.

**4.4** **The decisive split: does it generalise?** Reproducing codon dwell on Chothani (training
tissues) is partly just fitting the training data. Reproducing it on MOUSE, where no transcript was
ever seen and the species differs, is evidence of a learned codon rule. Report these separately and
do not average them.

**4.5** Relate to attribution: does per-codon model dwell correlate with per-codon saliency from
Phase 3? Two independent routes to the same quantity agreeing is much stronger than either alone.

**4.6** Figure + tutorial page. If codon dwell is NOT recovered, say so -- a negative result here
bounds what the model can be claimed to have learned, and is worth as much as a positive.

## Phase 5 — conditional follow-ups

**5.1** ISM around start codons ONLY IF Phase 3.2 shows `internal`/`dORF` fail on sequence rather
than coverage. If they fail on coverage, ISM is the wrong experiment and this is skipped, with the
reason recorded.

**5.2** `#75` matched-depth arm for the mouse 3x3. The factorial predicts the depth axis is flat
above ~7 RNA samples; this tests that directly by subsampling GSE243134's 19 RNA samples to 7 and 2.

## Phase 6 — consolidate

**6.1** Assemble `tutorial/source/learned/` into a coherent section, rebuild, verify served.
**6.2** Refresh `results.md`, `methods.md`, `manuscript/FIGURES_PLAN.md` (append, never overwrite).
**6.3** Final commit + push.

## Reporting rule for this run

Report at the end of each PHASE, not each step, unless something is blocked or a result contradicts
an earlier claim. Contradictions get reported immediately -- the point of this plan is to find them,
and a plan that only reports confirmations is not worth running.
