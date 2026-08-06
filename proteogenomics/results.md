# Proteogenomics pilot -- results

Headline question: does the Ribo-seq P-site model, scoring ORFs from a query cell type's own RNA-seq (no
Ribo-seq for that sample), select a translated-ORF search space that finds novel peptides at a higher rate
than a naive all-candidate background, and competitively with an experimental Ribo-seq ORF DB?

All novel-peptide counts are at 1% CLASS-specific FDR (novel = peptide maps only to non-canonical ORFs,
scored against REV_nuORF| decoys only). See `methods.md` for the pipeline and `DATA_PROVENANCE.md` for data.

## UC-A: A549 whole-proteome (cross-project RNA-seq -> CCLE TMT proteome) -- SUPERSEDED

> Old checkpoint, f0-threshold DB builder, unfrozen `calibrate_mass = 2`. The 2.00x / 2.31x ratios
> below are replaced by "UC-A CURRENT" further down. Retained for the two-universe (in-universe vs
> fresh) comparison, which has not been re-run.

Two universes were run. "In-universe" reuses the Fibroblast training universe (covers only the A549 tx that
overlap it); "fresh universe" scores A549's OWN expressed transcriptome (the cell-type-specific setting, the
intended use). Both show the model selecting a ~2x higher novel-peptide discovery rate than the null at a
3-4x smaller DB, with less canonical-ID churn.

### In-universe (Fibroblast universe)

| DB | novel ORFs | novel peptides @1% class FDR | rate (pep/1k ORFs) | canonical PSMs |
|----|---:|---:|---:|---:|
| canonical | 0 | -- | -- | 294,903 |
| model | 53,253 | 6 | 0.11 | 279,404 |
| null | 194,881 | 11 | 0.06 | 253,131 |

Model/null discovery-rate ratio **2.00x**. Canonical churn: model -5.3%, null -14.2%.

### Fresh universe (A549's own expressed transcriptome)

| DB | novel ORFs | novel peptides @1% class FDR | rate (pep/1k ORFs) | canonical PSMs |
|----|---:|---:|---:|---:|
| canonical | 0 | -- | -- | 294,903 |
| model | 98,047 | 10 | 0.10 | 264,237 |
| null | 294,866 | 13 | 0.04 | 238,646 |

Model/null discovery-rate ratio **2.31x**. Canonical churn: model **-10.4%**, null **-19.1%** (the model's
smaller selected DB displaces far fewer canonical IDs). Novel-peptide overlap: model-only 1, shared 9,
null-only 4 (the model recovers 9 of the null's 13 and adds 1 the null's ranking buried; the null's 4 extra
come at 3x the DB and 2x the canonical-ID loss).

**Read:** on a deep tryptic whole-proteome the absolute cryptic-ORF yield is small (10-13 novel peptides),
as expected -- whole-proteome MS is not enriched for non-canonical products. The model's value shows as a
consistently ~2x higher discovery RATE and roughly half the canonical-ID churn of the naive background. The
whole-proteome cap on absolute yield is exactly why UC-B uses an immunopeptidome.

## UC-B: HBL-1 MHC-I immunopeptidome (matched RNA-seq, + PRICE comparison) -- OLD MODEL/METHOD

> Retained for the PRICE comparison, which has not been re-run. The model arm here is the OLD
> checkpoint with f0-threshold DB selection and unfrozen search params. For current model-vs-null
> numbers see "UC-B CURRENT" below.

The cryptic-ORF-enriched substrate (non-tryptic 8-14mer, label-free, ~29,800 MS2 across 2 reps). Four DBs
share one canonical base; novel-ORF counts model 70,843 < PRICE 136,238 < null 220,190.

| DB | novel ORFs | novel peptides @1% class FDR | rate (pep/1k ORFs) | canonical PSMs |
|----|---:|---:|---:|---:|
| canonical | 0 | -- | -- | 17,943 |
| model | 70,843 | 10 | 0.14 | 17,805 |
| null | 220,190 | 18 | 0.08 | 17,074 |
| price | 136,238 | **25** | **0.18** | 17,019 |

Discovery-rate ranking: **PRICE 0.18 > model 0.14 > null 0.08**. Ratios: model/null 1.73x, model/PRICE 0.77x,
null/PRICE 0.45x. Canonical churn: model -0.8%, null -4.8%, PRICE -5.2%.

Novel-peptide set overlap: model vs null = model-only 0 / shared 10 / null-only 8 (model is a clean subset
of null). model vs PRICE = model-only 4 / shared 6 / PRICE-only 19. null vs PRICE = shared 11.

**Honest read.** On the matched immunopeptidome the authors' PRICE DB WINS: it uses actual same-sample HBL-1
Ribo-seq, and beats the model on both absolute novel peptides (25 vs 10) and discovery rate (0.18 vs 0.14).
The model, a sequence + RNA-seq translation prior with NO Ribo-seq for the query sample, does two useful
things: (1) it beats the naive 3-frame null by 1.73x rate with a clean-subset selection (its 10 are all
inside the null's 18) and near-zero canonical churn (-0.8% vs the null's -4.8%); (2) it recovers 4
immunopeptides PRICE's Ribo-seq missed (complementary signal). But it does NOT replace matched experimental
Ribo-seq -- it recovers ~40% of PRICE's yield. The model's niche is therefore the cross-project / no-matched-
Ribo-seq scenario (UC-A), and as a complement that surfaces a few ORFs even Ribo-seq misses -- not as a
substitute for a matched Ribo-seq experiment when one exists.

## BOTH released models on all 4 MS datasets: mamba4 vs attn (2026-08-05)

Both models ship (2026-07-31): mamba4 primary/GPU-only, attn the maintained CPU release. The MS
applications were initially run on mamba4 ONLY, which left the whole utility claim resting on the
model a CPU-bound user cannot run. Both are now run on all four datasets. Sources:
`<dataset>_pilot/pgx_{mamba4,attn}/report/table.json`.

The comparison is clean by construction: gencode and the null arms depend only on
(species, universe, mzML, enzyme), not on the model, so the content-addressed cache served
BYTE-IDENTICAL null searches to both models (they completed in 1 s as cache hits). Any difference
below is attributable to the model arm alone.

Novel peptides / DB novel sequences, and density per 1,000 DB sequences:

| dataset | arm | mamba4 | attn | mamba4 dens | attn dens |
|---|---|---|---|--:|--:|
| HBL-1 | model_standard | 8 / 10,345 | 14 / 10,620 | 0.773 | **1.318** |
| HBL-1 | model_poisson | 14 / 1,837 | 12 / 1,819 | **7.621** | 6.597 |
| DoHH2 | model_standard | 32 / 8,198 | 28 / 5,448 | 3.903 | **5.140** |
| DoHH2 | model_poisson | 14 / 1,955 | 16 / 1,796 | 7.161 | **8.909** |
| SU-DHL-4 | model_standard | 22 / 10,351 | 17 / 7,072 | 2.125 | **2.404** |
| SU-DHL-4 | model_poisson | 10 / 1,621 | 13 / 2,948 | **6.169** | 4.410 |
| A549 | model_standard | 9 / 10,312 | 8 / 10,677 | **0.873** | 0.749 |
| A549 | model_poisson | 11 / 2,404 | 7 / 2,296 | **4.576** | 3.049 |

**The two models are indistinguishable for this purpose.** Poisson arm: mamba4 median density 7.161
vs attn 6.597, and the sign flips by dataset -- mamba4 wins HBL-1, SU-DHL-4 and A549, attn wins
DoHH2. Totals across all four: 49 vs 48 novel peptides from 7,817 vs 8,859 database sequences.
Standard arm: 71 vs 67 peptides, and attn has the HIGHER median density (2.404 vs 2.125).

This reproduces the macrophage finding ("more ORFs, not better ORFs", density 26.1 vs 26.5) on a
second substrate class and on a deep TMT proteome. mamba4's +0.0204 held-out profile Pearson
advantage does NOT translate into better ORF selection for MS discovery.

**Consequence for the paper:** the accessibility framing is supported. A CPU-only user running attn
loses nothing measurable in proteogenomic yield. Report mamba4 as primary on the profile-prediction
metrics where it genuinely wins, and state explicitly that the MS application is model-agnostic
across the two releases -- do not imply the GPU model is required for the discovery result.

CAVEAT: n=4 datasets, single-digit-to-low-double-digit peptide counts per arm. This supports
"indistinguishable", which is the weaker and safer claim; it is not powered to detect a small real
difference in either direction.

## UC-A CURRENT: A549 TMT whole-proteome, released mamba4 + pgx + frozen params (2026-08-04)

Supersedes the 2.00x / 2.31x model-vs-null ratios, which came from the old checkpoint, the
f0-threshold DB builder, and unfrozen `calibrate_mass = 2` params (the A549 template's canonical
baseline diverged from both other arms -- see methods.md). Source:
`A549_pilot/pgx_mamba4/report/table.md`. Tryptic, so this run carries the FULL four-arm design
including `null_nc`, which is unusable under nonspecific digestion.

Canonical baseline: 89,203 GENCODE PSMs / 65,405 peptides -- by far the deepest dataset in the
project (the immunopeptidomes are ~1,400-5,500 PSMs), so this is where the design comparison has the
most statistical room.

| arm | novel pep | DB novel seqs | density /1k | vs null_atg | vs null_nc | dPSM | ncStart |
|---|--:|--:|--:|--:|--:|--:|--:|
| null_atg | 10 | 368,908 | 0.027 | 1x | 3x | -4,109 | 0 |
| null_nc | 38 | 3,914,177 | 0.010 | 0.4x | 1x | -11,278 | 6 |
| model_standard | 9 | 10,312 | 0.873 | 32x | 90x | -265 | 0 |
| **model_poisson** | 11 | 2,404 | **4.576** | **169x** | **471x** | **+0** | 1 |

**The canonical-cost column is the cleanest result in the project.** Novel peptides bought per
canonical PSM destroyed:

| arm | canonical PSMs lost per novel peptide found |
|---|--:|
| null_atg | 411 |
| null_nc | 297 |
| model_standard | 29 |
| **model_poisson** | **0** |

A 2,404-sequence model database finds 11 novel peptides at **zero** canonical cost. The 3.9M-sequence
near-cognate null finds 38 -- 3.5x more raw -- but destroys 11,278 canonical PSMs (12.6% of the
baseline) to do it. That is the database-design argument in one row, and it is far more legible here
than on the shallow immunopeptidomes where dPSM sits in the tens.

**The nulls win on raw count and lose catastrophically on everything else.** Same pattern as
everywhere else in this project: never compare raw novel counts across arms of different size.

**Calibration matters more here than anywhere else**: model_poisson beats model_standard 5.2x on
density (4.576 vs 0.873) and eliminates the residual -265 canonical cost entirely.

**ncStart**: null_nc finds 6 peptides absent from the AUG-only null (expected -- it enumerates
near-cognates by construction); model_poisson finds 1. `absent vs nc` is 0 for every arm: nothing any
arm found lies outside the near-cognate null's amino-acid space.

## UC-B CURRENT: three immunopeptidomes, released mamba4 + pgx + frozen params (2026-08-04)

**This supersedes the 2026-07-21 section below**, which used the old checkpoint, the f0-threshold
DB builder, and unfrozen search params. Method now matches the macrophage work exactly: released
`orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1`, Poisson-calibrated RiboCode calls on the predicted
signal, frozen nonspecific template. Source: `<LINE>_pilot/pgx_mamba4/report/table.md`.

Novel columns are 1% CLASS-SPECIFIC FDR (novel targets vs `REV_nuORF|` decoys only); GENCODE and
dPSM columns are 1% GLOBAL FDR against the gencode baseline.

| line | arm | novel pep (uniq) | DB novel seqs | **density** /1k | dPSM | ncStart vs atg |
|---|---|--:|--:|--:|--:|--:|
| HBL-1 | null_atg | 20 | 276,698 | 0.07 | -369 | 0 |
| HBL-1 | model_standard | 8 | 10,345 | 0.77 | -71 | 0 |
| HBL-1 | **model_poisson** | 14 | 1,837 | **7.62** | **+0** | 0 |
| DoHH2 | null_atg | 13 | 242,512 | 0.05 | +13 | 0 |
| DoHH2 | model_standard | 32 | 8,198 | 3.90 | +140 | 1 |
| DoHH2 | **model_poisson** | 14 | 1,955 | **7.16** | +18 | 0 |
| SU-DHL-4 | null_atg | 28 | 259,701 | 0.11 | +43 | 0 |
| SU-DHL-4 | model_standard | 22 | 10,351 | 2.13 | +0 | 1 |
| SU-DHL-4 | **model_poisson** | 10 | 1,621 | **6.17** | +0 | 0 |

**1. The density result is now consistent, which it never was before.** model_poisson vs null_atg:
**105x / 134x / 57x**. Compare the superseded f0 numbers -- 15.20x / 1.95x / 2.37x -- where one
dataset was 8x the others and the "pattern" was really one outlier. Per-1k density itself is tight
across lines (7.62 / 7.16 / 6.17) where the old method gave 0.35 / 0.30 / 0.27 against a null that
moved by 7x. This is the difference the method change makes, and it is the number for Fig 2b.

**2. The Poisson arm costs essentially nothing in canonical IDs** (dPSM +0 / +18 / +0) while
null_atg costs -369 on HBL-1. A 1,837-sequence database barely perturbs the FDR; a 276,698-sequence
one does. That is the DB-design argument made directly.

**3. The null still wins on RAW peptide count in 2 of 3** (20 vs 14, 13 vs 14, 28 vs 10) because it
is ~140x larger. Unchanged from the old method and still the honest framing: the claim is discovery
efficiency per sequence, never raw yield.

**4. Calibration is doing real work.** model_poisson beats model_standard on density everywhere
(7.62 vs 0.77, 7.16 vs 3.90, 6.17 vs 2.13) despite finding FEWER raw peptides on two of three lines.
Both arms are reported per the standing two-arm rule so this is visible rather than asserted.

**5. ncStart is 1 for model_standard on DoHH2 and SU-DHL-4, 0 elsewhere.** Small, and honestly
reported as small. `absent vs nc` is 0 throughout: nothing the model found lies outside the
near-cognate null's amino-acid space -- that null just needs 2.9M sequences to cover it, versus
1,837, and is unsearchable under nonspecific digestion (see methods.md).

CAVEAT: n=3 datasets, single-digit-to-low-double-digit peptide counts. The density separation is
large and consistent, but the absolute counts are small and the panel is half the planned six
(decision D2). Do not over-read any single line.

## UC-B SUPERSEDED: three immunopeptidomes, MS2Rescore-rescored, mm1 coverage (2026-07-21)

HBL-1 was extended to DoHH2 and SU-DHL-4 (locked decision D2; the remaining three of the six planned
datasets -- B721.221, THP-1, one mouse -- are not yet acquired). All three re-run end-to-end on mm1
coverage, novel peptides at 1% class-specific FDR on the rescored mokapot score. Sources:
`data/<line>_pilot/db_comparison_rescored_mm1.md`.

| dataset | model ORFs | model pep | rate | null ORFs | null pep | rate | model/null | overlap (model-only / shared / null-only) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| HBL-1 | 72,439 | **25** | 0.35 | 220,190 | 5 | 0.02 | **15.20x** | 23 / 2 / 3 |
| DoHH2 | 60,814 | 18 | 0.30 | 190,619 | **29** | 0.15 | 1.95x | 1 / 17 / 12 |
| SU-DHL-4 | 67,670 | 18 | 0.27 | 205,041 | **23** | 0.11 | 2.37x | 1 / 17 / 6 |

**Read this carefully -- it is weaker than the HBL-1 number alone suggests.**

1. **The model's per-ORF rate is consistently better** (0.27-0.35 vs 0.02-0.15) and its rate is the
   stable quantity across datasets. That is the defensible claim.
2. **On RAW novel peptides the null WINS in two of three** (29 vs 18, 23 vs 18). The model only leads
   on raw count in HBL-1. Any headline phrased as "the model finds more novel peptides" is false for
   DoHH2 and SU-DHL-4; it finds FEWER, from a database ~3x smaller.
3. **The 15.20x ratio is an outlier, not the pattern.** It comes from the null finding only 5 peptides
   in HBL-1 against 29 and 23 elsewhere. The reproducible effect is ~2x, not ~15x. Quoting 15.20x as
   the headline would not survive review.
4. **Selection behaviour differs qualitatively by dataset.** In HBL-1 the model finds 23 peptides the
   null misses; in DoHH2 and SU-DHL-4 it finds exactly 1, and is otherwise a near-subset of the null
   (shared 17 of 18 in both). So "the model surfaces ORFs the null cannot reach" holds for HBL-1 and
   essentially fails for the other two.
5. Canonical churn is consistently lower for the model (+1.3% to +3.9%) than the null (+2.3% to +5.1%),
   which is the size effect and is consistent across all three.

**CAVEAT 1 (blocking for Fig 2): these ORFs came from the OLD model, not either released model.**
`dump_line.sbatch` and `dump_hbl1.sbatch` pin
`RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes` with the bare (heuristic-Kozak)
`orf_track_v2.npy`. The released models are
`orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` (primary) and the matching
`orf_v2_attn_...` (CPU). The gap is five changes: architecture, union universe, Brain dropped,
`--kozak none`, and mm1 training coverage. The `_mm1` in the comparison filenames refers to the
COVERAGE INPUT at inference, not the checkpoint.

Note the pairing IS internally consistent -- the old model was trained on the heuristic Kozak track
and is fed the heuristic track at inference (`feedback_kozak_never_set`), so this is a stale model,
not a train/inference mismatch. And Task 25 showed mm20 vs mm1 coverage is near-identical on DoHH2
(85.7% of transcripts Pearson >= 0.999), so the coverage-posture part of the gap is small.

By contrast the macrophage results (12 populations) DO use both current release models. So the two MS
applications in the paper currently sit on different model generations, which is not defensible in a
figure. Re-run A549 + the three immunopeptidomes on the released models.

**CAVEAT 2 (CONFIRMED 2026-08-04, blocking): the arms were not scored under the same rules.**
`fragger_hbl1_hla.params` sets `calibrate_mass = 2`, and a divergence test over the full fraction set
shows the optimizer picking DIFFERENT parameters for the null database:

| DB | fragment tol | use_topN | intensity_transform | cutoff @1% FDR |
|---|---|---|---|---|
| canonical | 7 ppm | 150 | 1 | 18.415 |
| model | 7 ppm | 150 | 1 | 18.483 |
| **null** | 7 ppm | **100** | **0** | **16.742** |

The model-vs-null comparison in the table above IS the claim, and its two arms used different peak
selection and intensity handling, with the null's score cutoff 1.67 lower. So the 15.20x / 1.95x /
2.37x ratios are not measuring database quality alone.

I predicted this effect might be small here because 0.02 Da looked like a sensible Orbitrap HCD
tolerance. That was wrong twice over: the optimizer moved it to 7 ppm (0.02 Da is ~20 ppm at
m/z 1000, ~3x too wide), and the divergence that matters showed up in `use_topN_peaks` and
`intensity_transform` rather than the tolerance.

Fixed template: `msfragger/fragger_hbl1_hla_frozen.params`, frozen to the canonical DB's optimizer
output. Re-run required.

## MS2Rescore on A549 -- fully investigated, does NOT help (clean negative result)

MS2Rescore was fanned out across all 12 fractions x 3 DBs (basic + ms2pip CID; DeepLC dropped -- the TMT10
fixed label is unencodable and empties its RT calibration). It runs correctly (+8,639 PSMs/fraction at global
1% FDR in the test), but it does NOT improve novel/cryptic-ORF peptide detection on this tryptic whole-
proteome, and applied naively it hurts:

| approach | model novel | null novel | model/null | why |
|---|---:|---:|---:|---|
| **raw hyperscore class-FDR** | **10** | **13** | **2.31x** | the correct method here (stands) |
| global MS2Rescore + class-FDR | 1 | 6 | 0.50x | mokapot interleaves novel targets/decoys |
| class-aware MS2Rescore | untrainable | untrainable | -- | novel class too sparse to train |

Root cause (structural, not a bug): on a tryptic whole-proteome cryptic-ORF peptides are ultra-rare (~10
novel among ~300k PSMs). MS2Rescore/mokapot is semi-supervised and needs a populated target class.
- **Global rescoring** optimizes the canonical-dominated target/decoy separation and discards the raw
  spectral-match signal (hyperscore) that distinguished the few real novel peptides. Confirmed by the top-25
  novel PSMs: by hyperscore they are `TTTTTTTTTD...` (targets on top), by mokapot score `TDDDDTTTDD...`
  (targets/decoys interleaved) -- so the class-specific FDR collapses 10 -> 1. (My class-FDR logic is sound:
  re-running it on the hyperscore feature reproduces 9 ~= the raw 10.)
- **Class-aware rescoring** (mokapot re-trained on novel_t vs novel_d only, the methodologically correct fix)
  cannot even bootstrap: with ~10-25 real targets among ~35k noise + ~26k decoys, mokapot's (decoys+1) FDR
  finds no confident training set and the SVM folds contain a single class -- untrainable at train_fdr 0.05
  / 0.1 / 0.25 for both model and null.

This is exactly why MS2Rescore is valuable in immunopeptidomics (cryptic peptides abundant) but adds nothing
on a tryptic whole-proteome. The A549 headline stands on raw-hyperscore class-specific FDR (model 10 / null
13, 2.31x). Details: db_comparison_rescored.md (global) + db_comparison_rescored_classaware.md (class-aware).

## MS2Rescore on the HBL-1 immunopeptidome -- WORKS (the opposite of A549)

Fanned MS2Rescore across the 8 immunopeptidome pepXML (4 DBs x 2 reps; basic + ms2pip Immuno-HCD; DeepLC
dropped -- couldn't calibrate on the per-rep non-tryptic search). Here the novel class is populated enough
(non-tryptic, cryptic-enriched) that mokapot's rescored score separates novel targets from decoys, so the
class-specific FDR on the rescored score works and rescoring HELPS:

| DB | novel ORFs | raw novel | rescored novel | rescored rate (pep/1k ORFs) | canonical churn |
|----|---:|---:|---:|---:|---:|
| model | 70,843 | 10 | **26** | 0.37 | -0.1% |
| null | 220,190 | 18 | **4** | 0.02 | -4.5% |
| PRICE (authors') | 136,238 | 25 | **91** | 0.67 | -4.9% |

- **model/null: 1.73x (raw) -> ~20x rate (rescored).** Rescoring cleaned the naive null (18 -> 4: ms2pip
  demotes spurious matches in the huge 3-frame background) and recovered the model (10 -> 26: real novel
  peptides with good fragment match). The model's precise translation-informed DB gives dramatically cleaner
  FDR-controlled discovery; its edge over the naive background widens ~10x.
- **PRICE still leads (91)** but it is the authors' DB from the abandoned low-quality Ribo-seq -- caveated.
- Caveats: 2 reps, small absolute counts; the 20x is partly the class-FDR DB-size/decoy dependence (rate
  ratio ~18x is fairer, still large); ms2pip-only (DeepLC unavailable) so this is a floor on the lift.

**Cross-substrate conclusion:** MS2Rescore adds nothing on the A549 tryptic whole-proteome (cryptic class too
sparse) but substantially sharpens the immunopeptidome comparison, widening the model's advantage over the
naive null -- textbook behavior for where semi-supervised rescoring helps. Details: HBL1_pilot/
db_comparison_rescored.md.

## Macrophage cross-subtype, 12 populations x 6 arms, FROZEN search params (2026-08-03)

Supersedes the earlier run of the same design. That run used `calibrate_mass = 2`, which
re-derives six search parameters per database and made arms non-comparable -- a +699-sequence arm
appeared to *gain* +11,017 canonical PSMs over the GENCODE-only baseline. See methods.md
"Search parameters are FROZEN". Searches: `scripts/pgx/frozen_multiarm.sbatch`, job 36462358,
output `pgx_xsubtype/frozen/`, table `pgx_xsubtype/frozen_reports/crosssubtype_frozen.md`.

Medians across the 12 populations, all columns 1% FDR filtered (novel class-specific vs
`REV_nuORF|` decoys, canonical global):

| metric (median) | model predicted | riboNT | riboALL | null ATG | null near-cognate |
|---|--:|--:|--:|--:|--:|
| novel peptides (unique) | 26 | 27 | 32 | 39 | 126 |
| novel PSMs | 80 | 103 | 139 | 127 | 447 |
| DB novel sequences | 1,578 | 699 | 2,173 | 245,311 | 2,615,086 |
| **discovery density** (pep/1,000 seqs) | 14.0 | **38.6** | 15.0 | 0.1 | 0.0 |
| dPSM vs GENCODE-only | -140 | -70 | -171 | -8,493 | -27,549 |
| ncStart vs null_ATG | **2** | 0 | 0 | -- | -- |

**dPSM is now interpretable.** Every arm costs canonical identifications, and the cost scales
monotonically with database size (-70 for 699 sequences up to -27,549 for 2.6M). That is what
target-decoy competition looks like. The two positive outliers (LargeIntestinal riboNT +34,
SpleenRecruited riboNT +30) are ~0.015% of a ~214,000 baseline, i.e. noise around zero. Under
`calibrate_mass = 2` this column was positive and up to +11,017, which was an artifact.

**The discovery conclusions survived the re-run**, which is the useful cross-check: riboNT density
31.5 -> 38.6, model predicted 17.3 -> 14.0, riboALL 14.5 -> 15.0, null_ATG 0.15 -> 0.1, null_nc
0.05 -> 0.0. Raw novel-peptide medians moved by <=5. The three-orders-of-magnitude gap between the
real-Ribo-seq arms and the near-cognate null was never at risk from a tolerance change.

**Nulls win on raw counts and lose on density, by construction.** null_nc finds the most novel
peptides (median 126) because it is 3,700x larger than riboNT; per sequence it is the worst arm by
two to three orders of magnitude. Never compare raw novel counts across arms of different size.

**ncStart remains the model-only capability**: model predicted finds a median of 2 (max 3)
peptides that the AUG-only null cannot reach. Both real-Ribo-seq arms are structurally 0 because
RiboCode reports N-terminal extensions without testing whether the upstream start is used, so
those databases are ATG-only. That is a property of the caller, not evidence about non-AUG
initiation.

**One asymmetry to read carefully:** `absent_from_null_nc` is 0 for model predicted in all 12
populations but 10 (riboNT) / 14 (riboALL) at median. This is not a quality signal. The
model-predicted database is built from each population's OWN expressed universe, so it is fully
contained in that population's null; the Ribo-seq databases are BMDM-derived and shared across
rows, so some of their ORFs sit on transcripts below the 1 TPM cut in the test population and
therefore fall outside that population's null. It measures universe mismatch, not novelty.

Caveat: the frozen tolerance (200 ppm) was derived on Microglia against GENCODE and applied to all
12 populations on the assumption of one instrument and one study. Consistent per-population
behaviour is compatible with that; it has not been re-derived per population.

## Status

- Macrophage cross-subtype 12 x 6 under frozen params: DONE (supersedes the calibrate_mass=2 run).
- UC-A A549 in-universe (2.00x) + fresh (2.31x): DONE.
- UC-B HBL-1 immunopeptidome four-way (model/null/PRICE): DONE -- PRICE (matched Ribo-seq) wins; model beats
  null, recovers 40% of PRICE + 4 unique, requires no Ribo-seq.
- UC-B extended to DoHH2 + SU-DHL-4 (3 of the 6 planned datasets): DONE. Model leads on per-ORF RATE in
  all three, but loses on RAW novel-peptide count in two of three -- see the honest read above.
- MS2Rescore path validated on A549 (single fraction); full fan-out optional.

### Open / blocking

- Immunopeptidome panel is 3 of 6 (missing B721.221, THP-1, one mouse) -- FIGURES_PLAN flags Fig 2's
  statistical power as the key weakness. The density separation is now large and consistent
  (105x / 134x / 57x), but on 3 datasets with single-digit-to-low-double-digit peptide counts.
- No CPAT/CPC2 coding-potential baseline yet (locked decision D3). Still the biggest scientific gap:
  "beats a naive AUG enumeration" is much weaker than "beats CPAT-selected ORFs".
- The PRICE comparison (UC-B original) has not been re-run on the released model.

### Resolved since 2026-08-03

- `calibrate_mass = 2` verified divergent on all three templates; frozen templates now exist and are
  the pgx default.
- Immunopeptidomes moved from the OLD checkpoint + f0 threshold to the released mamba4 + pgx
  RiboCode-called + Poisson (DBs ~44x smaller, method now identical to the macrophage work).
- `null_nc` under nonspecific digestion established as infeasible (2.86e9 peptides vs MSFragger's
  ~2e9 cap); `null_atg` is the immunopeptidome null and `pgx.run` enforces it.
- **A549 (UC-A) re-run complete** on the released mamba4 + pgx + frozen TMT template, full four-arm
  design. model_poisson: 169x null_atg density at ZERO canonical cost. Supersedes 2.00x / 2.31x.
- **Figure D13 regenerated** from pgx results, both models, 4 datasets:
  `figures/D13_discovery_errorbars/D13_discovery_tradeoff_pgx.png`. The f0-era version is retained
  and labelled SUPERSEDED in its FIGURE_DATA_INPUTS.md.
- **Both released models now run on all 4 MS datasets.** mamba4 and attn are indistinguishable for
  proteogenomic discovery (49 vs 48 novel peptides total, sign of the density difference flips by
  dataset), reproducing the macrophage result. The CPU release is a genuine equivalent here.
