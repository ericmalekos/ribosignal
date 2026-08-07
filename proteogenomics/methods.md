# Proteogenomics pilot -- methods

Goal: test whether the per-nucleotide Ribo-seq P-site model (trained on Chothani matched-tissue
Fibroblast data, GENCODE v49 / GRCh38) is useful as a *translation prior for proteogenomic search-space
selection*. The claim under test: ORFs the model calls translated, scored directly from a cell type's own
RNA-seq (no Ribo-seq for the query sample), yield novel peptides at a higher rate than a naive all-candidate
(3-frame-like) background, and competitively with an experimental Ribo-seq ORF DB, while stealing fewer
canonical IDs.

All external data is logged in `DATA_PROVENANCE.md` before use (project rule). No Ribo-seq is used for any
query sample here: the model predicts translation from sequence + RNA-seq coverage alone.

## Two use cases

- **UC-A (cross-project, A549 whole-proteome).** Deep A549 raw RNA-seq FASTQ (ENCODE ENCSR000CON,
  a DIFFERENT project than the proteome) -> model translation calls -> search the CCLE A549 TMT10
  whole-proteome (MassIVE MSV000085836, Nusinow 2020). RNA and proteome are the same cell line from
  unrelated projects = the realistic unmatched-omics scenario.
- **UC-B (matched, HBL-1 immunopeptidome).** A single DLBCL study (Ruiz-Cuevas 2021, PXD020620 +
  PRJNA647736): same-sample raw RNA-seq FASTQ + a deep MHC-I immunopeptidome + the authors' own PRICE
  Ribo-seq ORF DB. The immunopeptidome is the higher-yield cryptic-ORF substrate (non-tryptic 8-14mers,
  label-free), and the PRICE DB gives a direct experimental-Ribo-seq comparison point the model must match.

## Pipeline (identical for both use cases except the search params)

1. **RNA-seq -> two derived tracks** (same cell line, one alignment):
   - **STAR** (GRCh38 GENCODE v49 index `/private/groups/.../STAR_indexes/star_index_grch38_v49`) ->
     per-nt transcriptome coverage hd5 on the model's universe (the model's RNA-seq input channel).
   - **salmon** (decoy-aware `/private/groups/.../Salmon_indexes/salmon_index_decoy_v49`) -> per-tx TPM
     (defines the EXPRESSED universe; salmon always decoy-aware per project rule).
   RNA-seq strandedness confirmed per dataset before alignment (A549 ENCSR000CON = reverse/ISR).
2. **Universe = the cell type's OWN expressed transcriptome** ("fresh universe"): every tx with TPM above
   threshold, capped at <= 10,000 nt (attention memory; the one-hot deployment model is O(L^2) in the
   mixer). Predicting on the cell type's own expressed tx, not the Fibroblast training universe, is the
   point of a cell-type-specific translation prior -- A549/HBL-1-specific transcripts (and their novel
   ORFs) are otherwise skipped. A coverage-only pack (`build_heldout_pack.py`, `coverage_only=True`) carries
   zeros for the (absent) Ribo-seq target so the pack format matches training.
3. **Model prediction** (`dump_pred_profiles.py --tx_list`): the deployed one-hot `orf_v2_attn` model emits
   a per-nt predicted P-site profile for every universe tx from sequence + RNA-seq coverage (no Ribo-seq).
4. **ORF enumeration + scoring** (`enumerate_score_orfs.py`): enumerate all AUG..stop ORFs per tx, classify
   (canonical / uORF / dORF / lncRNA_orf / nterm_ext / internal), score each by `pred_frame0` = fraction of
   the ORF's predicted P-site mass in frame 0 (in-frame periodicity the model predicts). Emits
   `candidates.faa` (headers carry class + f0).
5. **DB construction** (`build_a549_dbs.py`, `build_price_db.py`): three-to-four search DBs sharing an
   IDENTICAL canonical GENCODE v49 base (160,646 unique proteins >= 7 aa) + inline REV_ decoys, differing
   ONLY in the novel-ORF set, so the search isolates ORF-SET provenance:
   - `db_canonical` : canonical only (baseline).
   - `db_model`     : canonical + model-SELECTED novel ORFs (pred_frame0 >= 0.5).
   - `db_null`      : canonical + ALL novel candidates (naive 3-frame-like background).
   - `db_price`     : canonical + the authors' PRICE Ribo-seq ORFs (UC-B only).
   Novel proteins are deduplicated by sequence, dropped if identical to a canonical protein, tagged
   `nuORF|<class>|<id>` (PRICE = `nuORF|price|<id>`) so the class-specific FDR can isolate the novel class.
6. **MSFragger search** (4.2, conda `msfragger` env):
   - UC-A: TMT10 params (`fragger_a549_tmt.params`) -- add_Nterm/K = 229.162932, add_C = 57.021464,
     fragment 0.6 Da (ion-trap MS2), precursor +-20 ppm, tryptic, decoy_prefix REV_.
   - UC-B: HLA params (`fragger_hbl1_hla.params`, adapted from the sibling immunopeptidome params) --
     enzyme nonspecific (num_enzyme_termini 0), digest 8-14 aa, NO TMT (label-free), fragment 0.02 Da
     (Q Exactive HF Orbitrap HCD MS2, high-res), precursor +-10 ppm, variable Ox-M + N-term acetyl.
   NOTE: the `msfragger` conda build ships the Thermo `.raw` reader as `.exe` (needs mono/.NET, absent),
   so HBL-1 `.raw` are converted to indexed centroided mzML with ThermoRawFileParser 2.0 (proteomics env)
   before search. CCLE A549 arrived as CCMS-converted mzML already.
7. **Class-specific FDR + comparison** (`compare_model_vs_null.py`, generalized `compare_dbs.py`): a NOVEL
   peptide is a rank-1 PSM whose peptide maps ONLY to nuORF| ORFs (not shared with any canonical protein);
   its FDR is estimated against REV_nuORF| decoys ONLY (the honest per-class denominator, not the global
   decoy pool -- per the expression_context_human "crazy FDR" lesson: a global FDR inflates novel discovery
   ~10x). Reports per DB: novel peptides at 1% class FDR, discovery rate (pep / 1k ORFs), canonical target
   PSMs (PC-churn -- do novel ORFs displace canonical IDs?), and the novel-peptide SET OVERLAP between DBs
   (model-only vs PRICE-only vs shared = does the model recover ORFs the experimental Ribo-seq DB missed?).

## Optional rescoring (MS2Rescore)

`ms2rescore` (mokapot engine) adds fragment-intensity (ms2pip) and, where usable, retention-time (DeepLC)
features to boost PSM sensitivity at fixed FDR. Configs: `ms2rescore_a549.json` (basic + ms2pip CID; DeepLC
DROPPED -- the TMT10 fixed label mass 229.16 is unencodable by DeepLC and empties its q<=0.01 calibration
set), `ms2rescore_hbl1.json` (basic + ms2pip immuno-HCD + DeepLC; label-free, so DeepLC calibrates fine).
Rescoring reads the kept per-fraction pepXML. This is a sensitivity refinement; the headline model-vs-null
result stands on raw-hyperscore class-specific FDR.

## Search parameters are FROZEN for cross-database comparison (2026-08-03)

`calibrate_mass = 2` makes MSFragger re-derive six search parameters from a first-pass search
against **whichever database it is given**. Arms of the same population are then scored under
different rules, so canonical counts move for reasons unrelated to the sequences added. This was
diagnosed after a +699-sequence arm appeared to gain +11,017 canonical PSMs over GENCODE alone --
not target-decoy competition at all.

What the optimizer was silently changing (read off the `New <param> =` lines in the MSFragger log):

| param | template | optimizer |
|-------|----------|-----------|
| `fragment_mass_tolerance` | 0.6 Da | **200 ppm** |
| `use_topN_peaks` | 150 | 75 |
| `max_fragment_charge` | 3 | 1 |

The tolerance is the load-bearing one. 0.6 Da is the right kind of setting for this low-res
ion-trap MS2 but too wide; 200 ppm is ~0.1 Da at m/z 500 and ~0.2 Da at m/z 1000, and is worth
~4x the identifications. The apparent "mass calibration triples sensitivity" effect was the
optimizer compensating for that one value: `calibrate_mass = 1` (calibration WITHOUT parameter
optimization) gives 2,291 targets against 2,333 for no calibration at all, i.e. calibration
itself contributes nothing here.

Fix: `msfragger/fragger_macro_lfq_frozen.params` fixes all six values and sets
`calibrate_mass = 0`, so every arm is scored identically. Verified on Microglia, 18 fractions,
GENCODE vs riboNT:

| | canonical @1% | dPSM | delta rank1 | novel pep |
|---|---|---|---|---|
| `calibrate_mass = 2` | 205,348 -> 216,365 | **+11,017** | +28,195 | 16 |
| FROZEN | 214,214 -> 214,279 | **+65** | +174 | **28** |

Frozen is *more* sensitive (+8,866 canonical) and finds *more* novel peptides (28 vs 16) while
reducing dPSM 170-fold. Re-derive the values if the instrument or study changes: run one arm with
`calibrate_mass = 2` and read the log.

### Confirmed on ALL THREE templates (2026-08-04)

The macrophage template was not a one-off. Every production template was running
`calibrate_mass = 2`, and a divergence test on the full fraction set of each pilot -- canonical vs
model vs null, `_diag_optdiv.sbatch` -- shows the optimizer landing on different parameters per
database in both remaining use cases:

| pilot | DB | fragment tol | use_topN | intensity_transform | remove_precursor | cutoff @1% |
|---|---|---|---|---|---|---|
| HBL-1 | canonical | 7 ppm | 150 | 1 | 1 | 18.415 |
| HBL-1 | model | 7 ppm | 150 | 1 | 1 | 18.483 |
| HBL-1 | **null** | 7 ppm | **100** | **0** | 1 | **16.742** |
| A549 | canonical | 150 ppm | 100 | **1** | **0** | 19.354 |
| A549 | **model** | 150 ppm | 100 | **0** | **1** | 18.061 |
| A549 | **null** | 150 ppm | 100 | 0 | 1 | 18.347 |

Which arm diverges differs by dataset, and both patterns are damaging:

- **HBL-1: the NULL diverges.** The model-vs-null discovery comparison IS the immunopeptidome claim,
  and the two arms were scored under different rules with the null's cutoff 1.67 lower.
- **A549: the CANONICAL baseline diverges** from both model and null, so the dPSM reference itself
  was not comparable to the arms measured against it.

Same underlying cause each time -- a fragment tolerance too wide for the data, which the optimizer
was silently correcting: 0.6 Da -> 200 ppm (macrophage), 0.02 Da -> 7 ppm (HLA), 0.6 Da -> 150 ppm
(A549 TMT). Three templates, three too-wide tolerances, three silent rescues.

Frozen templates now exist for all three: `fragger_macro_lfq_frozen.params`,
`fragger_hbl1_hla_frozen.params`, `fragger_a549_tmt_frozen.params`. Each is frozen to its CANONICAL
database's optimizer output, since canonical is the baseline every arm is compared against.

**Consequence:** the UC-A ratios (2.00x / 2.31x) and the immunopeptidome ratios
(15.20x / 1.95x / 2.37x) are all confounded and must be re-run frozen. Combine that with the model
re-run -- those results also used the OLD checkpoint -- so the searches happen once, not twice.

Two routes that do NOT work, recorded so they are not retried: `write_calibrated_mzml = 1` emits
MS1 spectra only (6,051 of 77,289 on a test fraction) and is not searchable; `write_mzbin_all = 1`
emits UNCALIBRATED spectra (2,336 vs 2,333 targets against `calibrate_mass = 0`, identical
cutoff) and is a format cache only. The mzBIN cache is still used by
`scripts/pgx/frozen_multiarm.sbatch` for speed -- convert once per population (~10 min), then
each arm searches the cache in ~7 min instead of re-parsing mzML -- and it guarantees every arm
sees byte-identical spectra.

**Consequence for earlier results:** every search predating this ran with per-database
auto-optimized parameters, so canonical counts and dPSM columns are not comparable across arms
and are being re-run. Novel-peptide counts and discovery densities are far more robust (the
riboNT-vs-null_nc gap is three orders of magnitude), but the re-run supersedes them too.

## CORRECTION: a NaN p-value silently zeroed four extension arms (found 2026-08-06)

`pgx.seqtools.bh` sorted p-values with `np.argsort`, which places NaN LAST, then applied a REVERSE
cumulative minimum. Since `np.minimum(NaN, x)` is NaN, the accumulate started on that NaN and
propagated it through every element: **one NaN p-value nulled every q-value in the array**, so
nothing could clear `q <= 0.05`.

The NaN comes from a degenerate Wilcoxon on a flat extension region -- a legitimate edge case, not
bad data. On B721/mamba4/standard: 4,473 tests, 2,785 with p < 0.05, exactly ONE NaN, and all 4,473
q-values returned NaN.

**This is the most dangerous failure shape in the project.** It does not crash or warn. It reports
`passing: 0` beside `would_pass_whole_orf_f0_0.5: 4434`, which reads as a striking finding -- "the
extension-region test rejects everything the crude criterion accepts" -- and is pure arithmetic.
It surfaced only because two models disagreed on the same input (0 vs 4,122 extensions), so an
internal inconsistency was visible.

SCOPE, checked rather than assumed. `bh()` has exactly ONE caller (`pgx/extensions.py`), so only
extension calling is affected: **4 of 48 arms, all `standard` (theta = 1)**. All 44 Poisson arms are
clean, which means every headline result stands (density table, A549 zero-canonical-cost, CPAT/CPC2
comparison, 12-population cross-subtype table -- all use the Poisson arm). The ORF-call evaluations
are also unaffected: `ribocode_dropin.py` passes `pval_adj="fdr_bh"` to RiboCode's own
`detectORF.main`, a separate statsmodels implementation, so drop-in F1 0.931 / 0.929, the 9-fold
LOTO spread, depth crossover and the hepatocyte F1 0.90 never touched this code.

| arm | tested | passing before | after | DB novel seqs before -> after |
|---|--:|--:|--:|---|
| DoHH2 / attn / standard | 5,261 | 0 | 2,840 | 5,448 -> 8,269 |
| SU-DHL-4 / attn / standard | 7,757 | 0 | 3,805 | 7,072 -> 10,854 |
| B721 / mamba4 / standard | 4,473 | 0 | 2,649 | 6,290 -> 8,921 |
| macrophage BMDM / attn_union / standard | 9,788 | 0 | 4,445 | 5,991 -> 10,411 |

Fixed in `bh()` (non-finite p-values are excluded from the correction and get NaN q-values back, so
they can never pass a threshold and never poison neighbours), pinned by `test_bh_is_nan_safe`, and
re-run via `pgx/redo_nan_arms.sh`. Post-fix pass rates sit in the normal band everywhere
(standard 45-61%, poisson 68-85%); no arm reports `tested > 0, passing = 0`.

**The correction makes the calibration argument STRONGER, which is worth stating explicitly rather
than quietly folding in.** The bug had been flattering the uncalibrated arm: a truncated database
pays a smaller FDR penalty. With the correct databases the theta = 1 arms get worse -- SU-DHL-4
model_standard fell from 17 novel peptides to 5 while its database grew 7,072 -> 10,854, so density
dropped 2.404 -> 0.461. model_poisson now beats model_standard on density by 2.0x on DoHH2
(7.16 vs 3.63) and 14.3x on SU-DHL-4 (6.17 vs 0.46), against 1.8x and 2.9x before the fix. That is
the same decoy-load mechanism the model-vs-null comparison rests on, operating within the model arms.

## Immunopeptidome runs: pgx configuration (2026-08-04)

### ORF selection is RiboCode-called, not f0-thresholded

The A549 and immunopeptidome pilots originally selected ORFs by enumerating every candidate in the
universe and keeping `pred_frame0 >= 0.5` (`enumerate_score_orfs.py` + `build_a549_dbs.py`). No
RiboCode, no Poisson calibration, no significance test. That is what made the databases balloon:
80,608 novel ORFs for HBL-1.

They now use the same pgx path as the macrophage work -- Poisson-calibrate the predicted depth, run
RiboCode `detectORF` on the PREDICTED signal, keep calls passing the q-value cut:

| line | f0 threshold | pgx standard (theta=1) | pgx **poisson** |
|---|--:|--:|--:|
| HBL-1 | 80,608 | 10,345 | **1,837** |
| DoHH2 | 65,631 | 8,198 | **1,955** |
| SU-DHL-4 | 72,560 | 10,351 | **1,621** |

~44x smaller than f0, and on the same scale as the macrophage model arm (median 1,578), so the two
MS applications are finally on one method. Calibration selected theta 0.05-0.10.

Extension testing behaves as designed on this substrate: HBL-1 tested 7,770 candidate N-terminal
extensions and 3,509 passed q<=0.05 (45%), against the retired whole-ORF `f0>=0.5` criterion which
would have admitted 7,696 of 7,770 (99%).

### null_nc is NOT usable with nonspecific digestion -- null_atg is the null

Nonspecific digestion emits every 8-14mer substring, so the ~2.5-3M-sequence near-cognate null
generates **2,857,862,951 peptides** against MSFragger's ~2 billion hard cap. All three immunopeptidome
`null_nc` searches died with "Too many peptides were generated". This is an ENZYME limit, not a
database-size limit -- the same database is fine tryptic, which is why the macrophage runs never hit
it.

DECISION: immunopeptidome (nonspecific) runs use **gencode / model_standard / model_poisson /
null_atg**. Splitting the database into sub-2e9 chunks was considered and declined. `pgx.run` now
drops `null_nc` automatically when `--enzyme nonspecific` and says so, so this fails loudly at
configuration time rather than an hour into a search.

**What this does NOT cost: the ncStart columns.** `absent vs nc` / `ncStart vs nc` are computed by
substring-scanning the null's AMINO-ACID SPACE built from `db_null_nc.fasta` -- they never needed the
null to be searched (verified: `null_nc AA space: 236.3 Mb` still loads). Dropping the search costs
only null_nc's own novel-peptide and dPSM row. The question "would a naive near-cognate enumeration
have found this peptide" is still answered.

### Template

`pgx.search.DEFAULT_TEMPLATE` now points at the FROZEN templates for both enzymes
(`fragger_macro_lfq_frozen.params`, `fragger_hbl1_hla_frozen.params`). pgx exists to compare
databases, so an unfrozen template -- which re-derives search parameters per database -- is never the
right default there. `pgx.run --template` overrides it to reproduce a pre-2026-08-04 run.

## Key parameters

| param | UC-A (A549) | UC-B (HBL-1) |
|-------|-------------|--------------|
| RNA-seq | ENCODE ENCSR000CON (cross-project) | SRR12285172 (same-sample) |
| universe cap | <= 10,000 nt | <= 10,000 nt |
| model | one-hot orf_v2_attn (deployed) | same |
| pred_frame0 model-select cut | 0.5 | 0.5 |
| min ORF length | 7 aa | 7 aa |
| MS instrument | Orbitrap Fusion TMT10 SPS-MS3, ion-trap MS2 | Q Exactive HF, Orbitrap HCD MS2 |
| enzyme | trypsin | nonspecific (num_enzyme_termini 0), 8-14 aa |
| labels | TMT10 (fixed 229.16 K/N-term) | label-free |
| fragment tol | 0.6 Da | 0.02 Da |
| decoy | REV_ inline reverse | REV_ inline reverse |
| FDR | 1% class-specific (novel = nuORF|-only) | same |
