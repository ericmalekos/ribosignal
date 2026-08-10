# `pgx` -- model-informed proteogenomics search pipeline

Takes the RNA preprocessing output and the model's predicted per-nt Ribo-seq signal, calls ORFs
with Poisson-calibrated RiboCode, builds matched search databases, runs MSFragger, and reports
FDR-controlled discovery. Works for human and mouse. No hardcoded paths.

```
predicted signal (pred_profiles.npz)
        |
        v
  calibrate  -- Poisson/theta sweep, CDS-anchored          -> theta.json
        |
        v
  call_orfs  -- RiboCode on the predicted signal, TWO arms -> calls_<arm>.tsv
        |
        +-- extensions -- N-terminal CDS extension scan     -> extensions.tsv
        |
        v
  build_dbs  -- gencode | model | null_atg | null_nc
        |
        v
  search     -- MSFragger, tryptic or nonspecific, hash-cached
        |
        v
  report     -- global FDR canonical + class-specific FDR novel -> table.md
```

## Quick start

```bash
cd proteogenomics/scripts
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

$PY -m pgx.run \
  --species mouse --label BMDM \
  --profiles <...>/pred_profiles.npz \
  --salmon <...>/rep1/quant.sf <...>/rep2/quant.sf <...>/rep3/quant.sf \
  --mzml '<...>/mzml/BMDM_Rep*_F*.mzML' \
  --enzyme tryptic \
  --out <...>/pgx/BMDM \
  --submit-searches

# once the MSFragger jobs finish
$PY -m pgx.run --species mouse --label BMDM --profiles <...> --out <...>/pgx/BMDM --report-only
```

Run `pgx.run` itself on a compute node: the calibration sweep is several RiboCode passes.
Every step is skip-if-done, so a rerun resumes rather than repeats.

## The databases

All share one GENCODE proteome base plus inline `REV_` decoys, so they differ ONLY in novel
content and the comparison isolates the selection rule. Every arm uses the SAME expressed universe
(TPM >= 1, protein_coding + lncRNA, no chrM, <= 10 kb).

| database | novel content | role |
|---|---|---|
| `db_gencode` | none | baseline detection; the reference for gained/lost canonical PSMs |
| `db_model_<arm>` | significant RiboCode calls + N-terminal extensions | the model's selection |
| `db_null_atg` | every ATG..stop ORF on the universe | naive AUG-only pipeline |
| `db_null_nc` | every ATG-or-near-cognate ORF, same universe | naive near-cognate pipeline |
| `db_cpat` / `db_cpc2` | the null_atg pool filtered by a coding-potential classifier | sequence-only comparator |

`db_null_atg` is a strict subset of `db_null_nc`, so the increment between them is exactly the
cost of naive non-AUG enumeration. Measured on BMDM: 1,533 / 218,715 / 2,306,644 novel sequences.

**The coding-potential arms are the comparator that makes the model's selection falsifiable.**
`null_atg` only establishes "better than enumerating everything"; CPAT and CPC2 shrink the space by
the same order of magnitude as the model from sequence composition alone, with no RNA-seq and no
translation model. `coding_potential.py` therefore enumerates from **exactly** the `null_atg` pool
(same universe, start codons, minimum length, and the `or prot in canon_seqs` clause) and then
applies each tool's own default human classifier, so the arms differ only in WHICH ORFs are kept.
That equality is pinned by `tests/test_invariants.py::test_coding_potential_pool_matches_null_arm`;
without it a silent drift would turn a selection-rule comparison into a pipeline comparison with no
error raised. Built separately from `run.py` via `coding_potential.sbatch`.

## Why RiboCode instead of a `pred_frame0` threshold

`f0` is a ratio over the ORF interval, so prepending a zero-signal upstream region changes neither
numerator nor denominator: it is *invariant* to N-terminal extension. Measured on BMDM, 96.2% of
51,646 upstream non-AUG ORFs sat within 0.05 of the true CDS's `f0` (median |delta f0| = 0.0026).
No threshold separates them at any value, which is why enabling non-AUG enumeration inflated the
old model DB from 13,835 to 142,248 sequences without adding evidence.

## N-terminal extensions are a separate step, by necessity

RiboCode cannot supply them:

* `orf_finder.orf_find` sets `alt_flag = 0` as soon as any in-frame ATG precedes the stop, so a
  **near-cognate extension of an annotated CDS is unreachable at any parameter setting**.
* For ATG extensions it does emit a call (`classfy_orf` labels same-stop ORFs `annotated`), but
  `start_check` with `only_longest_orf=True` just takes the most 5' start and tests the WHOLE ORF,
  never asking whether that start is used.

`pgx.extensions` therefore applies RiboCode's own frame statistic (`extract_frame` + `test_frame`,
Wilcoxon f0>f1 and f0>f2, Stouffer-combined, then BH) to the **extension region alone** -- the
segment from the candidate start to the annotated CDS start. That is where the two hypotheses
differ: real upstream initiation puts periodic P-sites there, a downstream true start leaves it
empty. On BMDM this admitted 564 of 30,243 candidates, where the whole-ORF `f0` criterion admitted
100% of everything it scored.

## Calibration anchor

theta scales the predicted effective depth. The dial is anchored on **relative** CDS recall
(`hits / hits at the saturating theta`), not absolute recall of the annotation: only ~53% of
expressed CDS-bearing transcripts are detectably translated on BMDM, so an absolute 0.90 target is
unreachable by construction and the dial would always saturate. The relative anchor is
self-calibrating across sequencing depths and reproduced theta = 0.05, the operating point
validated as transferable in `docs/count_head_calibration_checks.md` (Check 5).

## Significance

RiboCode filters per-ORF on `--pval` BEFORE BH adjustment, so on the surviving set BH moves p by
at most 4e-4 and `--orf-qvalue 0.05` is a no-op. **`--pval` is the real significance lever**;
`--orf-qvalue` only tightens further.

## Enzyme modes

| `--enzyme` | settings | use |
|---|---|---|
| `tryptic` | trypsin KR, `num_enzyme_termini = 2`, 1 missed cleavage, 6-50 aa | shotgun proteome |
| `nonspecific` | `num_enzyme_termini = 0`, 8-14 aa | HLA / immunopeptidomics |

Tolerances, modifications and instrument settings are inherited from a validated template in
`proteogenomics/msfragger/`; override with `--template`.

MSFragger 4.2 sizes its own fragment-index slices from the JVM heap and rejects `num_slices`, so
the only lever for a large database is `XMX` in `search.sbatch`. The 2.36M-target near-cognate
null built 672,835,480 fragments as 6.27 GB in a single slice.

## Avoiding duplicated work

* The Poisson arm at theta* is **reused from the calibration sweep**, never recomputed.
* Each search is keyed by `sha256(database) + sha256(settings)` in `<out>/.pgx_hash`; a matching
  digest with outputs present is skipped. The GENCODE baseline is therefore searched once per
  (mzML set, enzyme) and reused by every arm's reporting. **`database_name` is excluded from the
  settings digest**: it holds the FASTA's absolute path, which is per-run, so including it made two
  searches over byte-identical databases hash differently and defeated `--shared-search-root`
  entirely. The database's content is already covered by `sha256(database)`; where the file lives
  is not part of what the search is.
* `--shared-search-root` puts the model-independent arms (`gencode`, `null_atg`, `null_nc`) in a
  common location. Those depend only on (species, universe, mzML, enzyme), so comparing two model
  checkpoints on the same data re-searches only the model arm. This is the expensive one to get
  right: a `null_nc` search over ~2.5M targets runs for hours and is byte-identical across
  checkpoints. Shared arms are symlinked back under `<out>/search` so reporting still sees one
  directory.
* Every step is skip-if-done on its output file.

## Reporting rules

| column group | FDR |
|---|---|
| canonical / dPSM | **global** 1% (all targets vs all decoys) |
| novel | **class-specific** 1% (novel targets vs `REV_nuORF|` decoys only) |

A class-specific FDR on the canonical side is confounded: a large novel space cannibalises
canonical decoys faster than canonical targets, so a junk database can appear to gain canonical
identifications. A global FDR on the novel side inflates non-canonical discovery ~10-13x.

**No double counting against GENCODE.** MSFragger reports every protein containing a peptide, so a
peptide occurring in ANY GENCODE protein is counted as GENCODE, never as novel. The protein list is
`;`-delimited, not `,`; `pgx.report` splits on both, so correctness does not depend on MSFragger
listing the canonical protein first.

`ncStart` counts passing novel peptides absent from a null's amino-acid space (substring scan,
since a peptide is a fragment). Reported against both nulls: vs `atg` answers "would an AUG-only
pipeline ever have found this", vs `nc` answers "would even naive near-cognate enumeration".

## Modules

| module | role |
|---|---|
| `refs.py` | species registry + path resolution (CLI > `PGX_*` env > registry) |
| `universe.py` | salmon quant.sf -> expressed universe tx list + FASTA |
| `rc_io.py` | RiboCode `*_collapsed.txt` reader, coordinate conversion, class mapping |
| `seqtools.py` | FASTA, translation, ORF enumeration, CDS geometry, BH |
| `calibrate.py` | Poisson/theta sweep + CDS anchoring (`sweep_theta.sbatch` fans it out) |
| `call_orfs.py` | two-arm RiboCode calling on the predicted signal |
| `extensions.py` | N-terminal extension scan (runs in the `ribocode` env) |
| `build_dbs.py` | the gencode / model / null databases + class maps |
| `coding_potential.py` | the CPAT and CPC2 comparator arms (same pool as `null_atg`, different filter) |
| `search.py` / `search.sbatch` | MSFragger templating, hash cache, staging |
| `report.py` | FDR control + the output table |
| `run.py` | end-to-end orchestrator |
