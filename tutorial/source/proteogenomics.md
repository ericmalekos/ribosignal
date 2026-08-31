# Proteogenomics: from predicted signal to peptide evidence

The ORF calls are a prediction. Mass spectrometry is the independent test: if a predicted non-canonical
ORF is really translated, its peptides should be in the proteome. The `pgx` pipeline
(`proteogenomics/scripts/pgx/`) turns a predicted profile into a search database, runs MSFragger, and
reports FDR-controlled discovery against matched naive baselines.

```{mermaid}
flowchart LR
    P["Predicted profile<br/><small>pred_profiles.npz</small>"]
    C["Poisson calibration<br/><small>theta*</small>"]
    R["RiboCode calls<br/><small>two arms</small>"]
    E["Extension scan<br/><small>non-AUG starts</small>"]
    D["Search databases<br/><small>model / null / GENCODE</small>"]
    M["MSFragger<br/><small>tryptic or non-specific</small>"]
    T["FDR table"]
    P --> C --> R --> D --> M --> T
    R --> E --> D
    classDef s fill:#eaf2fb,stroke:#3a6ea5,color:#16314a;
    class P,C,R,E,D,M,T s;
```

## The databases

Every database shares one canonical base (the GENCODE proteome) plus inline reversed `REV_` decoys, so the
arms differ **only** in novel content and the comparison isolates the selection rule. Both nulls use the same
expressed universe as the model arm (TPM >= 1, protein_coding + lncRNA, no chrM, <= 10 kb).

```{list-table}
:header-rows: 1
:widths: 22 14 64

* - Database
  - Novel seqs
  - What it represents
* - `db_gencode`
  - 0
  - baseline detection; the reference for gained/lost canonical PSMs
* - `db_model_poisson`
  - 1,533
  - calibrated calls ($\theta$*) + extensions
* - `db_model_standard`
  - 10,454
  - uncalibrated calls ($\theta$=1) + extensions
* - `db_null_atg`
  - 235,686
  - every ATG..stop ORF on the universe -- the naive AUG pipeline
* - `db_null_nc`
  - 2,512,388
  - every ATG-or-near-cognate ORF -- the naive non-AUG pipeline
* - `db_cpat` / `db_cpc2`
  - 7.6k-21k
  - the `null_atg` pool filtered by a coding-potential classifier
```

`db_null_atg` is a strict subset of `db_null_nc`, so the increment between them is exactly the cost of naive
non-AUG enumeration. Numbers are BMDM, mouse, `mamba4` union model.

### The coding-potential arms are what make the selection rule falsifiable

Beating `null_atg` only establishes "better than enumerating everything". CPAT and CPC2 are the standard
sequence-only coding-potential selectors, and they shrink the search space by the *same order of magnitude* as
the model using sequence composition alone -- no RNA-seq, no translation model. That makes them the honest
comparator.

`pgx/coding_potential.py` enumerates from **exactly** the `null_atg` pool (same universe, start codons,
minimum length, and the `or prot in canon_seqs` clause), then applies each tool's own default human
classifier. So the arms differ only in *which* ORFs are kept. That equality is pinned by
`tests/test_invariants.py::test_coding_potential_pool_matches_null_arm` -- without it a silent drift would
turn a selection-rule comparison into a pipeline comparison, with no error raised.

**The result splits by assay, and the split is the finding.** Novel peptides at 1% class-specific FDR:

| dataset | substrate | CPAT | CPC2 | model $\theta$=1 | model Poisson |
|---|---|--:|--:|--:|--:|
| A549 | tryptic | **25** | 23 | 9 | 11 |
| HBL-1 | HLA-I | 3 | 10 | 8 | **14** |
| SU-DHL-4 | HLA-I | 7 | 14 | **22** | 10 |
| DoHH2 | HLA-I | 4 | 8 | **32** | 14 |

CPAT/CPC2 win the tryptic whole proteome; the model wins all three HLA-I immunopeptidomes. The mechanism is
visible in what each keeps: coding-potential tools score a *transcript's* composition, so they retain long,
codon-biased, ORF-like sequences -- exactly what a tryptic digest samples well, and the least novel population
available. The model scores *per-nucleotide translation* from the sample's own RNA-seq, so it retains short,
non-canonical ORFs that are translated in **this** cell type, which is what HLA-I presentation samples.

```{admonition} Where the cell-type specificity actually comes from
:class: important
The input ablation (Task 61) measured this rather than assuming it, and the answer is narrower than
"the model is cell-type-specific". Zeroing the RNA-seq channel and retraining the deployed recipe
costs almost nothing in profile SHAPE -- sequence alone recovers 98.5% of it (pc profile Pearson
0.6601 of 0.6699) -- but costs a lot in MAGNITUDE: count Pearson falls 0.8990 to 0.7063.

So the specificity runs through the **count head**. A transcript not expressed in the query cell type
gets low predicted depth, fails the caller's significance test, and never enters the database. That is
a real cell-type-specific mechanism for deciding *which* ORFs are called, and it is what this
comparison rests on. It is not a claim that the predicted profile shape differs between cell types --
it largely does not, and that should not be written anywhere.
```

```{figure} /img/B6_input_ablation/B6_input_ablation.png
:alt: Input-modality ablation on the deployed recipe
:width: 100%

Retraining the deployed recipe with one input zeroed. **Left:** zeroing RNA-seq costs almost nothing
in profile shape (0.660 of 0.670, 98.5%), and RNA-seq alone cannot produce a shape at all (0.123).
**Right:** the count head loses 0.193 Pearson without RNA-seq, and RNA-seq alone (0.808) beats
sequence alone (0.706) at per-transcript depth. Generator `figures/B6_input_ablation/`.
```

Figure `figures/D12_cpat_cpc2/`. Caveat: single-digit to low-double-digit counts throughout; the direction is
consistent across three independent immunopeptidomes but no single dataset carries the claim.

## Why the ORF caller, and not a score threshold

The first version of this database selected ORFs with a threshold on `pred_frame0`, the in-frame fraction of
predicted density. That criterion **cannot work for N-terminal extensions**, and the reason is arithmetic
rather than empirical:

$$f_0 = \frac{\sum_{i \equiv 0 (3)} \hat{y}_i}{\sum_i \hat{y}_i}$$

is a **ratio** over the ORF interval. Prepending a zero-signal upstream region adds nothing to either the
numerator or the denominator, so $f_0$ is **invariant** to extension, not merely diluted. An upstream CUG
that happens to be in frame with a translated CDS inherits that CDS's score exactly.

Measured on BMDM over 12,719 stop-codon groups containing an annotated CDS start plus upstream in-frame
starts: **96.2%** of 51,646 upstream non-AUG ORFs sat within 0.05 of the true CDS's $f_0$, with a median
$|\Delta f_0|$ of **0.0026**. No threshold separates them at any value.

```{admonition} The symptom that exposed it
:class: warning
Enabling non-AUG enumeration multiplied the candidate pool ~10x (625,698 -> 5,956,246 ORFs). The $f_0$ cut
passed the same ~7.4% of it, inflating the model database from 13,835 to 142,248 sequences without adding
any evidence. A score that measures reading frame cannot answer a question about start position.
```

## N-terminal extensions need their own test

RiboCode cannot supply these, for two independent reasons in its source:

- `orf_finder.orf_find` sets `alt_flag = 0` as soon as **any** in-frame ATG precedes the stop, then skips the
  alternative-start list entirely. An annotated CDS always has an ATG, so a **near-cognate extension of an
  annotated CDS is unreachable at any parameter setting**.
- For ATG extensions it does emit a call (`classfy_orf` labels same-stop ORFs `annotated`, and `start_check`
  with `only_longest_orf=True` takes the most 5' start), but it tests the **whole ORF** and so never asks
  whether that start is used.

`pgx.extensions` applies RiboCode's own statistic (`extract_frame` + `test_frame`: Wilcoxon $f_0 > f_1$ and
$f_0 > f_2$, Stouffer-combined, then BH) to the **extension region alone** -- the segment from the candidate
start to the annotated CDS start. That is where the two hypotheses differ: real upstream initiation puts
periodic P-sites there, a downstream true start leaves it empty.

```{list-table}
:header-rows: 1
:widths: 60 20 20

* - BMDM, calibrated arm
  - n
  - 
* - candidate upstream in-frame starts
  - 30,243
  - 
* - dropped: too short
  - 5,396
  - 
* - dropped: < 5 nonzero frame-0 codons in the extension region
  - 23,963
  - 
* - testable
  - 884
  - 
* - **passing q <= 0.05**
  - **564**
  - 
* - what the retired whole-ORF $f_0 \ge 0.5$ would have admitted
  - 884 / 884
  - **100%**
```

Passing starts are led by CTG (152), then GTG (95) -- consistent with the mammalian near-cognate literature
and with the ORF track's own weights ({doc}`model/architecture`). Only **18** ATG extensions survive, against
the 236 that RiboCode labelled: its most-5'-start choice is unsupported roughly 92% of the time.

## FDR: which correction where

```{list-table}
:header-rows: 1
:widths: 24 20 56

* - Column group
  - FDR
  - Why
* - canonical, $\Delta$PSM
  - **global** 1%
  - what the search reports; keeps target-decoy competition intact
* - novel
  - **class-specific** 1%
  - novel targets vs `REV_nuORF|` decoys only
```

A class-specific FDR on the canonical side is confounded: a large novel space cannibalises canonical
**decoys** faster than canonical targets, deflating the estimate until a junk database appears to *gain*
canonical identifications. Raw rank-1 counts are never reported -- they measure spectra *absorbed*, which
scales with database size.

```{admonition} A global cut does not merely inflate the novel class -- it fails in BOTH directions
:class: important
The project carried "a global FDR inflates non-canonical discovery ~10x" as a rule of thumb for a long time.
Measured directly on the four shipped searches, it is only half right, and which way it errs depends on how
much of the search the **canonical** class occupies:

| dataset | digestion | canonical peptides | global / class-specific novel |
|---|---|--:|--:|
| A549 | tryptic, TMT | 65,405 | **37.5x over-report** |
| HBL-1 | nonspecific, HLA-I | 2,429 | 0.6x to 2.3x |
| SU-DHL-4 | nonspecific, HLA-I | 813 | 0.4x to 1.0x |
| DoHH2 | nonspecific, HLA-I | 1,507 | **0.2x under-report** |

On A549 the canonical class is numerous and high-scoring, so it drags the global cut down to hyperscore 19.5
-- inside the novel **decoy** bulk (the class-specific cut is 30.3). In the immunopeptidomes the canonical
class is too sparse to dominate, the global cut lands *stricter*, and real identifications are discarded.

Database size does **not** predict the error: the four `null_atg` arms span 242k-369k novel sequences and give
37.5x / 1.9x / 0.8x / 2.2x. Canonical-class size does. So a global cut is not conservative-by-default here; it
is uncontrolled, and no direction of bias can be assumed. Figure: `figures/S_fdr_rigor/`.
```

**No double counting.** MSFragger reports every protein containing a peptide, so a peptide occurring in any
GENCODE protein is counted as canonical, never novel. This is not rare: 3,489 of 191,358 BMDM rank-1 PSMs map
to both a novel ORF and a GENCODE protein.

```{admonition} MSFragger delimits the protein list with a semicolon, not a comma
:class: warning
Splitting on the wrong character collapses the protein list into a single token. The classification then
survives only because MSFragger happens to list the canonical protein first -- correct by luck, not by
construction. `pgx.report` splits on both.
```

## Is the advantage reproducible, or one lucky dataset?

The single-dataset version of this result was the weakest part of the story. It is now **5 datasets x 2
released models x 2 threshold arms = 20 points, every one above 1.0**, median **93.6x**, range 10.7x to
2,744x.

```{figure} /img/F2b_discovery_forest/F2b_discovery_forest.png
:width: 100%

Discovery-density ratio, model / naive AUG null. Every point favours the model. The Poisson arm leads
$\theta$=1 in every dataset -- the calibration result restated on a fifth independent axis.
```

The metric is a **rate ratio, not a count**, and that choice is load-bearing. On raw count the unselected null
actually *wins* HBL-1 (20 vs 14) and SU-DHL-4 (28 vs 22) -- while carrying 25x to 180x more sequences and
paying a canonical-ID cost the model arms do not. Density asks the question that matters for search-space
selection: per sequence you commit to searching, how much do you find?

The density ratio is a proxy for something now measured directly. On the four HLA-I immunopeptidomes, the
null's *achieved* novel-class FDR is **8.74% [8.04, 9.50] against 1.34% [1.08, 1.67]** for the model arms at
the same nominal 1% -- so the null's extra raw-count discoveries are roughly one in eleven wrong. See
[](immunopeptidomics.md).

Two honest notes for reading the figure:

- **THP-1's 2,744x should not anchor the claim.** Its null found a single novel peptide from 301,855
  sequences, so the denominator is one count away from undefined. THP-1 also has the strictest FDR cut in the
  panel (23.5 vs 17.4-18.2), because a fast-scanning Orbitrap Fusion produced many marginal MS2 spectra that
  generate decoy hits. Quote the median.
- **No error bars, deliberately.** These are single-search point estimates, not replicated measurements. The
  honest uncertainty statement is the peptide count printed beside each point (7 to 16). A Poisson interval on
  those counts would imply a replication structure that does not exist.

  This holds for the deterministic hyperscore searches on this page. It does **not** hold for the rescored
  HLA-I analysis in [](immunopeptidomics.md), where mokapot's stochastic fit gives a genuine replication
  structure over 15 seeds and the pooled decoy counts do support intervals. Do not carry this note across to
  that page.

## The ORF-length floor is set by the ASSAY, not by a default

Before any of the numbers below mean anything, the search database has a minimum ORF length, and it is
**not one value**:

```{list-table}
:header-rows: 1
:widths: 34 12 54

* - assay
  - floor
  - why
* - tryptic whole-cell lysate
  - **30 AA**
  - tryptic peptides come from digesting whole proteins, so sub-30-aa ORFs are not the population the
    assay samples. Same floor the ORF-call side applies (90 nt, stop excluded).
* - MHC / HLA immunopeptidomics
  - **7 AA**
  - HLA-I peptides are 8-11 aa. A 30-aa floor would delete the microproteins the assay exists to find.
```

**Classify from the search parameters, never from the dataset name.** `search_enzyme_name_1 = trypsin`
with `num_enzyme_termini = 2` is tryptic; `nonspecific` with `termini = 0` is MHC.

`pgx.build_dbs --assay {tryptic,mhc}` sets the floor and refuses a conflicting `--min-aa`, so it cannot
drift.

:::{admonition} The floor is load-bearing for HLA-I and free for tryptic
:class: important

Novel peptides at 1% class-specific FDR whose supporting ORFs are **all** sub-30-aa -- the ones a 30-aa
floor deletes -- measured by `proteogenomics/scripts/orf_length_audit.py`:

| arm | mouse macrophage (tryptic) | human (mostly HLA-I) |
|---|--:|--:|
| model (Poisson) | **0 / 295 (0.0%)** | **14 / 52 (26.9%)** |
| null AUG | 19 / 453 (4.2%) | 22 / 72 (30.6%) |
| CPAT | -- | **0 / 39 (0.0%)** |
| CPC2 | -- | **0 / 55 (0.0%)** |

On tryptic data the floor costs the model nothing. On HLA-I it would cost ~30% -- and because CPAT and
CPC2 have **zero** short-only discoveries (their databases are 0-2% short by construction, since
composition scoring selects long ORF-like sequences), the loss is not symmetric. Applying 30 aa to HLA-I
would turn "the model wins all three HLA-I immunopeptidomes" into "wins 2 of 3".

That asymmetry is the point, not an embarrassment: sub-30-aa microproteins are exactly what a
translation model finds and a coding-potential selector rejects.
:::

**Status note (2026-08-14):** the macrophage numbers on this page were computed at 7 AA. They are
tryptic and are scheduled to be rebuilt at 30 AA; because model short-only is 0.0% there, the headline
numbers are expected to move very little. See `docs/STATUS_CURRENT_VS_ARCHIVED.md`.

## Result: a bigger database destroys discovery it contains

BMDM, mouse, `mamba4` union model, tryptic, 18 fractions.

```{list-table}
:header-rows: 1
:widths: 20 10 10 12 12 12 10 14

* - Arm
  - novel PSMs
  - novel pept
  - DB novel seqs
  - GENCODE PSMs
  - $\Delta$PSM
  - ncStart
  - $\Delta$TOTAL pept
* - `gencode`
  - 0
  - 0
  - 0
  - 335,548
  - +0
  - 0
  - +0 (55,156)
* - **`model_poisson`**
  - 121
  - 40
  - **1,533**
  - 335,187
  - **-361**
  - **3**
  - **-21**
* - `null_atg`
  - 177
  - 50
  - 235,686
  - 325,254
  - -10,294
  - 0
  - -1,998
* - `null_nc`
  - 657
  - 177
  - 2,512,388
  - 319,544
  - -16,004
  - 9
  - -4,457
```

:::{admonition} The last column is the one to read first, and it is negative
:class: warning

$\Delta$TOTAL is the change in **total unique peptides** (canonical + novel) against searching
GENCODE alone -- the number an experimentalist actually takes home. Every other column can favour an
arm that nonetheless hands back fewer peptides than not searching a novel database at all.

On this run the model is **-21**, not a gain. Two things follow, and neither is optional:

1. **It flips sign on search parameters.** The frozen-parameter BMDM run (used by
   `figures/P11_bmdm_proteomics` and `F2b`) gives **+31** for the same databases and the same
   spectra. The model's effect on total peptides is within noise of zero, and quoting either number
   without its search configuration is quoting noise as a result.
2. **BMDM is the most favourable of 12 macrophage populations.** Across all 12 (frozen parameters)
   the model is above the GENCODE-only baseline in **5**, median **-8** of ~55,000 (0.01%). The
   naive AUG null is above it in **0 of 12**, median **-1,784**.

The defensible claim is therefore **cost-neutrality**, not a gain: the model adds novel discoveries
at a cost indistinguishable from zero, while the nulls cost ~2,000-4,500 peptides every time. Do not
write "the model finds more peptides".
:::

The model database costs **28x less** canonical signal than the AUG null (0.11% vs 3.07% of baseline) at
**154x** fewer sequences. (Medians across all 12 macrophage populations: **156x** smaller database,
**111x** higher discovery density. The BMDM ratios are representative; its $\Delta$TOTAL is not.) But the decisive number is the peptide overlap: of 40 model peptides, **17 are
found by neither null**, and every one is explained:

```{list-table}
:header-rows: 1
:widths: 66 10 24

* - Cause
  - n
  - 
* - absent from the null's sequence space (non-AUG start)
  - 3
  - unreachable at any FDR
* - present in the null database but below its FDR threshold
  - 13
  - cut rises 20.20 -> 28.10
* - present but not rank-1 in the larger space
  - 1
  - 
```

Database size raises the class-specific hyperscore cut from **20.20** to **28.10**, so 13 real peptides sit
in the null's own FASTA and cannot be reported at 1% FDR. A bloated database does not merely cost canonical
identifications -- it **destroys non-canonical discovery it structurally contains**.

```{admonition} The honest caveat: this is a tradeoff, not a clean win
:class: important
`null_nc` finds 177 novel peptides to the calibrated model's 40, and its 128 unique peptides **replicate
across fractions as well as the model's own** (18.8% single-fraction vs the model's 20.0%; 63.3% seen in
>= 3 fractions vs 52.5%). An early guess that they were chance matches in a 2.5M-sequence space was
**refuted** by that check. So naive near-cognate enumeration is very likely finding real translation the
calibrated caller misses -- at a cost of 4,573 canonical peptides. The `standard` arm ($\theta$=1, 10,454
sequences) exists to probe that middle ground. Report both arms; do not present the calibrated model as
dominant.
```

Replication across fractions (`pgx.replication`) is the one credibility check here that does not depend on
the decoy model at all, which is exactly what is in question for a multi-million-sequence search space.
Read it as a **comparison between arms over identical spectra**, never against an absolute cutoff: genuine
low-abundance peptides do legitimately appear once.

## Running it

```bash
cd proteogenomics/scripts
python -m pgx.run \
  --species mouse --label BMDM \
  --profiles <...>/pred_profiles.npz \
  --salmon <...>/rep{1,2,3}/quant.sf \
  --mzml '<...>/mzml/BMDM_Rep*_F*.mzML' \
  --enzyme tryptic \
  --out <...>/pgx/BMDM \
  --shared-search-root <...>/pgx_shared \
  --submit-searches

# once the MSFragger jobs land
python -m pgx.run ... --report-only
```

`--enzyme tryptic` (fully tryptic, 6-50 aa) or `nonspecific` (`num_enzyme_termini = 0`, 8-14 aa, for HLA
immunopeptidomics). `--db-arms` defaults to `standard,poisson` so calibrated and uncalibrated are always
side by side. `--shared-search-root` keeps the model-independent arms (`gencode`, both nulls) in one place,
so comparing two checkpoints re-searches only the model arms -- a `null_nc` search runs for hours and is
byte-identical across checkpoints.

```{admonition} Never blank a validated enzyme field
:class: warning
`search_enzyme_nocut_1 = P` is classic trypsin, which does **not** cut before proline. Setting it empty
silently switches the search to trypsin/P: measured on BMDM, the GENCODE-only baseline moved
335,548 -> 352,307 PSMs (+5.0%) from that one character, on byte-identical databases. A blank **is** an
override. State every enzyme field explicitly at its validated value.
```
