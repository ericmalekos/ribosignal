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
```

`db_null_atg` is a strict subset of `db_null_nc`, so the increment between them is exactly the cost of naive
non-AUG enumeration. Numbers are BMDM, mouse, `mamba4` union model.

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
canonical identifications. A global FDR on the novel side inflates non-canonical discovery ~10-13x. Raw rank-1
counts are never reported -- they measure spectra *absorbed*, which scales with database size.

**No double counting.** MSFragger reports every protein containing a peptide, so a peptide occurring in any
GENCODE protein is counted as canonical, never novel. This is not rare: 3,489 of 191,358 BMDM rank-1 PSMs map
to both a novel ORF and a GENCODE protein.

```{admonition} MSFragger delimits `proteins` with `;`, not `,`
:class: warning
Splitting on the wrong character collapses the protein list into a single token. The classification then
survives only because MSFragger happens to list the canonical protein first -- correct by luck, not by
construction. `pgx.report` splits on both.
```

## Result: a bigger database destroys discovery it contains

BMDM, mouse, `mamba4` union model, tryptic, 18 fractions.

```{list-table}
:header-rows: 1
:widths: 22 12 12 14 14 14 12

* - Arm
  - novel PSMs
  - novel pept
  - DB novel seqs
  - GENCODE PSMs
  - $\Delta$PSM
  - ncStart
* - `gencode`
  - 0
  - 0
  - 0
  - 335,548
  - +0
  - 0
* - **`model_poisson`**
  - 121
  - 40
  - **1,533**
  - 335,187
  - **-361**
  - **3**
* - `null_atg`
  - 177
  - 50
  - 235,686
  - 325,254
  - -10,294
  - 0
* - `null_nc`
  - 657
  - 177
  - 2,512,388
  - 319,544
  - -16,004
  - 9
```

The model database costs **28x less** canonical signal than the AUG null (0.11% vs 3.07% of baseline) at
**154x** fewer sequences. But the decisive number is the peptide overlap: of 40 model peptides, **17 are
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
