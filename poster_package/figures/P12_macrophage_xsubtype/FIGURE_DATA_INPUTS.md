# P12 -- macrophage sweep: PREDICTED vs REAL Ribo-seq ORF database, 12 populations

**What it shows:** the harder question P11 does not ask. P11 compares the model against naive
enumeration; this compares it against a database built from **actual ribosome profiling**. The data
existed as `pgx_xsubtype/crosssubtype_table.md` with no figure until 2026-08-14.

## Panel -> data

| element | source |
|---|---|
| all three panels | `proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30/reports/report_<POP>.json` (12) |
| arm `model_predicted_attn` | **attention** model ORFs called PER POPULATION from that population's own data |
| arm `model_ribocode_bmdm_nt` | REAL Ribo-seq, RiboCode on BMDM untreated -- **506 seqs, fixed** |
| arm `model_ribocode_bmdm` | REAL Ribo-seq, RiboCode on BMDM NT + LPS -- **1,265 seqs, fixed** |
| baseline | `gencode` arm of the same report |

**Source tree is the 30-AA attention build** (2026-08-15), matching P11 and P13. Rule 5 puts the
tryptic floor at 30 aa and macrophages are trypsin/termini=2. `--source` / `--arm` reproduce the
superseded 7-aa mamba panel.

**Database sizes are measured at plot time, not written into the labels.** The 7-aa panel hardcoded
"699 seqs" and "2,173 seqs"; the SAME two databases are 506 and 1,265 at the 30-aa floor, so a fixed
label would have silently misreported them the moment the source tree moved. The suptitle likewise
derives its comparison from the numbers and announces a reversal if one ever occurs.

Frozen search parameters, matching P11 and F2b so the panels sit on one footing. Novel columns at
1% class-specific FDR; canonical at 1% global FDR.

**Depth is balanced**: every one of the 12 populations has 18 mzML fractions, so cross-population
comparison is not confounded by acquisition depth.

## THE ASYMMETRY THAT DEFINES THE COMPARISON

The two Ribo-seq databases are **BMDM-derived and SHARED by every row** -- one fixed database applied
to all 12 populations. The model's database is **rebuilt per population**. So:

- **BMDM is the MATCHED case** for the Ribo-seq arms (their home ground), drawn after a dotted rule
  and excluded from every summary statistic.
- **The other 11 are a TRANSFER test**: can a fixed BMDM Ribo-seq database serve a different
  macrophage population?

Including BMDM in the medians would flatter the Ribo-seq arms on their own data. All numbers below
are transfer-only, n=11.

## Result (transfer populations, n=11)

30-AA attention build.

| arm | novel peptides (median) | delta TOTAL peptides (median) | above zero | discovery density (median) | DB size |
|---|--:|--:|--:|--:|--:|
| model (predicted, per population) | 24 | **-3** | **5/11** | **21.3** | median 1,091 (388-1,483) |
| real Ribo-seq, BMDM NT | 26 | **+18** | **8/11** | **51.4** | **506** |
| real Ribo-seq, BMDM NT+LPS | 34 | +12 | 6/11 | 26.9 | 1,265 |

**A fixed 506-sequence database derived from BMDM Ribo-seq matches or beats the per-population
predicted database on every axis**: comparable novel yield, better on total unique peptides, and
2.4x the discovery density on less than half the database size.

BMDM itself (matched, NOT in the medians): model 34 novel / **+11** total; real NT 24 / +10; real
NT+LPS 37 / +17.

### Correcting the floor did not rescue the model here

Worth stating because the 30-AA rebuild *helped* the model in P11. The superseded 7-aa mamba panel
read 26 vs 27 novel, -8 vs +4 total (4/11 vs 7/11), 13.1 vs 38.6 density (2.9x). Moving to 30 aa and
the attention model **narrowed the density gap** (2.9x -> 2.4x) but **widened the total-peptide gap**
(-8/+4 -> -3/+18). The conclusion is unchanged and this panel remains the honest counterweight to
P11.

## What this DOES and does NOT license

**Do NOT write "the predicted database beats Ribo-seq."** It does not, on the transfer populations.

**The defensible claim is that the model does not REQUIRE Ribo-seq.** The predicted arm needs only
RNA-seq plus sequence; every Ribo-seq arm here needed a ribosome profiling experiment in some
macrophage. That is the distinction a reader should leave with, and it is a claim about cost and
applicability, not about accuracy.

## Caveats a caption must carry

- **The delta-total differences are not a ranking.** They are the same size as the frozen-vs-unfrozen
  search-parameter swing measured on BMDM (+31 vs -21 for identical databases and spectra, results.md
  2026-08-14). Panel B says all three arms straddle zero; it does not say which is best.
- **The real-Ribo-seq arms are ATG-only by construction.** RiboCode reports N-terminal extensions
  without testing whether the upstream start is used, so their ncStart is structurally zero. That is
  a property of the caller, not evidence about non-AUG initiation, and it is the one axis where the
  predicted arm has a structural advantage (3 ncStart peptides on BMDM).
- **Discovery density is the robust axis here**, being a ratio insensitive to the search-parameter
  swing that moves the peptide totals.
- Null arms (`null_atg`, `null_nc`) are in the same reports but are NOT plotted here -- they are
  P11's comparison. Their database sizes are per-population and not cross-comparable.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/P12_macrophage_xsubtype && $PY make_macrophage_xsubtype.py
```

## ORF-length floor: RESOLVED -- this panel now uses 30 AA, matching the ORF-call panels

**Rule 5 (standing): the ORF floor is set by the assay.** Tryptic whole-cell lysate is 30 aa; MHC /
HLA immunopeptidomics stays at 7 aa because HLA-I peptides are 8-11 aa. Macrophages are tryptic
(`fragger.params`: `trypsin`, `termini=2`), so 30 aa applies. Classify from the search parameters,
never from a dataset name.

This panel was built at 7 aa until 2026-08-15 and is now rebuilt at 30 aa, so it sits on the same
floor as the ORF-call track (`build_loader` filters at 90 nt = exactly 30 aa, applied identically to
model and real Ribo-seq calls). The 7-aa version is reproducible via `--source` / `--arm` and is
retained only for the comparison recorded above.

The floor moved every database, not just the model's: `null_atg` fell from a median 245,312 to
84,226 sequences (the ~65% that were sub-30-aa), the two Ribo-seq databases from 699/2,173 to
506/1,265, and the model's own from a median 1,578 (mamba) to 1,091 (attention).

Measured by `proteogenomics/scripts/orf_length_audit.py`; declared in `docs/PIPELINE_POLICY.md`.
Re-run it if the databases are rebuilt -- do not assume the earlier 0% holds.
