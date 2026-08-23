# Miniplan: RiboCode-called proteogenomics search pipeline (`pgx`)

Scope-limited plan, NOT a replacement for the `prepare/` CLI plan (that one resumes after this).
Status key: PLANNED / IN-PROGRESS / DONE. Progress log at the bottom.

## Why this exists

The current DB build (`enumerate_score_orfs.py` + `build_a549_dbs.py`) selects novel ORFs with a
single scalar cut on `pred_frame0`, bypassing every ORF-calling component this project built and
validated. Two measured consequences:

1. **`f0` scores reading frame, not initiation.** `f0 = pm[0::3].sum()/pm.sum()` is a ratio over
   the ORF interval, so appending a zero-signal upstream region changes neither numerator nor
   denominator: `f0` is *invariant* to N-terminal extension, not merely diluted. Measured on BMDM
   near-cognate candidates, over 12,719 stop-codon groups containing an annotated CDS start plus
   upstream in-frame starts: **96.2%** of 51,646 upstream non-AUG ORFs land within 0.05 of the true
   CDS's `f0`; **median |delta f0| = 0.0026**. No threshold on `f0` can separate them at any value.
2. **No Poisson, no periodicity test, no significance.** The only gate is `f0 >= thresh`
   (`build_a549_dbs.py:103`). `pred_density` is written to the TSV and never read. Poisson injection
   is documented as the winning calibration lever (+0.10 novel precision, Check 4) and was never in
   this path.

Enabling non-AUG enumeration exposed this by multiplying the candidate pool ~10x (625,698 ->
5,956,246), which `f0` passed at the same ~7.4% rate, ballooning the model DB 13,835 -> 142,248 and
destroying the size advantage that made it interesting.

**Fix:** call ORFs with the validated caller (Poisson-calibrated RiboCode on the predicted signal)
and take the DB from significant calls.

## Design decisions (settled before coding)

| Decision | Choice | Rationale |
|---|---|---|
| Calibration anchor | **annotated CDS of expressed transcripts** | self-contained: the model predicts from sequence + RNA-seq, so the pipeline must not require real Ribo-seq calls to calibrate |
| theta selection | **smallest theta with CDS recall >= target** (default 0.90) | tightening monotonically improves novel precision (Check 3); take the most stringent point still meeting the recall floor |
| DB arm | **Poisson arm** by default (`--db-arm`) | the recommended de novo recipe; standard arm still called + reported, so the two-arm rule is met at the calling level without doubling MS searches |
| Significance | RiboCode `adjusted_pval` (BH), `--orf-qvalue 0.05` | RiboCode already runs `pval_adj="fdr_bh"` |
| Universe | TPM >= 1, `{protein_coding, lncRNA}`, len <= 10000, no chrM | user spec; identical for model and BOTH nulls so the arms are comparable |
| Biotype field | `transcript_type` (default), `gene_type` available | OPEN QUESTION, see below |
| null_nc composition | ATG + the 9 near-cognates (superset of null_atg) | makes `null_atg` a strict subset of `null_nc`, a clean nesting for interpretation |
| GENCODE double-count | drop at DB build (seq == GENCODE protein) AND reclassify at report (peptide is a substring of any GENCODE protein) | a peptide is a fragment, so set membership alone is insufficient |

### OPEN QUESTION: which biotype field the pc + lncRNA filter applies to

`build_line_universe.py` (and therefore every existing universe in this project, including the
mouse macrophage ones) filters on **`transcript_type`**, keeping only canonical coding isoforms.
Filtering on **`gene_type`** instead keeps every isoform of a protein-coding or lncRNA gene, which
additionally admits `nonsense_mediated_decay`, `retained_intron`, `TEC` and `non_stop_decay`
transcripts. Measured on A549: **45,166 vs 41,139** transcripts, the 4,027 extra being 3,747 NMD,
106 retained_intron, 52 TEC, 21 non_stop_decay.

The distinction matters here more than usual, because NMD and retained-intron isoforms are exactly
where non-canonical ORFs are enriched, so `gene_type` is the more inclusive choice for discovery
while `transcript_type` is the conservative one that keeps comparability with prior results.
The user's phrasing was "pc-genes", which reads as `gene_type`; the default is currently
`transcript_type` for continuity, exposed as `--biotype-field`. **Flag for decision.**

### The extension problem (why extensions need their own code path)

RiboCode's alternative-start mode is **fallback-only per common-stop**. Verified in
`orf_finder.py:orf_find`: `alt_flag` is set to 0 as soon as any in-frame ATG exists between the
previous stop and the current one, and the alt-start list is then skipped entirely. An annotated CDS
always has an ATG, so **RiboCode can never emit a near-cognate N-terminal extension of an annotated
CDS**, at any parameter setting.

However `classfy_orf` returns `"annotated"` when `orfiv.end == cdsiv.end`, and `start_check` with
`only_longest_orf=True` takes `start_list[0]` (the most 5' in-frame start). So RiboCode *does*
already emit **ATG** extensions, labelled `annotated`, detectable as `ORF_tstart < annotated_tstart`.

Therefore extensions come from two sources:
- **ATG extensions**: harvested from RiboCode `annotated` rows (free, already significance-tested).
- **Near-cognate extensions**: purpose-built scan, tested with RiboCode's own `test_frame`
  (Wilcoxon f0>f1, f0>f2, Stouffer-combined) applied to the **extension region alone** (proposed
  start -> first in-frame ATG), then BH-adjusted. Scoring the extension region alone is exactly what
  defeats the `f0` inheritance above: it asks whether *this start* is used, not whether the ORF sits
  in a translated frame.

Both land in a distinct accession namespace (`nuORF|nterm_ext|...`) so they are separable in every
table.

## Layout

New package `proteogenomics/scripts/pgx/`, argparse CLI per step (so steps fan out as sbatch arrays
and resume independently), plus one orchestrator. No hardcoded paths: every reference resolves
CLI arg > env var > species registry.

```
pgx/
  refs.py         species registry + path resolution (human/mouse)
  universe.py     salmon quant.sf -> expressed universe tx list + FASTA
  calibrate.py    Poisson/theta CDS-anchored calibration (implements Check 6)
  call_orfs.py    two-arm RiboCode calling on the predicted signal + q-filter
  extensions.py   CDS N-terminal extension scan (ATG harvest + near-cognate test)
  build_dbs.py    db_gencode / db_model / db_null_atg / db_null_nc
  search.py       MSFragger runner, --enzyme tryptic|nonspecific, hash-cached
  report.py       FDR + the mandated output table
  run.py          end-to-end orchestrator
```

## Steps

### S1. `refs.py` + package skeleton    Status: PLANNED
Species registry for human (v49) and mouse (vM38): GENCODE proteome, `tx2cds`, `tx2biotype`,
pc+lncRNA transcript FASTAs, RiboCode annotation dir. Resolution order CLI > env
(`PGX_<KEY>`) > registry. Verified present on disk:
`gencode.{v49,vM38}.pc_translations.fa.gz`, `data/{,mouse_}tx2cds.tsv`,
`data/tx2biotype{,_mouse}.tsv`, mouse annot `data/mouse_ribocode_annot/`, human annot
`biotype_probe/expression_context_human/data/ribocode_annot/`.
Exit check: `python -m pgx.refs --species mouse --print` lists existing files only.

### S2. `universe.py`    Status: PLANNED
Wrap the validated `build_line_universe.py` rule (TPM >= `--min-tpm`, pc+lncRNA, <= `--max-len`, no
chrM) behind the package, multi-quant mean-TPM merge included (from `build_macro_universes.py`).
ONE universe object feeds the model arm and both nulls.
Exit check: reproduces an existing macrophage universe tx list byte-for-byte.

### S3. `calibrate.py`    Status: PLANNED
Implements Check 6 (`--cds-recall`), previously PLANNED and never built. Sweep theta over
`ribocode_dropin.py --variant pred_preddepth --pred_poisson --pred_scale <theta>`; at each theta
compute recall of **annotated CDS on expressed transcripts**; pick the smallest theta with
recall >= target. Emit `theta_curve.tsv` + `theta.json`.
Exit check: monotone curve; chosen theta reproduces the Check 3 ballpark on Hepatocytes.

### S4. `call_orfs.py`    Status: PLANNED
Both arms through `ribocode_dropin.py`: standard (`--pred_scale 1.0`, deterministic) and Poisson
(`--pred_scale theta* --pred_poisson`). Filter on `adjusted_pval <= --orf-qvalue`. Map RiboCode
`ORF_type` -> `{canonical, uORF, dORF, lncRNA_orf, internal}` reusing `build_ribocode_db.py`'s
`TYPE_MAP`/`SKIP_TYPES`. Emit per-arm call tables + a comparison summary.
Exit check: standard arm reproduces an existing drop-in run's call count.

### S5. `extensions.py`    Status: PLANNED
(a) Harvest ATG extensions from arm output: `ORF_type == annotated` and
`ORF_tstart < annotated_tstart`.
(b) Near-cognate scan: for each expressed tx with an annotated CDS, enumerate upstream in-frame
near-cognate starts before the first upstream in-frame stop; extract the **extension region only**
(start -> annotated CDS start); require a minimum in-frame P-site count; apply RiboCode
`extract_frame` + `test_frame`; BH-adjust across all candidates; keep q <= `--ext-qvalue`.
Emit `extensions.tsv` with `start_codon`, `ext_len_aa`, `pval`, `qval`, plus the extended AA seq.
Exit check: the fraction of near-cognate extensions surviving is far below the 88.7% that bare `f0`
passed, and the 96.2% inheritance cohort collapses.

### S6. `build_dbs.py`    Status: PLANNED
Four DBs, shared GENCODE base + inline `REV_` decoys (`decoy_prefix = REV_`):
- `db_gencode`  GENCODE proteome only (baseline detection)
- `db_model`    GENCODE + significant RiboCode calls + extension entries
- `db_null_atg` GENCODE + all ATG ORFs on the expressed universe
- `db_null_nc`  GENCODE + all ATG + near-cognate ORFs on the same universe
Dedup by AA sequence; drop any novel sequence identical to a GENCODE protein. Write a class map
(accession -> class, start codon, source arm) per DB for class-specific FDR.
Exit check: `db_null_atg` novel set is a strict subset of `db_null_nc`; no novel sequence equals a
GENCODE protein in any DB.

### S7. `search.py`    Status: PLANNED
MSFragger runner. `--enzyme tryptic` -> `search_enzyme_name_1 = trypsin`, `num_enzyme_termini = 2`,
`digest 6-50`; `--enzyme nonspecific` -> `nonspecific`, `num_enzyme_termini = 0`, `digest 8-14`
(the two existing param files, templated rather than copied). Content-hash each DB FASTA; skip a
search whose `<out>/.pgx_hash` already matches. The `db_gencode` baseline is searched ONCE per
(mzML set, enzyme) and reused by every arm's reporting.
Exit check: rerunning the orchestrator launches zero redundant searches.

### S8. `report.py`    Status: PLANNED
Peptide classification: build the GENCODE AA space once; any peptide that is a substring of any
GENCODE protein counts as canonical, never novel (the no-double-count rule). FDR per the standing
rule: **global** 1% for canonical/dPC, **class-specific** 1% (novel targets vs `REV_nuORF|` decoys)
for novel. Columns, per arm:

| column | definition |
|---|---|
| `novel PSMs` | PSMs assigned to novel targets at class-specific 1% FDR |
| `novel seqs w/ PSM` | distinct non-GENCODE DB sequences with >= 1 passing PSM |
| `GENCODE matches` | canonical peptides at global 1% FDR |
| `dPC` | change in GENCODE PSMs vs the `db_gencode` baseline (gained/lost to FDR) |
| `DB seqs` | novel target sequences in that DB (the denominator) |
| `ncStart` | passing model peptides absent from the null AA space (substring scan) |

Exit check: `db_gencode` row has `dPC == 0` by construction; totals reconcile across arms.

### S9. `run.py`    Status: PLANNED
End-to-end orchestrator: `--species`, `--profiles`, `--salmon`, `--mzml`, `--enzyme`, `--out`,
`--cds-recall`, `--orf-qvalue`, `--ext-qvalue`, `--db-arm`, plus reference overrides. Resumable;
prints the sbatch commands for the heavy steps rather than assuming a scheduler.
Exit check: a dry run on BMDM emits a complete, path-correct step list.

### S10. Validation run    Status: PLANNED
BMDM (mouse), tryptic. Compare against the existing `f0`-threshold DBs: DB sizes, novel yield,
dPC in both directions, ncStart. Write results to `results.md` and record the method in `methods.md`.

## Results (BMDM, mouse, mamba4 union model, tryptic)

### Calibration (S3) -- theta* = 0.05

The absolute anchor had to be corrected after measurement. Recall of ALL annotated CDS in the
universe plateaus at **0.533** (10,236 of 19,190), so a 0.90 absolute target is unreachable by
construction. The dial is now anchored on the **achievable** CDS set (`--recall-mode relative`,
hits / hits at the saturating theta; plateau confirmed, top two thetas within 2%). theta* = 0.05
reproduces the operating point Check 5 independently validated as transferable across datasets.

**Why the plateau sits at 0.533 -- annotation redundancy, NOT detection failure.** An earlier
version of this section claimed "roughly half of expressed CDS-bearing transcripts are not
detectably translated"; that is wrong. The universe carries **1.90** CDS-bearing isoforms per gene,
Ribo-seq signal concentrates on essentially one of them, and the caller finds it. At the ceiling
those 10,236 transcripts span **9,777 of 10,080 genes**: transcript-level recall 53.3% is
**gene-level recall 97.0%**, at 1.05 called isoforms per gene. The ceiling therefore measures how
many isoforms the annotation lists, not what the sample translates -- which is also why adding
2,071 NMD transcripts under `gene_type` moved the denominator but barely moved the ceiling (+72):
those are extra isoforms of genes already covered.

| theta | CDS recall (abs) | CDS recall (rel) | calls | uORF | dORF | lncRNA_orf |
|---|--:|--:|--:|--:|--:|--:|
| 0.02 | 0.4483 | 0.8404 | 8,999 | 178 | 28 | 183 |
| **0.05** | **0.5069** | **0.9504** | 10,786 | 509 | 58 | 455 |
| 0.10 | 0.5261 | 0.9863 | 11,996 | 918 | 107 | 822 |
| 0.20 | 0.5319 | 0.9973 | 13,417 | 1,614 | 185 | 1,342 |
| 1.00 | 0.5334 | 1.0000 | 17,436 | 3,445 | 822 | 2,793 |

At theta* the dial keeps 95% of achievable CDS while cutting uORF calls 6.8x and dORF calls 14x.

### Two-arm calling (S4)

| arm | theta | calls | novel | uORF | dORF | lncRNA_orf | annotated | of which nterm_ext |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| standard | 1.00 | 16,960 | 6,623 | 3,496 | 486 | 2,641 | 10,239 | 256 |
| poisson | 0.05 | 10,786 | 1,022 | 509 | 58 | 455 | 9,730 | 236 |

The poisson arm at theta* was reused from the calibration sweep, not recomputed.

**Significance gate.** RiboCode filters per-ORF on `--pval` BEFORE BH, so on the surviving set BH
moves p by at most 4e-4 and `--orf-qvalue 0.05` is a no-op. `--pval` is the real lever; the
q-filter only tightens.

### Extension scan (S5) -- the diagnosis confirmed and fixed

| stage | n |
|---|--:|
| translated CDS targets | 9,730 |
| candidate upstream in-frame starts | 30,243 |
| dropped, too short | 5,396 |
| dropped, < 5 nonzero frame-0 codons in the extension region | 23,963 |
| testable | 884 |
| **passing q <= 0.05** | **564** |

Passing by start codon: CTG 152, GTG 95, AAG 72, AGG 67, TTG 46, ACG 37, ATC 34, ATT 27, ATG 18,
ATA 16. CTG leading matches both the mammalian near-cognate literature and the ORF track's own
weighting.

**The contrast that justifies the step:** of the 884 candidates that reached the test, the retired
whole-ORF `f0 >= 0.5` criterion would have admitted **884 of 884 (100%)**, and by inheritance it
would have admitted essentially all 30,243. The extension-region test admits 564, a **54x**
reduction, on evidence rather than on a threshold.

Note also that RiboCode itself reported 236 `nterm_ext` among its `annotated` calls, but only 18
ATG extensions survive the region-specific test: `only_longest_orf=True` picks the most 5' start
without ever testing it, so that label is unsupported roughly 92% of the time.

### Databases (S6), after the classifier fix

| database | targets | novel seqs | vs model |
|---|--:|--:|--:|
| `db_gencode` | 57,226 | 0 | baseline |
| `db_model_poisson` | 58,759 | **1,533** (972 calls + 561 extensions) | 1x |
| `db_null_atg` | 292,912 | 235,686 | 154x |
| `db_null_nc` | 2,569,614 | 2,512,388 | **1,639x** |

Model DB by class: nterm_ext 561, uORF 478, lncRNA_orf 442, dORF 52. By start codon: ATG 987 plus
546 non-AUG. For comparison the retired f0 build gave 13,835 novel at ATG-only and 142,248 once
near-cognates were enumerated.

Invariants verified on the built FASTAs:

| check | result |
|---|---|
| `null_atg` is a strict subset of `null_nc` | PASS |
| no novel sequence identical to a GENCODE protein (any DB) | PASS (0 / 0 / 0) |
| model novel seqs absent from `null_nc` | **0** (was 117 before the classifier fix) |
| model novel seqs absent from `null_atg` | 546, of which non-AUG **546 / 546 = 100%** |

The last row is the point of the whole design: every model sequence an AUG-only pipeline cannot
reach is unreachable *because of its start codon*, not because of any bookkeeping difference
between the arms. That is what makes `ncStart` interpretable.

### Correction: the null enumerator and RiboCode must classify ORFs by the SAME rule

The first build left the nulls missing an entire ORF class the model arm contains. RiboCode's
`classfy_orf` treats an out-of-frame ORF that overlaps the CDS as `Overlap_uORF` (opens 5' of the
CDS start) or `Overlap_dORF` (ends 3' of the CDS end), both of which map to novel classes. The
first `seqtools.cds_relationship` collapsed every out-of-frame overlap to `internal` and dropped
it, so those ORFs entered the model DB but could never enter a null.

Detected by an invariant check, not by inspection: **117 of 1,533 model sequences were absent from
the 2.3M-sequence near-cognate null**, and all 117 were **ATG**, from RiboCode calls, classed uORF
(107) and dORF (10). A near-cognate null that a model sequence can escape on an AUG start is not a
null. `cds_relationship` now splits out-of-frame overlaps three ways to mirror `classfy_orf`
exactly, and only a fully-contained out-of-frame ORF remains `internal`.

This also forced a sharper definition of `ncStart`. Absence from a null has two possible causes:
the start codon (what the column is meant to measure) and classifier disagreement (an artifact).
The table now reports both -- `absent vs <null>` for raw absence, and `ncStart vs <null>` for the
subset whose supporting DB entry opens at a non-AUG codon, which is literally "absent by virtue of
initiating at a non-canonical start" and is immune to any future classifier drift.

### Peptide-to-GENCODE assignment (S8)

MSFragger delimits its protein list with `;`, not `,`. The pre-existing reporting scripts split on
`,`, which collapses the list to a single token; the classification then survives only because
MSFragger happens to emit the canonical protein first. Measured: 3,489 of 191,358 rank-1 PSMs in
the BMDM model search map to BOTH a novel ORF and a GENCODE protein, and all 3,489 are ordered
canonical-first, so the old code was correct by luck rather than construction. `pgx.report` splits
on both characters, so the no-double-count rule holds regardless of ordering.

### Correction: never override a validated enzyme setting with an empty value

The first `pgx.search` ENZYME table set `search_enzyme_nocut_1 = ""` for the tryptic mode. The
validated template has `P`, i.e. classic trypsin, which does NOT cut before proline; the empty
value silently switches the search to trypsin/P and changes the peptide space.

Caught by comparing baselines rather than by reading the code: the GENCODE-only baseline came out
at **352,307** PSMs against the legacy final-recipe arm's **335,548**, a +5.0% gap on databases that
were verified byte-identical (57,226 target sequences, same sequence set). A `diff` of the rendered
parameters against the template isolated it to that one character. All four first-round searches
were discarded and rerun. The tryptic block now renders identical to
`fragger_macro_lfq.params` apart from `database_name`, which is asserted in the progress log.

Lesson for the enzyme table generally: state every enzyme field explicitly and equal to the
validated value, rather than leaving fields blank and assuming the template shows through -- a
blank IS an override.

### Correction: the search cache key must not contain the database PATH

`--shared-search-root` silently did nothing. The cache digest was
`sha256(database) + sha256(rendered params)`, and the rendered params contain `database_name`, the
FASTA's absolute path -- which is per-run. So the attn run re-searched `gencode`, `null_atg` and
`null_nc` even though `md5sum` proved all three FASTAs byte-identical to the mamba4 ones
(`c0d66cb6...`, `7c0f8010...`, `d87a6cf8...`). Worse, each run would then rewrite `.pgx_hash` with
its own path, so the two checkpoints would have ping-ponged, re-searching a 2.5M-target database
every time.

`database_name` is now excluded from the settings digest; the database's CONTENT is already covered
by `sha256(database)`, and where the file lives is not part of what the search is. Verified: the
three shared arms now hash identically from either run's copy of the FASTA.

### Consensus over Poisson draws (the randomness the calibrated arm carries)

The calibrated arm samples `Poisson(rate)`, so one run's call list is one sample from a
distribution. Five seeds at theta* = 0.05 on BMDM:

| quantity | value |
|---|---|
| total calls | 10,823 +- 23 (**CV 0.21%**) |
| pairwise Jaccard (10 pairs) | **0.822** |
| called in all 5 draws | 8,869 = **66.9%** of the 13,250 union |
| called in exactly one draw | 1,673 = **12.6%** |

Counts reproduce; membership does not. The instability is concentrated where it matters least for
calibration and most for discovery: `annotated` varies 0.16% (which is why theta* itself is stable),
while dORF varies ~17%, lncRNA_orf ~7%, uORF ~4%.

`pgx.consensus` keeps ORFs recurring in >= K of N draws. On BMDM at >= 3 of 5 it removes **189 of
1,022 novel calls (18.5%)** as draw-specific -- calls that would otherwise have entered a search
database as findings. Use it whenever the call list is an end product.

### Baseline for comparison: the LEGACY f0-threshold databases (same BMDM spectra)

Produced by `pgx.report` against the existing `search/BMDM` output, so the FDR treatment is
identical to the new arms and only the database construction differs.

| arm | novel PSMs | novel pept | novel seqs w/PSM | GENCODE PSMs | dPSM | dPept |
|---|--:|--:|--:|--:|--:|--:|
| canonical | 0 | 0 | 0 | 335,548 | +0 | +0 |
| model (f0 >= 0.5) | 205 | 46 | 104 | 333,022 | -2,526 | -480 |
| null (all ATG candidates) | 140 | 41 | 95 | 326,434 | -9,114 | -1,814 |

The legacy model DB already beat its null on both axes (more novel, smaller canonical cost). The
question S10 answers is whether calling ORFs properly improves on that, at 9x fewer novel sequences
(1,533 vs 13,835).

### S10 RESULT (BMDM, corrected enzyme): the model DB recovers what the null destroys

| arm | novel PSMs | novel pept | novel seqs w/PSM | DB novel seqs | GENCODE PSMs | dPSM | dPept | ncStart |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| gencode | 0 | 0 | 0 | 0 | 335,548 | +0 | +0 | 0 |
| **model_poisson** | 121 | 40 | 20 | **1,533** | 335,187 | **-361** | -61 | **3** |
| null_atg | 177 | 50 | 104 | 235,686 | 325,254 | **-10,294** | -2,048 | 0 |

The GENCODE baseline reproduces the legacy final-recipe arm EXACTLY (335,548 both), confirming the
enzyme fix and putting both builds on identical footing.

154x smaller DB, 28x lower canonical cost, 123x higher discovery density (26.1 vs 0.21 novel
peptides per 1,000 DB sequences). Using the null instead buys +10 novel peptides for -1,987
canonical peptides.

> **These are BMDM numbers and BMDM is the most favourable of the 12 macrophage populations.** Do
> not generalise them. On TOTAL unique peptides BMDM gains +31, the best of the 12; the median
> across all populations is **-8**, and the model is above the GENCODE-only baseline in only **5 of
> 12** (the naive AUG null in 0 of 12, median -1,784). Medians across the 12: DB 156x smaller,
> discovery density 111x higher. The defensible cross-population claim is COST-NEUTRALITY, not a
> gain. See `figures/P11_bmdm_proteomics/FIGURE_DATA_INPUTS.md` and results.md 2026-08-14.

**Peptide overlap: 23 shared, 17 model-only, 27 null-only.** Every model-only peptide is explained:

| cause | n |
|---|--:|
| absent from the null's sequence space (non-AUG start, unreachable at any FDR) | 3 |
| present in the null DB but below its FDR threshold (cut rises 20.20 -> 28.10, +7.91) | 13 |
| present but not seen as rank-1 in the larger space | 1 |

A larger database does not merely cost canonical identifications, it **destroys non-canonical
discovery it structurally contains**: 13 real peptides sit in the null's own FASTA and cannot be
reported at 1% FDR because the null's own size raised the bar.

## Progress log

- 2026-08-01 -- plan written. Diagnosis measured (96.2% f0 inheritance, median |delta f0| 0.0026);
  RiboCode alt-start fallback-only semantics confirmed in `orf_finder.py:orf_find`; ATG extensions
  confirmed available via `classfy_orf` + `only_longest_orf`; all species reference files verified
  present for human and mouse.
- 2026-08-01 -- S1 DONE (`refs.py`, both species fully resolved). S2 DONE (`universe.py`,
  reproduces the existing BMDM universe BYTE-IDENTICAL, 21,393 tx). S3 DONE (`calibrate.py`,
  implements Check 6; absolute anchor corrected to relative after measurement). S4 DONE
  (`call_orfs.py`, both arms, poisson arm reused from the sweep). S5 DONE (`extensions.py`).
  S6 DONE (`build_dbs.py`, all four DBs). S7 DONE (`search.py` + `search.sbatch`, hash-cached).
  S8 DONE (`report.py`). S9 DONE (`run.py` + `pipeline.sbatch`). S10 IN-PROGRESS.
- 2026-08-01 -- corrections found by running, not by inspection: (a) the absolute CDS-recall anchor
  is unreachable, replaced by the relative anchor; (b) `num_slices` is not an MSFragger 4.2
  parameter, removed rather than left as a warning-generating no-op; (c) the null enumerator
  disagreed with RiboCode on out-of-frame CDS overlaps, leaving the nulls missing a class the model
  had (117 ATG sequences escaped a 2.3M null) -- classifiers aligned, verified 0 escapes;
  (d) `ncStart` split into `absent vs <null>` + `ncStart vs <null>` so start-codon absence is never
  conflated with classifier drift. `report.py` validated end-to-end on the legacy search output.
  Nulls rebuilt and re-searched (jobs 35986865/66); `gencode` + `model_poisson` unaffected
  (35986855/56). Added `--shared-search-root` so model-independent arms are searched once across
  checkpoints, and `--biotype-field` for the open transcript_type vs gene_type question.
