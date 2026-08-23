# P13 -- per-population novel-PSM Venn: predicted DB vs real BMDM-NT Ribo-seq DB

**What it shows:** the question P12's counts cannot answer. P12 establishes that the two databases
deliver comparable novel-peptide *yield* across 12 macrophage populations. Two arms can each report
27 novel PSMs and either agree completely or not at all. This measures the overlap directly.

## Panel -> data

| element | source |
|---|---|
| all 12 Venns | `proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30/novel_psm_venn.json` |
| produced by | `proteogenomics/scripts/macro_novel_psm_venn.py` (~5 min) |
| upstream PSMs | `pgx_xsubtype_aa30/frozen/<POP>/{model_predicted_attn,model_ribocode_bmdm_nt}/rank1.tsv.gz` |
| arm A (blue) | `model_predicted_attn` -- **attention** model ORFs called PER POPULATION |
| arm B (orange) | `model_ribocode_bmdm_nt` -- REAL Ribo-seq, RiboCode on BMDM untreated, **506 seqs, fixed** |

**Source tree is the 30-AA attention build** (2026-08-15), matching P11 and P12. Rule 5 puts the
tryptic floor at 30 aa and macrophages are trypsin/termini=2. The superseded 7-aa mamba venn is
reproducible with `--source .../pgx_xsubtype/novel_psm_venn.json`; the figure's legend reads the
Ribo-seq database size out of the reports rather than asserting it, because that same database is
699 seqs at 7 aa and 506 at 30 aa.

## THE UNIT IS `(spec_id, peptide)`, NOT `spec_id`

A scan counts as shared only when **both arms assigned it the same peptide**. Scans that pass in both
arms but were assigned DIFFERENT peptides are a third category (`diff-pep`), reported under each
panel and deliberately **excluded from the intersection**: they are rank-1 competition changing the
winner as the search space changes, and folding them in would inflate agreement. A two-circle Venn
cannot draw that category, which is why it is annotated rather than dropped.

**`diff-pep` totals 1 across all 12 populations** -- a single spectrum, in LungResident; every other
population is 0. (This line previously read "is 1 across all 12 populations", which invites reading
it as 1 *per* population.) That is an important QC result in its own right: when both databases pass
a spectrum they essentially always agree on the peptide, so the disagreement this figure measures is
about *which spectra pass*, not about peptide assignment. The value is 1 in both the 7-aa mamba and
30-aa attention builds, so it is a property of the search, not of either database.

## FDR

Each arm carries **its own 1% class-specific cut** (novel targets vs `REV_nuORF|` decoys only). That
is correct -- each database controls its own novel-class error rate -- but it means a spectrum can
sit inside one circle and outside the other purely by threshold. That is a real difference in what a
database delivers at fixed error, not an artifact.

**No FDR code was reimplemented.** `macro_novel_psm_venn.py` imports the pipeline's own
`pgx.report.cut`. Its single-pass loader is asserted equal to `pgx.score_compare.passing_psms` on the
first population (BMDM: 90 PSMs, cut 20.68), and **all 12 x 2 arm counts were cross-checked against
`novel_psms` in the frozen reports with zero mismatches.**

## Result

30-AA attention build. Sorted by Jaccard.

| population | model | Ribo-seq | shared | model-only | Ribo-only | J |
|---|--:|--:|--:|--:|--:|--:|
| LargeIntestinal | 75 | 107 | 46 | 29 | 61 | 0.338 |
| SmallIntestinal | 156 | 107 | 65 | 91 | 42 | 0.328 |
| Peritoneal | 114 | 107 | 50 | 64 | 57 | 0.292 |
| SpleenRecruited | 101 | 102 | 45 | 56 | 57 | 0.285 |
| LungRecruited | 120 | 104 | 49 | 71 | 55 | 0.280 |
| LiverRecruited | 59 | 55 | 23 | 36 | 32 | 0.253 |
| SpleenResident | 89 | 107 | 36 | 53 | 71 | 0.225 |
| RAW264 | 46 | 144 | 33 | 13 | 111 | 0.210 |
| LungResident | 41 | 102 | 24 | 16 | 77 | 0.202 |
| Kupffer | 54 | 36 | 10 | 44 | 26 | 0.125 |
| Microglia | 11 | 94 | 5 | 6 | 89 | 0.050 |
| **BMDM (MATCHED)** | 90 | 61 | 27 | 63 | 34 | 0.218 |

**Transfer totals (n=11, BMDM excluded): model 866 PSMs, Ribo-seq 1,065, shared 386 -- 25.0% of the
1,545-PSM union.** Every population has some overlap (11/11), but Jaccard never exceeds 0.34.

**THE HEADLINE SURVIVED THE FLOOR CHANGE AND THE MODEL SWAP.** The 7-aa mamba build gave 885 / 1,012
/ 367 = 24.0% of a 1,530 union, Jaccard median 0.247, overlap in 11/11. The 30-aa attention build
gives 25.0%, median 0.253, 11/11. Two different ORF-length floors and two different architectures
agree to within one percentage point, so the complementarity below is a property of predicted-vs-
measured databases and not of one particular build.

## What this licenses

**The two databases are largely COMPLEMENTARY, not redundant.** Across the 11 transfer populations
the model finds **479 novel PSMs a real BMDM Ribo-seq database misses**, and that database finds
**678 the model misses**. Neither is a subset of the other in any population.

(Read `a_only` / `b_only` from the JSON, not `A_n - shared`: the latter is 480 / 679 because it
folds the one `diff-pep` spectrum back into the exclusive regions. 479 + 386 + 1 = 866.)

This is the sharpest statement available about what the predicted database adds: not more discoveries
than Ribo-seq (P12 shows it does not), but *different* ones. The practical implication is that the
union outperforms either alone -- which has not been tested here and would need its own search.

## Caveats a caption must carry

- **BMDM is the matched case** (orange frame): the Ribo-seq database was built from BMDM Ribo-seq. It
  is excluded from every transfer statistic.
- **The Ribo-seq database is FIXED across all 12 rows** (506 seqs, BMDM-derived); the model's is
  rebuilt per population and varies 388-1,483 seqs. Overlap differences partly track that.
- **Microglia is the extreme case** (J=0.050, 5 shared): it also has the smallest model database
  (388 seqs) in P12. The two observations are the same fact seen twice.
- Do NOT read circle area as database size -- `matplotlib_venn` scales circles by PSM count.
- This is PSM-level. Peptide-level counts are in the values JSON (`shared_peptides`,
  `a_only_peptides`, `b_only_peptides`) and are smaller, since several PSMs can share a peptide.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
$PY proteogenomics/scripts/macro_novel_psm_venn.py          # ~5 min, writes novel_psm_venn.json
cd figures/P13_macrophage_psm_venn && $PY make_macrophage_psm_venn.py
```

## ORF-length floor: RESOLVED -- this panel now uses 30 AA, matching the ORF-call panels

**Rule 5 (standing): the ORF floor is set by the assay.** Tryptic whole-cell lysate is 30 aa; MHC /
HLA immunopeptidomics stays at 7 aa because HLA-I peptides are 8-11 aa. Macrophages are tryptic
(`fragger.params`: `trypsin`, `termini=2`), so 30 aa applies. Classify from the search parameters,
never from a dataset name.

This panel was built at 7 aa until 2026-08-15 and is now rebuilt at 30 aa, so it sits on the same
floor as the ORF-call track (`build_loader` filters at 90 nt = exactly 30 aa, applied identically to
model and real Ribo-seq calls). The 7-aa version remains reproducible via `--source` and is retained
only for the comparison recorded above.

**The earlier audit predicted this, and it held.** Before the rebuild, novel peptides at 1%
class-specific FDR whose supporting ORFs were ALL under 30 aa were: model 0 of 295 (0.0%), real
Ribo-seq 8 of 306 (2.6%), `null_atg` 19 of 453 (4.2%) -- so the floor change was expected to cost
the model nothing and the null the most. Measured outcome: shared fraction 24.0% -> 25.0%, Jaccard
median 0.247 -> 0.253, overlap still 11/11.

Measured by `proteogenomics/scripts/orf_length_audit.py`; declared in `docs/PIPELINE_POLICY.md`.
Re-run it if the databases are rebuilt -- do not assume the 0% holds.
