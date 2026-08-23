# Pipeline policy: one recipe, both species, enforced by construction

> Current vs superseded results: **`docs/STATUS_CURRENT_VS_ARCHIVED.md`**.

Status: POLICY. Adopted 2026-08-12. Applies to every Ribo-seq and RNA-seq dataset in this project.

This exists because the project accumulated **13 separate alignment scripts** whose processing
decisions diverged arbitrarily, and because the one step that mattered most was reachable but
optional. The result was a reproducibility scare that turned out to be self-inflicted. The policy is
written to make that specific failure impossible rather than to describe good intentions.

## What actually went wrong (the case this policy is built on)

The training packs' recipe was never lost. It was documented the whole time in
`build_psite_target.py`:

> "rRNA/tRNA/miRNA transcripts + cross-gene multimappers removed; within-gene isoform multimappers
> are retained and each counted (RiboCode default = multimapper posture A)"

and implemented in `filter_tx_heldout.py`. `prepare_from_bams.py` lists that filter as step 1 of its
own documented flow. But it was exposed as `--ncrna-tx ... optional`, so every dataset processed in
2026-08 silently skipped it, came out with ~5% excess P-sites, and the resulting 12.6% divergence was
initially misread as "the training data cannot be regenerated".

The lesson is narrow and general: **an optional flag for a mandatory step is a latent reproducibility
bug.** Documentation did not prevent it. A default did.

## The rules

### 1. One pipeline. Species is configuration, not code.

There is exactly ONE path from reads to pack:

```
FASTQ -> [adapter/UMI per samplesheet]
      -> STAR --alignEndsType EndToEnd --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM
      -> sort + index
      -> filter_tx_heldout.py   (ncRNA drop + cross-gene read drop)     <- CANONICAL, not optional
      -> RiboCode metaplots -> RiboCode -g -> per-nt *_psites.hd5
      -> prepare_pack.py  -> 7-file pack (+ ORF track, --kozak none)
```

RNA-seq differs only in that it goes to `rnaseq_coverage.py` instead of the RiboCode branch, and
never uses `--maximum-length` or `--discard-untrimmed`.

Mouse and human differ ONLY in which files are passed: STAR index, RiboCode annotation, tx2biotype,
ncRNA tx list, tx->gene map, universe. Never in which steps run. If a step is right for one species
it is right for both; if it is not, that is a finding to document, not a fork to add.

**The only legitimate per-dataset variation** is what the library actually is:

| knob | why it varies | example |
|---|---|---|
| adapter | submitters differ | TruSeq, poly-A `A{6}`, N-pad `--trim-n`, CAR-T-specific |
| UMI | protocol-dependent | GSE304796 has a 12 nt 5' UMI; nothing else does |
| read layout | SE vs PE | GSE208041 RNA is paired |

Everything else is fixed. Adapter choice must be VERIFIED against the reads (see rule 3), never
assumed from the submitter's description.

**Two STAR flags are absolute for Ribo-seq** (user instruction, 2026-08-12):
`--alignEndsType EndToEnd` and `--outFilterMultimapNmax 1`. These are not tunable per dataset.

EndToEnd is not a stylistic preference. Soft-clipping shifts the inferred P-site, which is the exact
quantity the model predicts. Measured on HUVEC against its reference pack, same tissue and same
pooling code, alignment recipe the only difference:

| alignment | per-nt r vs reference |
|---|--:|
| EndToEnd + bowtie2 depletion | **0.935** |
| soft-clipped, no depletion | 0.583 |

Soft-clipped alignments do not merely add noise, they move the signal.

**Audit 2026-08-12, RESOLVED 2026-08-13:** Janich, Wang and GSE243134 Ribo alignments all LACKED
EndToEnd, so the entire mouse-liver 3x3 was built on soft-clipped data. `align_gse243134.sbatch`
predated the rule; that recipe was then deliberately copied to Janich and Wang for factorial
consistency, which propagated the violation rather than containing it. Copying an existing recipe
for internal consistency does not make it correct -- check it against the rules first.

All three were re-aligned through `scripts/riboseq_align.sbatch` and the full chain rebuilt
(40 alignments -> 3 pools -> 9 packs -> 18 dumps -> rescore). Final-recipe pooling removed 3.5-4.9% of
P-sites and 3.8-4.3% of observed ORF calls; every qualitative conclusion held, including the
headline matched-vs-mismatched RNA effect (+0.0030 -> +0.0035 F1). Off-recipe outputs are archived
at `results/_archive_offrecipe_2026_08_13/` rather than deleted, so the comparison stays auditable.

Human datasets were compliant on EndToEnd but their packs were built without the ncRNA + cross-gene
filter. All 11 human Ribo runs were RE-ALIGNED canonically on 2026-08-13 (job 36765546, 11/11, 0
failed; filter drop 9.1-10.1%, matching the mouse rate of 8.86% median over 71 runs). The three human
PACKS have not yet been rebuilt from those alignments, because doing so changes every CAR-T-derived
number and that is a decision rather than a cleanup -- see "Consequences for existing artifacts".

### 2. Deviation by declaration, never by omission.

Enforced in code as of 2026-08-12: `prepare_from_bams.py --ribo-bam` now REQUIRES `--ncrna-tx` and
`--tx-to-gene`. Skipping the canonical filter requires `--no-ncrna-filter`, which is recorded in the
pack provenance. A silent skip is no longer possible.

Generalise this: any step that the final recipe includes must be on by default. An option to
turn it off is fine; an option to turn it on is not.

### 3. Inputs are listed, not globbed.

Added 2026-08-13, from a near-miss. The leukocyte rebuild selected RNA inputs with
`ls data/heldout_bam/<arm>/*.bam | grep -v bai`. Those directories are COMBINED: they hold the arm's
RNA alignments, the old off-recipe Ribo alignments, and for the GSE120762 arms genome-coordinate
BAMs beside transcriptome ones. The glob would have fed Ribo runs in as RNA coverage and
double-counted every sample present as both BAM types. It fails the usual way: a pack of exactly the
right shape, with wrong numbers in it.

- **Name the runs.** Per-dataset run lists (`data/heldout_refs/<arm>_rna_runs.txt`,
  `data/liver3x3/<ds>_{ribo,rna}_runs.txt`) are the input contract. A directory listing is not.
- **Match the specific suffix.** `*.toTranscriptome.out.bam`, never `*.bam | grep -v bai`. Reference
  space is not something to infer from what happens to be in a folder.
- **Assert the disjointness that matters.** The rebuild aborts if any Ribo run appears among the RNA
  inputs. Cheap, and it makes the assumption visible.
- A glob that is correct only because of what currently sits in the directory is not a guarantee,
  and it silently stops being correct when someone adds a file.

Scan by artifact, not by the arm you noticed (`feedback_blast_radius_scan_by_artifact`): checking
every RNA BAM directory and every script using the idiom is what showed `mouse_janich_liver` clean,
confirmed the 3x3 unaffected, and bounded the problem to two directories.

### 3b. A new script must inherit its environment from the working one it copies.

Added 2026-08-13, after this failure mode recurred three times in one day.

`scripts/train_canon.sbatch` was written fresh, alongside the working `train_loto_union.sbatch`,
reusing its ARGUMENTS but not its ENVIRONMENT. Both trainings died in under 70 s:

| omission | symptom | why it was silent |
|---|---|---|
| no `RIBO_ONEHOT_FASTA` | `KeyError: 'ENST00000869382.1'` | `dataset.py` DEFAULTS it to `fibroblast_universe.fa` (36,668 tx). The final-recipe packs are on the union universe (84,472). It failed on the first transcript outside the Fibroblast subset -- 47,817 of 84,472 were missing |
| SIF hardcoded to rinalmo | `ModuleNotFoundError: mamba_ssm` | `mamba_ssm` ships only in the orthrus image; the mixer must select the SIF |

Both are environment, not arguments, and neither appears in `args.json`, so the "every hyperparameter
was read back out of the deployed run" check that the script's own header boasts about would never
have caught them.

The same shape recurred earlier the same day: an RNA glob copied from a script whose directory
happened to be clean (rule 3), and a recipe copied from `align_gse243134.sbatch` to Janich and Wang
"for factorial consistency" that propagated a missing `EndToEnd` (rule 1). In every case the copy was
of the visible part and the divergence was in the part that has a default.

- **Copy the whole invocation**, `--env` lines included, then change what must change.
- **A default is not a decision.** Anything with a fallback (`RIBO_ONEHOT_FASTA`, `RIBO_PACK_SUFFIX`,
  `RIBO_ORF_TRACK`) must be set explicitly in a new script, even when the default is correct today.
- **Preflight the cheap invariants.** `train_canon.sbatch` now verifies the one-hot FASTA covers the
  pack's `tx_order` and that the SIF exists, before requesting work. Tested in both directions: it
  passes on `union_universe.fa` and aborts on `fibroblast_universe.fa`.

### 4. Reproducibility is enforced, not intended.

| mechanism | status | what it prevents |
|---|---|---|
| canonical filter required or declined | **DONE** | the silent-skip failure above |
| `filter_tx_heldout` refuses coordinate-sorted input | **DONE** | the false 0.30% filter reading |
| explicit run lists + Ribo/RNA disjointness assert | **DONE** (rebuild path) | Ribo alignments pooled as RNA coverage |
| provenance written for every pack | partial: 31 packs carry `inputs_known: false` | packs whose recipe is unknowable |
| recipe version + hash in provenance | TODO | comparing two packs built by different recipes |
| single entry point; new `align_*.sbatch` is a lint failure | TODO | the 13-script sprawl |
| end-to-end fixture hashed in CI | TODO | a recipe change silently altering output |
| adapter verified against reads, recorded | practice, should be automated | trimming for an adapter that is not there |

## RULE 5: the ORF-length floor is set BY ASSAY (30 AA tryptic, 7 AA MHC)

**The rule is 90 nt = exactly 30 amino acids.** RiboCode's `ORF_length` excludes the stop codon
(verified: 1911 nt / 3 = 637 aa for a 637-aa AAseq), so the two are the same number.

| assay | floor | why |
|---|--:|---|
| **tryptic whole-cell lysate** | **30 AA** | tryptic peptides come from digesting whole proteins; sub-30-aa ORFs are not the population the assay samples. Same floor as the ORF-call track. |
| **MHC / HLA immunopeptidomics** | **7 AA** | HLA-I peptides are 8-11 aa. A 30-aa floor would exclude the microproteins the assay exists to detect. |

**Classify by the SEARCH PARAMETERS, never by a name list.** `search_enzyme_name_1 = trypsin` with
`num_enzyme_termini = 2` is tryptic; `nonspecific` with `termini = 0` is MHC. In this project that
makes the 12 mouse macrophage populations and A549 TRYPTIC (30 aa), and HBL-1, SU-DHL-4, DoHH2,
THP-1 and B721.221 MHC (7 aa).

**Enforced, not just documented.** `pgx/build_dbs.py --assay {tryptic,mhc}` SETS `--min-aa` (30/7)
and aborts if a conflicting `--min-aa` is passed; omitting both warns loudly and falls back to 7 for
backward compatibility. `pgx/filter_db_min_aa.py` carries the same check.

**ORF-call evaluation applies 30 aa symmetrically, at scoring time.** Calls are made permissively
(`ribocode_dropin.py --min_aa 5`) and the floor is applied by `compare_dropin_calls.build_loader`
(`min_len=90`), which passes it to model calls and real Ribo-seq calls through one code path -- the
`is_pred` flag gates only the enrichment filter, never the length. `compare_dropin_ctg` and
`compare_dropin_allalt` default to 90 too. **No scorer may override it.** Filtering at scoring rather
than at calling is deliberate: it guarantees model and reference are cut at the same place, which a
call-time threshold cannot.

**What the floor costs, measured per assay** (`proteogenomics/scripts/orf_length_audit.py`, novel
peptides at 1% class-specific FDR whose supporting ORFs are ALL sub-30-aa):

| arm | MOUSE macrophage (tryptic) | HUMAN (mostly HLA-I) |
|---|--:|--:|
| model (Poisson) | **0 / 295 (0.0%)** | **14 / 52 (26.9%)** |
| model (theta=1) | -- | **25 / 76 (32.9%)** |
| null_atg | 19 / 453 (4.2%) | 22 / 72 (30.6%) |
| **CPAT** | -- | **0 / 39 (0.0%)** |
| **CPC2** | -- | **0 / 55 (0.0%)** |

**The two assays genuinely need different floors and neither number transfers to the other.** On
tryptic macrophage data the 30-aa floor costs the model nothing. On human HLA-I it would cost ~30%,
and because CPAT/CPC2 have 0.0% short-only (their databases are 0-2% short by construction) the loss
is NOT symmetric: **D12's "the model wins all three HLA-I immunopeptidomes" would become "wins 2 of
3"** (HBL-1: model 14 -> 8 vs CPC2 10 -> 10). That is precisely why HLA-I stays at 7 aa.

**Changing the floor requires a RE-SEARCH, not a post-hoc filter.** Shrinking the database changes
decoy competition, the FDR score cutoff, and which peptide wins rank 1. `orf_length_audit.py`'s
`kept_peptides` column is that post-hoc approximation and is not a substitute.

## Terminology: "canonical" is reserved for BIOLOGY. The pipeline is the FINAL RECIPE.

Extended 2026-08-15. The earlier version of this rule only banned "non-canonical" for provenance; it
then used "canonical" for the recipe on the very next line, which is the same collision one word
shorter.

**"Canonical" in this project means one thing: the annotated CDS.** `canonical CDS`, `canonical
peptides`, `canonical PSMs`, `canonical ORFs`, and their negation `non-canonical` (uORF, dORF, novel,
internal, Overlap_*) appear in every per-class table, ceiling figure and proteomics report.

**The processing pipeline is the FINAL RECIPE.** A pack, alignment, retrain or reference produced by
it is **on-recipe** or **final-recipe**; anything else is **off-recipe**.

  on-recipe / final-recipe   produced by the pipeline in rule 1, deviations = none
  off-recipe                 any declared or undeclared deviation from it

Reusing the biological term for a provenance property produces sentences that parse and mislead:
"non-canonical F1 from a non-canonical pack", or "the canonical retrain improved canonical CDS".
The second is worse because it looks like it says something.

**Directory names are NOT renamed.** `data/packed_canon*`, `results/loto_canon/`,
`results/mouse_liver_3x3_canon/`, `_canonical_<T>/` and the `_canon` pack suffix are load-bearing in
scripts, `RIBO_PACK_SUFFIX`, provenance JSONs and every results path -- the same argument that keeps
the dataset registry's internal keys unrenamed (`data/dataset_registry.tsv`). **`_canon` on a path
means "final recipe".** Rename prose, never paths.

**Renaming this in bulk needs a negative lookbehind.** `\bcanonical reference\b` matches INSIDE
`non-canonical reference`, because the hyphen is a word boundary -- a plain substitution turns it
into "non-final-recipe reference", inverting the project's most load-bearing term. Use
`(?<!non-)\bcanonical X\b`, and assert the count of `non-canonical` and `canonical CDS|peptide|PSM|
identification|ORF|cost|Kozak` is unchanged before writing any file.

## RESOLVED 2026-08-12: the training data IS reproducible

The reproducibility scare is closed. The recipe comparison below was run on HUVEC; the ACROSS-TISSUE
answer is in the next section, and it is not 3.04%.

> **Do not quote 3.04% as "the" reproducibility number.** It is HUVEC's, and HUVEC turned out to be
> the BEST of the eight tissues. Across all eight the median is 3.86% and the worst is 11.66%
> (2026-08-13, next section). The table below is a RECIPE comparison -- which alignment recipe
> reproduces the training data -- and HUVEC is the fixed substrate that makes those four rows
> comparable to each other. It is not an estimate of typical reproducibility.

| recipe (all on HUVEC) | divergence from `packed_union_HUVEC` | per-nt r |
|---|--:|--:|
| **CANONICAL** (EndToEnd + mm1 + filter on query-grouped bams) | **3.04%** | **0.957** |
| A: my reconstruction (bowtie2 + snoRNA + EndToEnd, no filter) | 11.69% | 0.935 |
| B: minimal (no bowtie2, no EndToEnd, no filter) | 12.82% | 0.583 |
| C: "with filter" -- filter fed COORDINATE-SORTED bams, silently inert | invalid | -- |
| reference: Janich pack via the current pipeline | 1.69% | -- |

Applying the documented recipe collapses the divergence 4x, into the same regime as the Janich gate.
The residual is consistent with RiboCode `metaplots` re-deriving periodic read lengths per sample,
which is the known irreducible stochastic component.

**The 12.6% reported earlier was recipe drift, not a property of the data.** Three compounding
mistakes produced it, all on the reconstruction side:

1. A recipe was reconstructed WITHOUT first reading `build_psite_target.py`, which documented it.
2. Stages the original never had were added (bowtie2 contaminant depletion, snoRNA depletion).
3. When the filter was finally tested, it was fed COORDINATE-SORTED bams. It requires query-name
   grouping; without it, cross-gene detection silently collapses and it reported 0.30% instead of
   the true ~9.5%. That wrong number was then used to dismiss the very step that mattered.

Guards now make (3) impossible: `filter_tx_heldout.py` refuses coordinate-sorted input, and
`prepare_from_bams.py` refuses to pool Ribo bams without either the filter arguments or an explicit
`--no-ncrna-filter` recorded in provenance. Both fired in production within a day of being added.

## ACROSS-TISSUE 2026-08-13: reproducibility is heterogeneous, and one tissue was not enough

The section above establishes WHICH RECIPE reproduces the training data, on one tissue. This one
establishes HOW WELL it does so across all eight. Every Chothani Ribo run (74) was re-aligned
canonically and re-pooled per tissue, then compared against the pack the model actually trained on.
Same metric as the HUVEC test, `sum|new - old| / sum(old)` over per-nt P-sites; the generalised
script reproduces the HUVEC number exactly (3.04%, r = 0.9567), which is why the other seven are
trustworthy.

| tissue | training P-sites | \|divergence\| | per-nt r |
|---|--:|--:|--:|
| Hepatocytes | 1,159,703,745 | 3.16% | 0.9679 |
| Fibroblast | 1,093,425,533 | 3.71% | 0.9374 |
| HCAEC | 714,826,272 | 3.20% | 0.9185 |
| Fat | 581,557,882 | **11.66%** | 0.8442 |
| ES | 508,576,671 | 4.01% | 0.9800 |
| VSMC | 486,734,365 | 5.85% | 0.8373 |
| HA_EC | 173,184,690 | 8.01% | **0.4320** |
| HUVEC | 95,959,146 | **3.04%** | 0.9567 |

**Median 3.86%, range 3.04-11.66%, per-nt r 0.432-0.980.**
Full table: `results/chothani_regeneration/divergence_by_tissue.{tsv,md}`; figure:
`figures/S_reproducibility/`.

Three things a reader of this policy needs:

1. **HUVEC is the best of the eight**, so the single-tissue number understated the typical case and
   understated the worst by nearly 4x. Quote the range or the median, never 3.04% alone.
2. **Divergence and correlation are partly decoupled.** Fat moved the most signal (11.66%) but moved
   it coherently (r = 0.844); HA_EC moved less (8.01%) but incoherently (r = 0.432). Report both.
3. **There is no depth relationship.** Hepatocytes (1.16B P-sites) and HUVEC (96M) both land at
   ~3.1% across a 12x range. Do not explain a tissue's divergence by its depth.

HA_EC (r = 0.432, every other tissue >= 0.837) is an unexplained outlier. Four hypotheses were tested
and rejected: uniform P-site shift (best-lag correlation peaks at lag 0), defective old pack (all 8
old packs' start-codon metagenes peak at +0 nt, frame0 0.732-0.825, HA_EC mid-range), changed P-site
calling (old vs new frame0 agree within tissue; HA_EC +0.028 vs HUVEC +0.027, and HUVEC scores
0.957), and metaplots offset disagreement (recorded as a prediction, then falsified -- ES has offset
disagreement and the HIGHEST r at 0.980). **No causal story about HA_EC goes in the manuscript until
it is settled.**

## What the filter actually costs

Measured on correctly ordered input, the ncRNA + cross-gene filter removes **9-22% of records**,
dominated by cross-gene paralog ambiguity rather than ncRNA:

| dataset | ncRNA dropped | cross-gene dropped | total |
|---|--:|--:|--:|
| mouse liver (typical) | 2.0% | 8.8% | ~10.8% |
| human HUVEC | 0.3% | 9.4% | ~9.7% |
| GSE120762 SRR7956053 | **49.9%** | 2.2% | 52.2% |

That last row is a library that is roughly half rRNA/tRNA/miRNA; the filter is the only thing
removing it, since the mouse path never had bowtie2 depletion. It still yields 9.1 M usable reads,
comparable to its peers.

Note the orthogonality that makes this step necessary at all: `--outFilterMultimapNmax 1` filters on
GENOMIC loci, while `--quantMode TranscriptomeSAM` emits one record per compatible TRANSCRIPT.
Verified empirically -- 58% of reads mapped uniquely to the genome yet averaged 2.68 transcript
records each. mm1 therefore never removes paralog ambiguity; only the cross-gene filter does.

## Consequences for existing artifacts

Audited per script, 2026-08-13. The three groups are in genuinely different states, and an earlier
version of this section got two of them wrong by asserting from memory instead of reading the
scripts:

| artifact | EndToEnd | mm1 | cross-gene filter | status |
|---|:--:|:--:|:--:|---|
| mouse-liver 3x3 (janich, wang, gse243134) | NO | yes | NO | **was off-recipe -> REBUILT + rescored** |
| human (THP-1/GSE208041, GSE39561, CAR-T) | yes | yes | **NO** | **off-recipe, open** |
| mouse leukocyte (gse120762 x2, gse155087) | yes | yes | yes | **on-recipe all along** |

**Mouse-liver 3x3 -- resolved.** Measured cost of the off-recipe alignment, canonical vs off-recipe
pooling of the same libraries:

| dataset | off-recipe P-sites | canonical | delta |
|---|--:|--:|--:|
| janich | 109,996,969 | 106,115,999 | -3.5% |
| wang | 78,257,663 | 74,504,278 | -4.8% |
| gse243134 | 150,531,643 | 143,221,071 | -4.9% |

Positional shifts exceed these totals, because EndToEnd changes where the 5' end sits and the P-site
offset is measured from it. Rebuilt and rescored 2026-08-13: observed calls fell 3.8-4.3%,
`pred_preddepth` F1 fell ~0.01, and every qualitative conclusion held. Off-recipe outputs archived at
`results/_archive_offrecipe_2026_08_13/`.

**Human -- still open, and worse than "missing the filter".** `process_thp1.sbatch` and
`process_cart_gse304796.sbatch` both pass EndToEnd and mm1 but never call `filter_tx_heldout`, so
ncRNA transcripts and cross-gene reads (~9.7% of records) remain.

CAR-T (GSE304796) has a second and larger problem, found 2026-08-13 while scoping the fix. That
pipeline runs `umi_tools dedup`, but it runs it on the GENOME BAM only. The transcriptome BAM is
copied straight from STAR scratch:

```
cp "$SCR/${RUN}.Aligned.toTranscriptome.out.bam" "$OUT/${RUN}.toTranscriptome.bam"
```

and the pack is built from exactly that file (`_work/psites/SRR34901743.toTranscriptome_psites.hd5`).
So the deduplication was performed on a BAM the pack never sees. Duplicate rate from the pipeline's
own `read_counts.txt`, genome arm:

| run | aligned | deduped | PCR duplicates |
|---|--:|--:|--:|
| SRR34901743 | 33,114,136 | 23,888,001 | 27.9% |
| SRR34901744 | 43,372,505 | 30,637,962 | 29.4% |
| SRR34901745 | 45,122,621 | 31,747,439 | 29.6% |

**The CAR-T pack's P-sites therefore include roughly 29% PCR duplicates**, on top of the ~9.7% the
filter would remove. Duplicates are not distributed evenly along a transcript -- they pile at
specific positions -- so this distorts profile shape, not just depth. Any CAR-T number should be
treated as provisional until the pack is rebuilt; the human codon-occupancy value (0.548) is
computed over THP-1 and CAR-T together and is flagged in the tutorial accordingly.

GSE208041 and GSE39561 have no UMIs at all (GSE208041's submitter stripped them before deposit), so
for those two the filter is the only gap.

**Re-aligned 2026-08-13 (job 36765546, 11/11, 0 failed).** `scripts/riboseq_align.sbatch` grew `umi`
and `cutadapt_extra` samplesheet columns -- rule 1 already listed UMI as a legitimate per-dataset
knob, the script simply had not implemented it. Order is forced by the tools: STAR ->
`filter_tx_heldout` (needs query-grouped input) -> `umi_tools dedup` (needs coordinate-sorted +
indexed), which the filter leaves behind.

Both predictions checked out, which is what makes the new alignments trustworthy rather than merely
new:

| | filter dropped | UMI duplicates | genome-arm reference |
|---|--:|--:|--:|
| gse208041 (5 runs) | 9.12-10.06% | n/a (no UMI) | -- |
| gse39561 (3 runs) | 9.29-9.55% | n/a (no UMI) | -- |
| cart SRR34901743 | 10.05% | **27.53%** | 27.9% |
| cart SRR34901744 | 9.91% | **29.05%** | 29.4% |
| cart SRR34901745 | 9.98% | **29.26%** | 29.6% |

Deduplicating in TRANSCRIPTOME space reproduces the genome-arm rate to within 0.4 points, so the
~29% duplicate estimate was right and the two coordinate spaces agree. The filter rate (9.1-10.1%)
also matches the mouse Chothani rate (8.86% median over 71 runs), i.e. the two species are now on one
recipe by measurement and not just by intent.

REMAINING: the three human PACKS have not been rebuilt from these alignments yet. Doing so will
change every CAR-T-derived number (including the human codon-occupancy 0.548, which pools THP-1 with
CAR-T), so it is held pending an explicit decision rather than done as a side effect.

A cheaper design for future UMI datasets, worth recording: dedup the GENOME bam (one record per read
under `--outFilterMultimapNmax 1`) and propagate surviving read names to the transcriptome bam. The
run above deduplicated 600-621M transcriptome records to resolve only ~33-45M molecules, roughly 15x
the necessary work, because the transcriptome bam carries ~10 records per read across isoforms. The
original CAR-T pipeline was in fact one propagation step away from correct: it deduplicated the
genome bam and then never carried the result across.

**Mouse leukocyte -- not off-recipe; the earlier claim here was wrong.** These were built by
`scripts/heldout/align_and_filter.sh`, which already applies EndToEnd, mm1, TranscriptomeSAM and
`filter_tx_heldout`. The filter received query-grouped input: STAR's `Aligned.toTranscriptome.out.bam`
is read-ordered regardless of `--outSAMtype BAM SortedByCoordinate`, which sorts only the genome BAM.
Their sole difference from canonical is two tightened STAR filters
(`--outFilterMismatchNoverLmax 0.05`, `--outFilterMatchNminOverLread 0.7` vs the defaults 0.3 / 0.66)
-- two of the four parameters the 13-arm HUVEC sweep found immaterial to periodicity.

Rebuilding them anyway turned the assumption into a measurement, and it is the cheapest kind of
reproducibility gate. GSE120762 LPS: P-sites 38,481,730 -> 38,636,019 (**+0.40%**), RNA coverage
**byte-identical** (`np.array_equal` true). The sign is informative: canonical is the LOOSER recipe
here, so it keeps slightly more reads -- consistent with the sweep, where tightening cost up to 5.3%
of reads for no gain in frame purity. That 0.40% sits alongside HUVEC 3.04% and the Janich 1.69%
gate as a third independent regeneration number.

The general lesson is the one this document already makes elsewhere: read the script, do not recall
it. "These packs are off-recipe" was inferred from the group they were filed under, not verified,
and it was wrong for three of the six artifacts named.

## For reviewers of this project

A pack is trustworthy iff its `provenance.json` records (a) every input file, (b) the recipe version,
and (c) any declared deviation. A pack with `inputs_known: false` is a pack whose numbers cannot be
regenerated, and that is now visible rather than silent -- which is the minimum bar, not the goal.
