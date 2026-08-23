# Human Ribo-seq datasets added 2026-08 (THP-1 x2, CAR-T)

Alignment QC as processed. ARM is determined by which recipe the sample actually took:
a bowtie2 contaminant-depletion log means the RPF path, its absence means the RNA path.
GEO/SRA `library_strategy` is wrong for GSE39561 and GSE304796 (both say RNA-Seq; both are Ribo-seq),
so arm was assigned from the submitters' sample titles, not from SRA metadata.

## GSE208041_thp1

**ALL SAMPLES HERE ARE VEHICLE (no-treatment). Verified against ENA 2026-08-15.** GSE208041
(PRJNA858047) contains THREE conditions x 5 reps -- Vehicle, LPS, LPSDex. The runs used here are
exactly the 5 Vehicle Ribo-seq reps and the 4 Vehicle RNA-seq reps; no LPS or LPSDex sample is in the
pack. Treated runs NOT used: Ribo SRR20106154-163, RNA SRR21228007-014.

The treatment condition was previously recorded NOWHERE in this repo -- only the RPF/RNA arm -- so
nothing stated that the reference was untreated. Added because that is the same provenance gap that
left Janich/Wang/Chothani without accessions.

| run | arm | input reads | STAR unique % | rRNA/tRNA depleted % |
|---|---|---|---|---|
| SRR20106164 | RPF (Vehicle rep5) | 52,023,766 | 60.17 | 67.72 |
| SRR20106165 | RPF (Vehicle rep4) | 32,692,290 | 62.57 | 66.32 |
| SRR20106166 | RPF (Vehicle rep3) | 47,764,149 | 58.58 | 72.06 |
| SRR20106167 | RPF (Vehicle rep2) | 45,770,613 | 57.97 | 71.63 |
| SRR20106168 | RPF (Vehicle rep1) | 48,016,958 | 61.45 | 70.86 |
| SRR21228005 | RNA (Vehicle rep2) | 81,472,331 | 90.53 | n/a |
| SRR21228006 | RNA (Vehicle rep4) | 82,273,565 | 91.64 | n/a |
| SRR21228015 | RNA (Vehicle rep1) | 74,560,134 | 90.71 | n/a |
| SRR21228016 | RNA (Vehicle rep3) | 91,298,138 | 93.07 | n/a |

## GSE304796_cart

| run | arm | input reads | STAR unique % | rRNA/tRNA depleted % |
|---|---|---|---|---|
| SRR34901743 | RPF | 52,426,665 | 63.16 | 21.55 |
| SRR34901744 | RPF | 67,391,487 | 64.36 | 18.87 |
| SRR34901745 | RPF | 71,546,979 | 63.07 | 20.87 |
| SRR34901749 | RPF | 65,046,960 | 72.74 | 3.85 |
| SRR34901750 | RPF | 59,962,970 | 71.87 | 3.93 |
| SRR34901751 | RPF | 60,436,502 | 73.33 | 3.26 |

## GSE39561_thp1

| run | arm | input reads | STAR unique % | rRNA/tRNA depleted % |
|---|---|---|---|---|
| SRR525273 | RPF | 50,720,051 | 29.50 | 43.28 |
| SRR525274 | RPF | 53,374,529 | 32.66 | 44.15 |
| SRR525275 | RPF | 57,018,418 | 32.03 | 39.16 |

## Notes

- GSE39561 has NO matched RNA-seq **of its own**, and it is the PREDICTION arm that survives, not the
  observed one. CORRECTED 2026-08-14 -- the previous wording here said the opposite and was wrong.
  Verified against the pack itself (`data/packed_heldout_human_gse39561/`):
    * `target_counts.npy` sums to **0** -- the library has no detectable 3-nt periodicity, so
      metaplots selects zero read lengths and NO P-sites can be placed. There is therefore no `real`
      arm and no `pred_obsdepth` arm.
    * `coverage.npy` sums to 177,105,838,826 with `n_rna_samples: 4` -- **byte-identical to
      GSE208041's**, on the same `tx_order`. It BORROWS GSE208041's RNA and universe.
  So only the standalone `pred_preddepth` arm is defined, and because the borrowed RNA and universe
  are the sole model inputs here, those predictions are IDENTICAL to GSE208041's. GSE39561 is a test
  of the DATA, not an independent test of the model. `human_02_dump.sbatch` and
  `score_human_orf_calls.py` both already state this correctly.
- GSE39561 unique% (29.5-32.7) is low but expected: 2012-era short single-end reads with a poly-A
  trim, after ~40% rRNA/tRNA depletion.
- CAR-T RPF carries 19-22% rRNA/tRNA vs 3-4% for its RNA arm, the normal footprint-vs-RNA split.
- CAR-T RNA is rRNA-depleted TOTAL RNA, not poly(A)-selected. The poly(A) gate
  (feedback_rnaseq_must_be_polya) is NOT yet applied to it; do not build a universe from this RNA
  until that check runs.

## Poly(A) gate (2026-08-10)

Salmon quant against the decoy-aware GENCODE v49 index, composition read by
`scripts/polya_gate_report.py`. Full table: `results/polya_gate/human_rna_polya_gate.tsv`.
GSE208041 was run alongside as a control -- an arm expected to pass makes a failure interpretable
instead of ambiguous.

| dataset | protein_coding TPM | lncRNA TPM | top-20 concentration | salmon mapping | ncRNA_host TPM |
|---|---|---|---|---|---|
| GSE208041 (control) | 92.1-92.5% | 2.6-2.8% | 21.9-26.2% | 77.6-80.9% | 15.3% |
| GSE304796 (CAR-T) | 54.4-55.9% | 38.8-40.4% | 42.6-44.8% | 55.9-58.3% | **43.1%** |

**GSE304796 RNA is confirmed NOT poly(A)-selected.** A single gene, ENSG00000282885, carries ~18-20%
of the library on its own. The four dominant transcripts are unnamed lncRNAs and three of them carry
the GENCODE tag `ncRNA_host` -- snoRNA host genes. In a Ribo-Zero library the abundant snoRNAs map
into their host transcripts, and because GENCODE types those hosts as lncRNA a biotype filter cannot
remove them. This is the documented total-RNA failure mode
(`feedback_rnaseq_must_be_polya`), and the name-prefix heuristic missed it entirely: these loci have
no gene symbol, so `ncRNA_host` is the discriminator, not the name.

**But the universe does not collapse, so this is not the GSE243134 case.** Excluding the 3,389
`ncRNA_host` genes leaves 35,825 transcripts at TPM>=1, against 36,070 for the clean control. The
GSE243134 liver universe crushed to 3,321. The contamination here is concentrated TPM *mass*, not
lost transcript *coverage*.

### Decision

GSE304796 may define a universe **only with `ncRNA_host` genes excluded**
(`data/ncrna_host_genes_v49.txt`, 3,389 genes from the v49 GTF). Rationale:

- Membership, not quantity, is what a universe needs, and 35,825 surviving transcripts is a full one.
- Excluding 43.1% of TPM mass renormalises the remainder upward, so a TPM>=1 cut on this library is
  effectively stricter than the same cut on a clean one (~TPM>=0.57 equivalent). That errs
  conservative, which is the safe direction for a gate.
- The distortion is in the TPM values, not in which transcripts are detected.

This judges universe definition only. It does not block the same RNA from being used as model
coverage input: salmon EM plus length normalisation and STAR mm1 plus ncRNA filtering disagree
strongly on multi-copy loci, and a TPM-contaminated library can still yield clean pack coverage
(`feedback_salmon_tpm_vs_pack_coverage`).

## CORRECTION (2026-08-10): snoRNA depletion does NOT explain the host-gene TPM

The read-level snoRNA/scaRNA depletion stage was added on the hypothesis that abundant snoRNAs were
mapping into their lncRNA host transcripts and inflating `ncRNA_host` TPM. **That hypothesis is
wrong, and the measurement says so.**

Depletion rate against `snorna_grch38_v49`, per sample:

| dataset | arm | reads removed |
|---|---|---|
| GSE208041 | RPF | 0.58-0.87% |
| GSE208041 | RNA | 0.05-0.06% |
| GSE304796 | RPF | 0.14-0.16% |
| GSE304796 | RNA | 0.32-0.34% |
| GSE39561 | RPF | 0.01% |

And the quantity it was supposed to fix did not move at all:

| sample | ncRNA_host TPM before | after |
|---|---|---|
| SRR21228005 | 15.3% | 15.3% |
| SRR21228006 | 15.2% | 15.2% |
| SRR21228015 | 15.4% | 15.4% |
| SRR21228016 | 15.9% | 15.9% |

Unchanged to the decimal. Removing every snoRNA read costs a third of a percent of the library and
leaves host-gene TPM exactly where it was, so **the TPM sitting in `ncRNA_host` genes is not
snoRNA-derived reads.** It is host-transcript signal: a Ribo-Zero total-RNA library retains
non-polyadenylated and nascent transcripts that poly(A) selection removes, and the SNHG-type host
lncRNAs are exactly that class. The 43.1%-vs-15.3% split between GSE304796 and the poly(A) control is
a genuine library-composition difference, not contamination that can be filtered out of the reads.

A second plausible-sounding explanation was also checked and rejected: snoRNAs are NOT simply
intronic and therefore absent from the quantified mature transcript. 252 of 972 snoRNA intervals
(25.9%) overlap a host-gene exon. The reason depletion has no effect is the read count, not the
annotation geometry.

### What this changes

- **The universe mitigation stands unchanged**: GSE304796 may define a universe only with
  `ncRNA_host` genes excluded. That operates on the quantification, which is where the problem
  actually lives, and 35,825 transcripts survive at TPM>=1 (control: 36,070).
- **The snoRNA stage is kept** -- it does what it says (snoRNA reads are now removed at read level,
  uniformly across all three human datasets, which is the consistency it was asked for) and costs
  0.01-0.87% of reads. It is simply not a fix for the universe, and must not be cited as one.

## OFF-RECIPE status (flagged 2026-08-13)

The GSE208041, GSE39561 and GSE304796 packs were built through `process_thp1.sbatch` /
`process_cart_gse304796.sbatch`, which DO use `--alignEndsType EndToEnd` but do NOT apply the
canonical ncRNA + cross-gene filter (`filter_tx_heldout`). Measured on comparable human data that
filter removes ~9.7% of records, dominated by cross-gene paralog ambiguity (9.4%) rather than ncRNA
(0.3%).

So these packs are OFF-RECIPE, but less so than the mouse-liver 3x3 was: they have EndToEnd and are
missing only the filter, where the mouse packs were missing both.

**UPDATE 2026-08-13: the ALIGNMENTS are now canonical; the PACKS are not yet rebuilt.** All 11 human
Ribo runs were re-aligned through `scripts/riboseq_align.sbatch` (job 36765546, 11/11, 0 failed) with
the filter applied to query-grouped output. Measured drop: 9.12-10.06% (GSE208041), 9.29-9.55%
(GSE39561), 9.91-10.05% (CAR-T) -- matching the mouse Chothani rate (8.86% median over 71 runs), so
both species are now on one recipe by measurement. Samplesheet + declared deviations:
`data/human_ribo_canon_samplesheet.tsv` and its sibling README.

CAR-T carried a SECOND and larger defect found while scoping that work: `umi_tools dedup` ran on the
GENOME bam while the pack was built from the raw STAR transcriptome bam, so ~29% PCR duplicates
reached the pack (27.9 / 29.4 / 29.6% across the three runs). The canonical re-alignment deduplicates
the bam the pack actually reads and reproduces those rates to within 0.4 points (27.53 / 29.05 /
29.26%), confirming both the defect and the fix.

Consequence for the reported human ORF-call numbers (`results.md`): comparisons among them are
internally consistent, but the absolute F1 values are provisional until the packs are rebuilt from
the new alignments. That rebuild is deliberately HELD rather than done as a side effect, because it
changes every CAR-T-derived number -- including the human codon-occupancy value of 0.548, which pools
THP-1 with CAR-T and is flagged provisional in `figures/C12_codon_occupancy/FIGURE_DATA_INPUTS.md`.
GSE39561 is separately low-value to rebuild: it is already closed as a negative result (40.7% frame
concentration vs a 33.3% random floor).

Note the arms differ in exposure: the RNA arms are unaffected (RNA-seq stays Local and is not
filtered), so the universes -- which are salmon/TPM-derived -- remain valid and do not need
regenerating. Only the Ribo-derived P-sites change.
