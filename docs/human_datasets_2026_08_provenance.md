# Human Ribo-seq datasets added 2026-08 (THP-1 x2, CAR-T)

Alignment QC as processed. ARM is determined by which recipe the sample actually took:
a bowtie2 contaminant-depletion log means the RPF path, its absence means the RNA path.
GEO/SRA `library_strategy` is wrong for GSE39561 and GSE304796 (both say RNA-Seq; both are Ribo-seq),
so arm was assigned from the submitters' sample titles, not from SRA metadata.

## GSE208041_thp1

| run | arm | input reads | STAR unique % | rRNA/tRNA depleted % |
|---|---|---|---|---|
| SRR20106164 | RPF | 52,023,766 | 60.17 | 67.72 |
| SRR20106165 | RPF | 32,692,290 | 62.57 | 66.32 |
| SRR20106166 | RPF | 47,764,149 | 58.58 | 72.06 |
| SRR20106167 | RPF | 45,770,613 | 57.97 | 71.63 |
| SRR20106168 | RPF | 48,016,958 | 61.45 | 70.86 |
| SRR21228005 | RNA | 81,472,331 | 90.53 | n/a |
| SRR21228006 | RNA | 82,273,565 | 91.64 | n/a |
| SRR21228015 | RNA | 74,560,134 | 90.71 | n/a |
| SRR21228016 | RNA | 91,298,138 | 93.07 | n/a |

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

- GSE39561 has NO matched RNA-seq. The model takes RNA coverage as an input, so this dataset
  supports the observed-P-site (`real`) ORF-calling arm only, not the model-prediction arms.
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
