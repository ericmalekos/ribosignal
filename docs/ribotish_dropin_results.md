# Ribo-TISH as a second ORF caller on predicted profiles

Written 2026-09-03. Companion to `scripts/orfcallers/ribotish_dropin.py`.

## Why a second caller

Every ORF-calling result in this project runs through RiboCode. A failure mode of the model and a
failure mode of RiboCode are indistinguishable from inside that setup. Ribo-TISH tests the same
predicted profiles with a different statistic, so a shared conclusion is evidence about the model
rather than about one caller's thresholds.

## How the drop-in works, and why no BAM is needed

`ribotish predict` accepts `--inprofile`, a per-transcript P-site profile, instead of reading bam
files (`run/predict.py:212-220`). Its profile object is indexed 0-based over the mature transcript
(`self.length = trans.cdna_length()`, `zbio/ribo.py:56-66`) with `nhead = ntail = 0`, which is the
same axis the packs use, so a density array drops straight in.

File format, from `predict.py:217` and `ribo.py:338-352`:

    Gid \t Tid \t Symbol \t {pos:count, ...} \t {pos:count, ...}

Column 3 is the TIS profile, written empty because there is no LTM/harringtonine data here, and
column 4 is the RPF profile. Both are sparse python dict literals over 0-based cDNA positions and
are `eval`ed on load.

`ribotish_dropin.py` **imports `build_density` from `ribocode_dropin`** rather than
re-implementing it. The two callers therefore see byte-identical density arrays and the only thing
that differs is the caller.

### Three things that had to be handled

**The GTF must be restricted to exactly the transcripts written.** The bam path is skipped only
when the profile covers every transcript of a gene (`predict.py:381-385`, `load = False`). Passing
the full annotation would leave `load = True` for any gene with an unprofiled transcript.

**A bam argument is mandatory even when unused.** `predict.py:82-84` exits with "No bam file
input!" unless `-t` or `-b` is non-empty. A **deliberately nonexistent** path is passed rather than
an empty bam: `multiRiboGene` sits inside `if load:` and `find_offset` only probes for a sibling
`.para.py`, so the file is never opened, and if `load` ever became True the run would crash loudly
instead of silently scoring every transcript at zero.

**cDNA lengths are checked, not assumed.** Ribo-TISH derives transcript length from GTF exons and
the pack derives it from the sequence it was built on. A silent disagreement would shift the whole
profile. Every transcript is compared and mismatches are dropped and counted. On the deployed
holdout: **70,883 written, 0 absent from the GTF, 0 length mismatches, 0 all-zero**.

## Collapsing: use `--longest`

Without it, 86% of smoke-test calls were `Truncated`, in-frame downstream ATGs sharing the
annotated stop. RiboCode's `*_collapsed.txt` keys on `(gene_id, ORF_gstop)`, one ORF per stop
codon, and `--longest` is the same rule. Both files are kept: `ribotish/longest/` and
`ribotish/all/`.

## Setup

Ribo-TISH 0.2.8, `conda_envs/ribotish`. Deployed checkpoint
`orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`, human hepatocytes
(GSE182371) held out from training, 70,883 scorable transcripts, GENCODE v49 primary annotation
(the same GTF RiboCode's annotation was built from), `--minaalen 5 --fpth 0.05`.

Runtimes on 8 cpus: 42 min (`real`), 87 min (`pred_obsdepth`), 25 min (`pred_preddepth`);
peak RSS 4.8 to 6.8 GB.

## Call counts

| | `--longest` | uncollapsed |
|---|--:|--:|
| `real` | 77,227 | 738,889 |
| `pred_obsdepth` | 89,535 | 776,592 |
| `pred_preddepth` | 85,800 | 762,331 |

Counts are **not** comparable to RiboCode's totals (18,964 / 24,991 / 22,338). Even under
`--longest`, Ribo-TISH enumerates per transcript while RiboCode collapses per gene, and it carries
`Truncated` and `Extended` categories that RiboCode folds into `annotated`.

## Per class, predicted arms scored against the `real` arm

Keyed on `(gene_id, chrom, strand, genomic stop)`. The `real` arm is Ribo-TISH's own calls on the
observed P-sites, so this measures how far the predicted profile diverges from the observed one
through one caller, exactly as the RiboCode drop-in does. It is not an external truth set.

### pred_obsdepth

| class | P | R | F1 | TP | FP | FN | n_real |
|---|--:|--:|--:|--:|--:|--:|--:|
| Annotated | 0.9504 | 0.9947 | 0.9720 | 17,405 | 909 | 92 | 17,497 |
| 5'UTR | 0.3245 | 0.7687 | 0.4563 | 1,957 | 4,074 | 589 | 2,546 |
| Truncated | 0.9373 | 0.9952 | 0.9654 | 1,868 | 125 | 9 | 1,877 |
| Novel | 0.2872 | 0.8159 | 0.4248 | 514 | 1,276 | 116 | 630 |
| 5'UTR:Known | 0.5255 | 0.8586 | 0.6520 | 340 | 307 | 56 | 396 |
| Extended | 0.9661 | 1.0000 | 0.9828 | 371 | 13 | 0 | 371 |
| Internal:CDSFrameOverlap | 0.8844 | 0.8947 | 0.8895 | 306 | 40 | 36 | 342 |
| 3'UTR:CDSFrameOverlap | 0.9163 | 0.9544 | 0.9350 | 230 | 21 | 11 | 241 |
| Internal | 0.1697 | 0.3544 | 0.2295 | 56 | 274 | 102 | 158 |
| Truncated:Known | 0.9304 | 0.9800 | 0.9545 | 147 | 11 | 3 | 150 |
| 3'UTR | 0.0489 | 0.4771 | 0.0887 | 52 | 1,011 | 57 | 109 |
| Extended:Known | 0.9839 | 0.9839 | 0.9839 | 61 | 1 | 1 | 62 |
| 5'UTR:CDSFrameOverlap | 0.4400 | 0.8049 | 0.5690 | 33 | 42 | 8 | 41 |
| Extended:CDSFrameOverlap | 0.9600 | 0.9600 | 0.9600 | 24 | 1 | 1 | 25 |
| Internal:Known | 0.6552 | 0.7917 | 0.7170 | 19 | 10 | 5 | 24 |
| 3'UTR:Known | 0.7273 | 0.8889 | 0.8000 | 8 | 3 | 1 | 9 |
| **ALL** | **0.7423** | **0.9556** | **0.8355** | 23,391 | 8,121 | 1,087 | 24,478 |

### pred_preddepth

| class | P | R | F1 | TP | FP | FN | n_real |
|---|--:|--:|--:|--:|--:|--:|--:|
| Annotated | 0.9512 | 0.9949 | 0.9725 | 17,407 | 894 | 90 | 17,497 |
| 5'UTR | 0.3625 | 0.7152 | 0.4812 | 1,821 | 3,202 | 725 | 2,546 |
| Truncated | 0.9363 | 0.9941 | 0.9643 | 1,866 | 127 | 11 | 1,877 |
| Novel | 0.2913 | 0.7587 | 0.4210 | 478 | 1,163 | 152 | 630 |
| 5'UTR:Known | 0.5533 | 0.8131 | 0.6585 | 322 | 260 | 74 | 396 |
| Extended | 0.9634 | 0.9946 | 0.9788 | 369 | 14 | 2 | 371 |
| Internal:CDSFrameOverlap | 0.8860 | 0.8860 | 0.8860 | 303 | 39 | 39 | 342 |
| 3'UTR:CDSFrameOverlap | 0.9303 | 0.9419 | 0.9361 | 227 | 17 | 14 | 241 |
| Internal | 0.1934 | 0.3354 | 0.2454 | 53 | 221 | 105 | 158 |
| Truncated:Known | 0.9430 | 0.9933 | 0.9675 | 149 | 9 | 1 | 150 |
| 3'UTR | 0.0616 | 0.3578 | 0.1051 | 39 | 594 | 70 | 109 |
| Extended:Known | 0.9839 | 0.9839 | 0.9839 | 61 | 1 | 1 | 62 |
| 5'UTR:CDSFrameOverlap | 0.4762 | 0.7317 | 0.5769 | 30 | 33 | 11 | 41 |
| Extended:CDSFrameOverlap | 0.9615 | 1.0000 | 0.9804 | 25 | 1 | 0 | 25 |
| Internal:Known | 0.7600 | 0.7917 | 0.7755 | 19 | 6 | 5 | 24 |
| 3'UTR:Known | 0.7273 | 0.8889 | 0.8000 | 8 | 3 | 1 | 9 |
| **ALL** | **0.7787** | **0.9469** | **0.8546** | 23,177 | 6,587 | 1,301 | 24,478 |

## What this establishes

**The two callers agree on the shape of the failure, which is the point of running the second
one.** Annotated calls are essentially flat between observed and predicted (RiboCode +1.8%,
Ribo-TISH +2.3%), so the model reproduces canonical CDS faithfully. Every non-canonical class
inflates, and the two callers rank the inflation identically: 3'UTR / dORF worst (RiboCode 6.2x at
obsdepth, Ribo-TISH 5.5x), then novel (2.5x and 2.2x), then uORF / 5'UTR (2.0x and 1.8x).
`internal` is the sole class that falls in both (0.59x and 0.95x).

**3'UTR is the worst class by a wide margin under Ribo-TISH**, precision 0.049 at obsdepth against
1,011 false positives on 109 real calls, roughly ten false 3'UTR ORFs per true one. This is the
same over-smoothing signature the project already records for RiboCode: the softmax leaks diffuse
probability into 3' UTRs that the observed data does not support.

**The depth dial behaves as it does elsewhere.** `pred_preddepth` cuts 3'UTR false positives from
1,011 to 594 and lifts overall precision 0.742 to 0.779, costing recall 0.956 to 0.947.

## Limits

Ribo-TISH and RiboCode are not yet on a common threshold: Ribo-TISH gates on `--fpth 0.05` frame p,
RiboCode on `--pval 0.05` plus the loader's 90 nt and 0.5x enrichment gates. No precision or recall
statement should be compared across the two callers until that is imposed.

The class vocabularies also differ. `Truncated` and `Extended` have no RiboCode equivalent, and
Ribo-TISH's `Novel` is not the same set as RiboCode's `novel`.
