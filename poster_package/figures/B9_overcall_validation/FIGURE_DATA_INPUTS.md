# B9 -- over-call validation: are the model's unsupported ORF calls false positives?

**What it shows:** the measurement behind P9's "model only" annotation. P9 carries the headline;
the four strata and the gene-coverage split that produced it had no figure until 2026-08-14. This is
a **rigor** panel -- its subject is a failure mode.

## Panel -> data

| element | source |
|---|---|
| all three panels | `results/merged_liver_ribocode/overcall/{overcall_meta.json,overcall_recovery_by_class.tsv}` |
| produced by | `scripts/score_merged_liver_overcall.py` |
| the arbiter | `results/merged_liver_ribocode/merged_liver_collapsed.txt` -- joint RiboCode over **28** aligned with the final recipe mouse-liver libraries (Janich 5 + GSE243134 21 + Wang 2), 26,388 raw calls |
| model arm | `pred_preddepth` (standalone), GSE243134 RNA -- matching P9 |
| transcript space | 22,974 tx / 12,374 genes; merged reference = 15,254 calls in it |

Filtering and genomic keying via `compare_dropin_calls.build_loader`, the single path behind every
drop-in number in this project.

## Panel A is the whole argument, and it needs all three bars

A recovery rate for the model's extras means nothing alone: a deeper reference recovers more of
everything. The two reference strata are what make 7.7% interpretable.

| stratum | n | corroborated | what it establishes |
|---|--:|--:|---|
| **ANCHOR** model calls that ARE experiment-supported | 11,441 | **99.3%** | the merged reference is not the limitation |
| **CONTROL** observed calls unique to ONE experiment | 1,054 | **74.5%** | how often "unsupported by the other arms" just means "too deep for them" |
| **THE QUESTION** model calls no experiment supports | **1,942** | **7.7%** | a ~10x shortfall -- false positives |

**Never plot or quote the 7.7% without the 74.5%.** The generator builds them into one panel so they
cannot be separated.

## Panel B: cross-architecture agreement is NOT evidence

| | n | corroborated |
|---|--:|--:|
| mamba4 alone | 1,942 | 7.7% |
| attn alone | 2,056 | 8.6% |
| **BOTH models, no experiment** | 1,504 | **8.7%** |

The 5-set UpSet family was built expecting the shared bar to identify real-but-shallow ORFs. It does
not. **Never use cross-architecture consensus as a confidence filter** -- it is a shared systematic
error. The falsified hypothesis is recorded beside `FAMILIES` in the P9 generator.

## Panel C: per class, and the gene-coverage split (in the suptitle)

The gap holds in every class, including canonical CDS (control 92-100% vs model-unique 8.5%).
`internal` is the one class where the CONTROL is itself weak (33-54%), so it is the one place a low
model number is genuinely ambiguous -- stated on the panel.

Gene-coverage split, which closes the "no data there" objection:

| | n | corroborated |
|---|--:|--:|
| extras on a gene the merged reference DOES call | 881 | **16.9%** |
| extras on a gene with NO merged call anywhere | **1,061** | 0% **by construction** |

**The second row's rate is definitional** -- recovery requires a merged call on the gene, and the
stratum is defined by their absence. Read its `n`, never its rate; the figure shows it only in the
suptitle for that reason. Since the universe is already salmon TPM >= 1, those 1,061 are
RNA-expressed transcripts with no ribosome signal in 28 pooled libraries: the lncRNA over-call mode.

## Caveats a caption must carry

- This is ONE tissue (mouse liver) and one model family. The same construction has not been done for
  any human dataset.
- The merged reference is a fourth, deeper validation arm; it does NOT replace the three independent
  observed arms the 3x3 factorial needs.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
$PY scripts/score_merged_liver_overcall.py
cd figures/B9_overcall_validation && $PY make_overcall_validation.py
```
