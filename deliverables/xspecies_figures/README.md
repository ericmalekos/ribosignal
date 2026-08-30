# Cross-species Ribo-seq model evaluation: figures, data and write-up

Self-contained bundle. Everything needed to read the results, inspect the numbers behind
them, or regenerate the figures from scratch is here; nothing points outside this
directory.

Generated 2026-08-30 from the `xspecies-expansion` branch of `ribosignal`.

## Layout

```
docs/METHODS.md        full methods (~4,600 words)
docs/RESULTS.md        results (~1,300 words incl. tables)
figures/Fig1-4         PDF (vector, for publication) + PNG (400 dpi, for preview)
figures/figure_values.json   every plotted value, with its source table per panel
scripts/make_xspecies_figures.py   regenerates all four figures from data/
data/                  the tables the figures and the write-up are built from
```

## Regenerating the figures

```
python3 scripts/make_xspecies_figures.py
```

Paths resolve relative to the script, so this works wherever the bundle is unpacked.
Requires `matplotlib`, `numpy`, `scipy`. Outputs overwrite `figures/`.

Verified: regenerating inside this bundle reproduces the published figures, with 232 of
234 plotted values identical and the remaining two being `NaN` in both (yeast has no uORF
or ncORF class, so those recalls are undefined rather than zero).

## The figures

| figure | content |
|---|---|
| Fig1_design | runs per species; transcript universes; the Drosophila periodicity failure that dropped that arm |
| Fig2_orf_calls | ORF-call precision/recall; recall by ORF class; uORF recall with class size and depth |
| Fig3_profiles | per-nucleotide profile correlation; position-shuffled null; the Spearman tie artifact |
| Fig4_mm25 | multimap posture effect on the call set, and on the model's own score |

## The data

| file | rows | content |
|---|---|---|
| `xspecies_orf_calls.tsv` | 28 | ORF-call precision/recall/F1 and per-type recall, both architectures |
| `xspecies_orf_classes_attn.tsv` | 28 | CDS / uORF / ncORF / other breakdown, transformer |
| `xspecies_orf_classes_mamba4.tsv` | 24 | same, Mamba |
| `xspecies_profile_eval.tsv` | 14 | per-nt profile agreement, shuffled null, zero fraction, both Spearman variants |
| `xspecies_mm1_vs_mm25_calls.tsv` | 7 | multimap posture comparison, with contributing-run counts |
| `xspecies_model_vs_mm25_reference_attn.tsv` | 14 | model scored against both truth sets |
| `xspecies_align_qc.tsv` | 136 | per-run trimming and alignment metrics |
| `xspecies_biotype_census.tsv` | - | per-species annotation composition |
| `xspecies_dataset_registry.tsv` | 144 | one row per sequencing run: accession, study, adapter, postures |
| `species_reference_registry.tsv` | 9 | one row per species: assembly, annotation source, biotype convention |
| `xspecies_universe_summary.tsv` | 7 | transcript universe sizes (summarised from 21 MB of per-transcript tables) |

## Caveats carried in the data

Three things a reader should know before quoting a number, all documented in `docs/METHODS.md`:

1. **`profile_r` is not the pipeline's `pearson_median`.** It applies log1p to the profile
   first; the pipeline's training and evaluation code does not. The two rank transcripts
   at Spearman 0.65 to 0.79 and share only 34 to 42% of their top decile, so they are not
   interchangeable. Which is the better criterion is unsettled.
2. **`profile_rho_ordinal` is invalid on sparse data.** It is the pipeline's shared
   Spearman, which assigns ordinal ranks without averaging ties; on profiles that are 74
   to 95% zeros this produces spurious negative values. `profile_rho` is the corrected
   average-rank statistic and is the one to use. Both are shipped so the difference is
   auditable.
3. **Chimpanzee's mm25 row is not like-for-like** and is flagged `NO_depth_differs`: five
   of its eight runs failed the caller's periodicity gate under that posture, so its two
   call sets rest on different amounts of data.

## Verification performed

- 50 headline numbers in `docs/RESULTS.md` checked programmatically against the bundled
  tables: 50 verified, 0 mismatched.
- Figure regeneration reproduces published values (232/234 identical, 2 NaN-vs-NaN).
- `SHA256SUMS` covers every file in the bundle.

## Not included

Raw FASTQ, BAMs and packs (932 GB) are in the group warm archive and indexed in
`data/ARCHIVE_INDEX.tsv` of the main project tree, with per-path restore commands. The
per-transcript universe tables (21 MB) are summarised here rather than shipped.
