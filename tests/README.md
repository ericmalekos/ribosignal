# Tests

```
conda_envs/cas12a/bin/python3 tests/test_regressions.py
```

About a second. No GPU, no SLURM, no real data.

## What they check

Each test covers a bug that has already happened here. In every case the pipeline ran to
completion and produced plausible numbers, so there was nothing to notice at the time. These tests
make those specific mistakes fail immediately.

| test | what it checks |
|---|---|
| `test_pack_contract` | offsets, lengths and dtypes are consistent, and row `k` of `tx_order` maps to `offsets[k]:offsets[k+1]` in coverage, target and the ORF track alike |
| `test_pack_rejects_missing_tx_instead_of_skipping` | a transcript missing from the one-hot FASTA raises instead of being skipped. Skipping once reduced a cross-species evaluation to a biased subset |
| `test_pack_rejects_length_mismatch` | a transcript whose length disagrees between sources raises. Annotation-version drift would otherwise misalign the ORF track against the labels |
| `test_orf_coordinate_conversion` | RiboCode's 1-based inclusive coordinates convert to 0-based half-open correctly (`start0 = tstart - 1`, `end0 = tstop`) |
| `test_orf_qval_falls_back_to_pval` | a missing `adjusted_pval` is not treated as passing |
| `test_cds_relationship_splits_out_of_frame_overlaps` | out-of-frame overlaps stay a separate class. Collapsing them left 117 of 1,533 BMDM model sequences absent from a 2.3M-sequence null |
| `test_novel_classes_cover_the_overlap_types` | the overlap classes are still mapped to NOVEL, so they do not drop out of the null again |
| `test_tryptic_enzyme_does_not_cut_before_proline` | `search_enzyme_nocut_1` is set. Leaving it blank gives trypsin/P and moved the BMDM baseline from 335,548 to 352,307 PSMs |
| `test_frozen_templates_disable_calibration` | `calibrate_mass = 0` in the frozen templates. At `2`, MSFragger re-derives up to six parameters per database, which added 11,017 PSMs to one macrophage column |
| `test_search_hash_ignores_database_path` | the search cache key ignores the absolute database path. Including it made `--shared-search-root` do nothing |

## Checked by breaking them

Each test was confirmed to fail when the thing it checks is broken:

| change made | test failed? |
|---|---|
| `search_enzyme_nocut_1` blanked | yes |
| `overlap_uorf` / `overlap_dorf` collapsed into `cds_offframe` | yes |
| ORF `start0` shifted by +1 | yes |
| `write_pack` made to skip a missing transcript | yes |

Worth repeating after editing a test, particularly if you loosen an assertion.

## What they do not check

These test code logic on synthetic fixtures. They do not look at real data, so they would not have
caught:

- **The Janich untrimmed FASTQs.** 95% of reads still carried adapter, salmon mapped 0.03%, and the
  universe came out at 2,434 transcripts. Needs a mapping-rate check on real input.
- **The too-wide MSFragger fragment tolerances.** Needs a comparison against real spectra.
- **A stale model checkpoint.** The A549 and immunopeptidome results ran on an old checkpoint for
  weeks with nothing reporting it.

Catching those needs checks on live inputs at pipeline entry: mapping rate, adapter fraction,
checkpoint identity. That does not exist yet. Most of the real damage in this project came from
that category rather than from the logic errors above, so it is worth adding.
