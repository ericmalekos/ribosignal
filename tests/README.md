# Tests

```
conda_envs/cas12a/bin/python3 tests/test_invariants.py      # ~1 s, no GPU, no SLURM, no real data
```

## What these are for

Every test here pins an invariant that has **already broken silently** in this project. The common
failure mode has never been a crash -- the pipeline ran, produced plausible numbers, and the mistake
surfaced days later. These convert that class of failure into a loud one.

| test | pins |
|---|---|
| `test_pack_contract` | offsets/lengths/dtypes and the axis invariant: row `k` of `tx_order` is `offsets[k]:offsets[k+1]` in coverage, target, AND the ORF track |
| `test_pack_rejects_missing_tx_instead_of_skipping` | the dump that silently skipped onehot/FASTA-mismatched transcripts and collapsed a cross-universe eval to a biased subset |
| `test_pack_rejects_length_mismatch` | annotation-version drift that would misalign the ORF track and embeddings against the labels |
| `test_orf_coordinate_conversion` | RiboCode 1-based-inclusive -> pgx 0-based half-open (`start0 = tstart-1`, `end0 = tstop`) |
| `test_orf_qval_falls_back_to_pval` | a missing `adjusted_pval` must not become a default-pass |
| `test_cds_relationship_splits_out_of_frame_overlaps` | the null-escape bug: 117 of 1,533 BMDM model sequences absent from a 2.3M near-cognate null because out-of-frame overlaps were collapsed to one class |
| `test_novel_classes_cover_the_overlap_types` | the same classes must map to NOVEL, or they drop out of the null again |
| `test_tryptic_enzyme_does_not_cut_before_proline` | blank `search_enzyme_nocut_1` -> trypsin/P; moved the BMDM baseline 335,548 -> 352,307 PSMs (+5.0%) |
| `test_frozen_templates_disable_calibration` | `calibrate_mass = 2` re-derives up to 6 params PER DATABASE; inflated the macrophage dPSM column by +11,017 and diverged on the HLA and A549 templates too |
| `test_search_hash_ignores_database_path` | absolute DB path in the cache digest made `--shared-search-root` silently do nothing |

## These were mutation-checked

A test that passes but cannot fail is decorative. Each core invariant was deliberately broken and the
suite confirmed to notice:

| mutation | result |
|---|---|
| `search_enzyme_nocut_1` blanked | caught |
| `overlap_uorf`/`overlap_dorf` collapsed into `cds_offframe` | caught |
| ORF `start0` shifted by +1 | caught |
| `write_pack` made to silently drop a missing tx | caught |

Re-run that check after editing a test, especially if you loosen an assertion.

## What this does NOT cover

Structural regressions only. It would NOT have caught:

- **the Janich untrimmed-FASTQ bug** -- 95% of reads carried adapter and salmon mapped 0.03%, giving
  a 2,434-tx universe. Needs a mapping-rate assertion on real data.
- **the too-wide MSFragger fragment tolerances** -- needed a divergence test against real spectra
  (`calibrate_mass = 2` run per DB, diffing the `New <param>` lines).
- **using a stale model checkpoint** -- the A549 and immunopeptidome results ran on the old
  pre-union, pre-nokozak model for weeks without anything complaining.

Those want assertions at pipeline entry on live inputs (mapping rate, adapter fraction, checkpoint
identity), which is a different and still-missing mechanism. Adding one is worth doing: most of this
project's real damage came from that category, not from the logic errors above.
