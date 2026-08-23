# P5 -- held-out human validation, four datasets: data inputs

Poster panel, and the poster's central validation claim. A: ORF-call F1 per dataset, both arms.
B: recall split into annotated CDS vs non-canonical.

## Panel -> data

| arm | kind of held-out | source | produced by |
|---|---|---|---|
| hepatocyte | LOTO holdout (tissue withheld, same study) | `results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/dropin_metrics.json` | the deployed LOTO run |
| iPSC-CM | cross-study, out-of-distribution CELL TYPE | `results/heldout/human_ruizorera/released_mamba4/dropin/dropin_metrics.json` | Ruiz-Orera 2024, ENA PRJEB65856 |
| THP-1 | cross-study, monocytic line | `results/human_orf_calls/scored/human_by_class.tsv` | `scripts/score_human_orf_calls.py` |
| CAR-T | cross-study, primary engineered T cells | same TSV (or `results/human_orf_calls_canon/...` when it exists) | same |

Model is **mamba4** throughout (locked decision D1b). Both calling arms are always shown
(`feedback_two_arm_orf_calling`); every bar states `n_ref`
(`feedback_orf_call_transcript_space`).

## The two source shapes are pooled differently

`dropin_metrics.json` already carries `recall_grouped` (annotated / non-canonical). The scored TSV is
per-class, so the generator pools `novel, uORF, internal, Overlap_uORF, dORF, Overlap_dORF` into the
non-canonical group itself. Both end at the same two numbers; the arithmetic lives in
`from_dropin()` / `from_scored()`.

## CAR-T PROVENANCE -- check this before quoting the panel

CAR-T's ORIGINAL pack carries ~29% PCR duplicates: `umi_tools dedup` ran on the genome BAM while the
pack was built from the undeduplicated transcriptome BAM (`docs/PIPELINE_POLICY.md`). The generator
therefore PREFERS `results/human_orf_calls_canon/scored/human_by_class.tsv` when present, prints
which source it used, and records `cart_source_is_canonical` in `P5_values.json`.

**If that flag is false, the CAR-T bars are provisional and must not go on a printed poster.**

## Caveats a caption must carry

- **Panel B is the honest half.** Overall F1 is ~0.9 on every arm because annotated CDS dominates
  each reference set. Split out, the model is at ~0.99 recall on annotated CDS and 0.55-0.60 on
  non-canonical -- and non-canonical is where the proteogenomics application operates.
- **The four arms are not equally hard.** Hepatocyte is the easiest (same study, same protocol, and
  ~99.7% of its transcripts appear in some training tissue -- memorisation is worth ~0.10 profile r
  at matched depth, see `C13_memorisation`). iPSC-CM is the hardest by construction: different lab,
  different protocol, and cardiomyocyte is absent from the training panel entirely.
- GSE39561 is deliberately excluded. It has no observed arm (no 3-nt periodicity, so no P-sites can
  be placed) and it borrows GSE208041's RNA and universe, so its predictions are identical to
  GSE208041's. It is not an independent fifth measurement.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/P5_heldout_human && $PY make_heldout_human_panel.py
```
