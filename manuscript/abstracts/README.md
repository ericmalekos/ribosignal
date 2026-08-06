# Abstracts

Moved here 2026-08-05 from `/data/tmp/emalekos/` (node-local scratch on mustard, periodically
cleaned and invisible from compute nodes -- wrong home for a deliverable). Copies verified
byte-identical before the originals were removed.

| file | what |
|---|---|
| `MALEKOS_Eric_AIandBioengineering_UCSC.md` | BIC 2026 abstract, written 2026-08-01, word-count constrained |
| `BIC2026_portal_fields.md` | the portal submission fields for the same |

## STALENESS: the macrophage numbers are from superseded searches

Written 2026-08-01. The MS numbers quoted come from `pgx_{attn,mamba4}_union`, which ran with
`calibrate_mass = 2` and `fragment_mass_tolerance = 0.6 Da` -- verified from the retained
`fragger.params`. On 2026-08-03 that configuration was shown to re-derive six search parameters PER
DATABASE, so the arms of a comparison were not scored under the same rules (see
proteogenomics/methods.md, "Search parameters are FROZEN").

Status of each quoted figure:

| abstract text | status |
|---|---|
| "1,244 sequence proteomic search database" | **VALID** -- database size is model-derived, not search-derived |
| "identified 33 novel peptides" | **stale** -- search-dependent |
| "costing only 390 canonical spectrum matches" | **stale, most affected** -- dPSM is the column the artifact distorted |
| "destroys 4,634 canonical peptides" | **stale, most affected** -- same column, null_nc arm |
| held-out hepatocyte F1 0.90 / uORF precision 0.766 / non-canonical 0.591 | **VALID** -- model-eval metrics, no MS search involved |

Scale of the correction, measured on Microglia: frozen vs unfrozen moved canonical IDs
205,348 -> 214,214 (+4.3%) and dPSM +11,017 -> +65. The direction for the abstract's cost figures is
therefore "probably smaller than stated", but they must be re-measured, not adjusted.

**To refresh:** re-run the BMDM five-arm, both-models comparison under the frozen template. That is
task #50 ("Rebuild macrophage churn on union attn + mamba4, two threshold arms"), still pending. The
2026-08-03 cross-subtype re-run does NOT cover it -- that used the RiboCode-called per-population
design (gencode / model_predicted / riboNT / riboALL / null_atg / null_nc), not
mamba4-vs-attn x standard-vs-poisson.

## Model description

The abstract describes the TRANSFORMER (~5M params), per an explicit instruction on 2026-08-01 to
focus on it and set Mamba aside. Note the tension: the 2026-07-31 decision makes `orf_v2_mamba4`
the primary/headline model, so an abstract describing the transformer describes what the paper
positions as supplemental.

This is defensible specifically for the MS claim -- the 2026-08-05 four-dataset comparison found the
two models indistinguishable for proteogenomic discovery (49 vs 48 novel peptides; sign of the
density difference flips by dataset). It would NOT be defensible for a profile-prediction claim,
where mamba4 genuinely wins (+0.0204 held-out Pearson, non-overlapping seed ranges).
