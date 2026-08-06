# FIGURE o2_venn -- Wang observed / Janich observed / model predicted, three-way ORF-call overlap

**What it shows:** the Wang-Janich overlap IS the between-experiment reproducibility ceiling, made
visual. What a pairwise F1 table hides is the region where the model recovers ORFs only ONE real
experiment saw -- those are not false positives, they are the model siding with one experiment where
the other lacked depth. That is the "as informative as a second experiment" claim made concrete.

Layout: rows = transcript space (INTERSECTION of TPM>=1, then UNION), columns = ORF class
(CDS, uORF, novel). One figure per (model arm, operating point).

## Panels / files

| file | model arm | operating point |
|---|---|---|
| `o2_venn_wang_theta1.png` / `.pdf` | Wang-trained (`wang_nokozak`) | theta = 1, uncalibrated |
| `o2_venn_wang_calib.png` / `.pdf` | Wang-trained | CDS-anchored calibrated theta |
| `o2_venn_janich_theta1.png` / `.pdf` | Janich-universe (`janich_wanguni`) | theta = 1, uncalibrated |
| `o2_venn_janich_calib.png` / `.pdf` | Janich-universe | CDS-anchored calibrated theta |

Two arms per the standing two-arm rule (`feedback_two_arm_orf_calling`): an uncalibrated panel
inflates the model-only lobe with calibration slack rather than biology, so both operating points are
always shown and neither is presented alone.

## Data inputs

Observed calls (RiboCode collapsed, real Ribo-seq):
- `data/heldout_psites/mouse_wang_liver/mouse_wang_liver_collapsed.txt` -- Wang mouse liver
- `data/ribocode_mouse_liver/Liver_5samp/Liver_5samp_collapsed.txt` -- Janich mouse liver, 5 samples

Predicted calls (RiboCode drop-in on the model's predicted profile, predicted depth):
- `results/o2_liver/wang_nokozak/dropin/pred_preddepth_collapsed.txt`
- `results/o2_liver/janich_wanguni/dropin/pred_preddepth_collapsed.txt`

## Regeneration

```
# cas12a env
python figures/o2_venn/make_o2_venn.py
```

## Conventions this figure obeys

- **Transcript space is stated per panel** (`feedback_orf_call_transcript_space`): every panel prints
  its own n and all three call sets are restricted to the same space. The INTERSECTION is the primary
  comparison -- all three parties had a fair shot at every transcript. The UNION is shown to expose the
  coverage asymmetry: a transcript expressed only in Wang cannot be called by Janich, so its absence
  there is not evidence of a missed call.
- Circles are **area-proportional** (matplotlib_venn 1.1.2 solves radii and offsets from the set sizes),
  so region area tracks count rather than being schematic.

## Caveat

The Janich arm uses the REBUILT Janich universe (10,431 tx). The earlier 2,434-tx universe was an
artifact of running salmon on raw untrimmed ARTseq fastqs (95% adapter-bearing, 0.03% mapping); stale
artifacts are quarantined under `data/_quarantine/janich_stale_universe_2026-07-30/`. Any panel
regenerated against the old universe is not comparable to these.
