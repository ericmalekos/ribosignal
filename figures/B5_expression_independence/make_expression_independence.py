#!/usr/bin/env python3
"""Figure B5: the translation signal EXCEEDS an expression-only baseline.

Framing note (important): this model deliberately takes per-nt RNA-seq coverage as input, so it IS
expression-aware -- that is what makes it cell-type-specific (the Fig 2 premise). The rigor question a
reviewer asks is therefore not "is it expression-independent" (it should not be) but "does it add a
genuine TRANSLATION signal on top of expression, or is a translated-ORF call just a high-expression
call?". This figure answers that: the model's predicted IN-FRAME FRACTION (pred_frame0, a periodicity
/ shape score) out-discriminates every MAGNITUDE / expression score by a wide margin, length-controlled,
on held-out Hepatocytes.

Bars, grouped by ORF stratum (all ORFs; non-canonical only), length-controlled AUROC of translated vs
decoy ORF:
  - Predicted in-frame fraction  (pred_frame0)          -- the model's periodicity/shape signal
  - Predicted density            (pred_density)          -- the model's magnitude head (~ expression)
  - Observed Ribo-seq density    (obs_density_ceiling)   -- raw experimental magnitude (expression proxy)
The shape score beats both magnitude scores by ~+0.12 (all) / +0.09 (non-canonical) -> the discrimination
is not explained by expression alone.

Caveat captured in FIGURE_DATA_INPUTS.md: the model's DENSITY head does track RNA expression (by design);
the FRAME head is the expression-orthogonal shape signal. All numbers read from localization_metrics.json,
so this regenerates if the model is retrained (e.g. mm1 RNA-coverage ablation). cas12a env.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/B5_expression_independence"
METR = NEW / "results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/localization_metrics.json"

MODEL_C, MAG_C, OBS_C = "#2C6FBB", "#C8A45B", "#B0B0B0"
STRATA = [("all", "All ORFs"), ("noncanonical", "Non-canonical ORFs")]
SCORES = [
    ("auroc_pred_frame0_lengthctrl", "Predicted in-frame\nfraction (shape)", MODEL_C),
    ("auroc_pred_density_lengthctrl", "Predicted density\n(magnitude)", MAG_C),
    ("auroc_obs_density_ceiling", "Observed Ribo-seq\ndensity (magnitude)", OBS_C),
]


def main():
    disc = json.load(open(METR))["discrimination"]
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ng, nb = len(STRATA), len(SCORES)
    w = 0.24
    x = np.arange(ng)
    for j, (key, lab, c) in enumerate(SCORES):
        vals = [disc[s][key] for s, _ in STRATA]
        off = (j - (nb - 1) / 2) * w
        bars = ax.bar(x + off, vals, w, color=c, edgecolor="white", lw=0.5,
                      label=lab)
        for xi, v in zip(x + off, vals):
            ax.text(xi, v + 0.006, f"{v:.3f}", ha="center", va="bottom", fontsize=7.3,
                    color=c if c != OBS_C else "#555")

    # delta annotation: shape - best magnitude, per stratum
    for xi, (s, _) in zip(x, STRATA):
        shape = disc[s]["auroc_pred_frame0_lengthctrl"]
        best_mag = max(disc[s]["auroc_pred_density_lengthctrl"], disc[s]["auroc_obs_density_ceiling"])
        d = shape - best_mag
        ax.annotate(f"+{d:.3f}\nover expression", xy=(xi - w, shape + 0.02),
                    xytext=(xi - w, shape + 0.075), ha="center", va="bottom", fontsize=7.2,
                    color=MODEL_C, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", lw=0.7, color=MODEL_C))

    ax.axhline(0.5, ls=":", lw=0.8, color="#999")
    ax.text(ng - 0.5, 0.505, "chance", fontsize=6.5, color="#999", va="bottom", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in STRATA], fontsize=9.5)
    ax.set_ylim(0.4, 1.02)
    ax.set_ylabel("Length-controlled AUROC\n(translated vs decoy ORF)", fontsize=9)
    ax.legend(frameon=False, fontsize=7.8, loc="upper right", ncol=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Translation discrimination exceeds an expression-only baseline\n"
                 "(held-out Hepatocytes; periodicity signal is orthogonal to expression magnitude)",
                 fontsize=9.8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"B5_expression_independence.{ext}", dpi=300, bbox_inches="tight")

    out = {s: {key: disc[s][key] for key, _, _ in SCORES} for s, _ in STRATA}
    for s, _ in STRATA:
        out[s]["shape_minus_best_magnitude"] = (
            disc[s]["auroc_pred_frame0_lengthctrl"]
            - max(disc[s]["auroc_pred_density_lengthctrl"], disc[s]["auroc_obs_density_ceiling"]))
    (HERE / "B5_values.json").write_text(json.dumps(out, indent=2) + "\n")
    for s, lab in STRATA:
        print(f"  {lab}: shape {disc[s]['auroc_pred_frame0_lengthctrl']:.3f} vs "
              f"magnitude {disc[s]['auroc_pred_density_lengthctrl']:.3f}/"
              f"{disc[s]['auroc_obs_density_ceiling']:.3f} "
              f"(+{out[s]['shape_minus_best_magnitude']:.3f})")


if __name__ == "__main__":
    main()
