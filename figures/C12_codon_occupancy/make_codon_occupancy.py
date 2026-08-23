#!/usr/bin/env python3
"""Figure C12: the model recovers codon-level occupancy it was never given.

Nothing in training supplies codon identity or dwell time. The model sees sequence, an ORF track and
RNA-seq coverage, and is trained on per-nucleotide P-site counts. If its predicted profile
nonetheless reproduces the measured per-codon occupancy, it has learned something about elongation
rather than only about where ORFs are.

The whole figure turns on plotting the right ceiling. Codon dwell is protocol-sensitive
(cycloheximide, digestion, depth), so a raw correlation of 0.26 is uninterpretable on its own. Two
empirical ceilings are measured first, from observed-vs-observed pairs with no model involved:

  within species   r = 0.944   codon occupancy is highly reproducible between independent experiments
  across species   r = 0.269   and is almost entirely species-specific, consistent with tRNA pools

Against those, the two model numbers say different things and must not be read as a ranking:

  human   0.548  = 58% of the within-species ceiling. Substantial, but these transcripts were ~99.7%
                   present in training, so it carries the memorisation component of figure C13.
  mouse   0.263  = 97% of the ACROSS-species ceiling. The transcripts were never seen, and a
                   human-trained model predicting mouse codon dwell cannot exceed how much codon
                   structure transfers between the species at all. It is at the transfer limit.

Same-dataset pairs are excluded from the ceiling: they share one observed vector and correlate at
1.000 by construction, which would inflate it.

The figure is drawn as points against horizontal ceiling lines rather than as bars, because the
claim is "fraction of an achievable maximum", and a bar chart invites reading 0.548 > 0.263 as the
model doing better on human.

Values are read from results/codon_occupancy_canon/_compare_asite.json at build time.

cas12a env. Regenerate: `python3 make_codon_occupancy.py`.
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
HERE = NEW / "figures/C12_codon_occupancy"
SRC = NEW / "results/codon_occupancy_canon/_compare_asite.json"


def main():
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} -- run scripts/codon_compare.py first")
    J = json.loads(SRC.read_text())

    ceil_w = J["ceiling_within_species_median"]
    ceil_a = J["ceiling_across_species_median"]
    within = [c["r"] for c in J["ceiling"] if c["same_species"]]
    across = [c["r"] for c in J["ceiling"] if not c["same_species"]]
    mouse = [m["r"] for m in J["model_vs_observed"] if m["species"] == "mouse"]
    human = [m["r"] for m in J["model_vs_observed"] if m["species"] == "human"]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.4, 3.6), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.0, 1.15]})

    # --- left: the empirical ceilings, model-free ---
    rng = np.random.default_rng(0)   # jitter only; seeded so the figure is reproducible
    for i, (vals, lab, col) in enumerate([(within, f"within species\nn={len(within)}", "#2b6cb0"),
                                          (across, f"across species\nn={len(across)}", "#a0aec0")]):
        x = i + (rng.random(len(vals)) - 0.5) * 0.22
        axL.scatter(x, vals, s=14, alpha=0.65, color=col, edgecolor="none", zorder=3)
        med = float(np.median(vals))
        axL.plot([i - 0.28, i + 0.28], [med, med], color="black", lw=2, zorder=4)
        axL.text(i, med + 0.045, f"{med:.3f}", ha="center", fontsize=8.5, fontweight="bold")
    axL.set_xticks([0, 1])
    axL.set_xticklabels(["within species", "across species"], fontsize=8.5)
    axL.set_ylabel("per-codon occupancy r\n(observed vs observed)", fontsize=8.5)
    axL.set_title("Empirical ceiling: no model involved", fontsize=9)
    axL.set_ylim(-0.1, 1.0)
    axL.axhline(0, color="0.8", lw=0.8, zorder=1)
    axL.grid(axis="y", alpha=0.25, zorder=0)
    axL.tick_params(labelsize=8)

    # --- right: the model against the ceiling that applies to it ---
    groups = [("mouse\n(never seen)", mouse, ceil_a, "across-species ceiling", "#c05621"),
              ("human\n(~99.7% seen)", human, ceil_w, "within-species ceiling", "#2b6cb0")]
    for i, (lab, vals, ceil, ceil_lab, col) in enumerate(groups):
        x = i + (rng.random(len(vals)) - 0.5) * 0.22
        axR.scatter(x, vals, s=26, alpha=0.85, color=col, edgecolor="white", lw=0.4, zorder=3)
        med = float(np.median(vals))
        axR.plot([i - 0.3, i + 0.3], [med, med], color="black", lw=2, zorder=4)
        axR.plot([i - 0.36, i + 0.36], [ceil, ceil], color=col, lw=1.6, ls="--", zorder=2)
        # Ceiling label ABOVE its own line, not to the right: on the mouse group the median sits
        # almost on the ceiling, so a right-hand label collided with the value annotation.
        axR.text(i - 0.36, ceil + 0.022, ceil_lab, va="bottom", ha="left", fontsize=7, color=col)
        # Label to the RIGHT of the swarm, not under the median: at these medians the text sat on
        # top of the points it describes.
        axR.annotate(f"{med:.3f}\n{100*med/ceil:.0f}% of ceiling",
                     xy=(i + 0.30, med), xytext=(i + 0.42, med), va="center", ha="left",
                     fontsize=8, fontweight="bold", color=col)
    axR.set_xticks([0, 1])
    axR.set_xticklabels([g[0] for g in groups], fontsize=8.5)
    axR.set_xlim(-0.6, 1.9)
    axR.set_ylim(-0.05, 1.0)
    axR.set_ylabel("per-codon occupancy r\n(predicted vs observed)", fontsize=8.5)
    axR.set_title("Model vs the ceiling that applies to it", fontsize=9)
    axR.grid(axis="y", alpha=0.25, zorder=0)
    axR.tick_params(labelsize=8)

    fig.suptitle("Codon occupancy: learned from sequence, never supplied in training", fontsize=10)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "C12_codon_occupancy.pdf")
    fig.savefig(HERE / "C12_codon_occupancy.png", dpi=200)
    plt.close(fig)

    out = {"source_json": str(SRC), "site": J["site"], "n_dumps": J["n_dumps"],
           "ceiling_within_species_median": ceil_w, "ceiling_across_species_median": ceil_a,
           "n_pairs_within": len(within), "n_pairs_across": len(across),
           "model_mouse_median": float(np.median(mouse)), "model_human_median": float(np.median(human)),
           "mouse_frac_of_across_ceiling": float(np.median(mouse) / ceil_a),
           "human_frac_of_within_ceiling": float(np.median(human) / ceil_w),
           "model_points": J["model_vs_observed"]}
    (HERE / "C12_values.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/'C12_codon_occupancy.pdf'}")
    print(f"  ceiling within {ceil_w:.3f} (n={len(within)}) / across {ceil_a:.3f} (n={len(across)})")
    print(f"  mouse {np.median(mouse):.3f} = {100*np.median(mouse)/ceil_a:.0f}% of across-species ceiling")
    print(f"  human {np.median(human):.3f} = {100*np.median(human)/ceil_w:.0f}% of within-species ceiling")


if __name__ == "__main__":
    main()
