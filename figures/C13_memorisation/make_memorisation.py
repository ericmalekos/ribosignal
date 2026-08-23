#!/usr/bin/env python3
"""Figure C13: how much of held-out performance is memorisation, measured at matched depth.

LOTO holds out a TISSUE, not transcripts. ~99.7% of the held-out tissue's transcripts also appear in
training, in a different cell type, so "held-out" performance is mostly performance on molecules the
model has already taken gradients on. This figure quantifies that instead of asserting it is small.

The honest comparison is the whole point, and it is not the obvious one:

  all seen         r = 0.616   the naive number, and it is badly inflated
  depth-matched    r = 0.248   seen transcripts matched to the unseen ones on log10 P-site count
                               and length (medians 67 vs 66 counts)
  never seen       r = 0.148   the 98 transcripts absent from every training tissue

Comparing 0.616 against 0.148 would attribute ~0.47 to memorisation. Almost all of that gap is
DEPTH: the unseen transcripts are shallow, and profile correlation rises steeply with counts. Once
depth is matched, the memorisation effect is ~0.10 with disjoint 95% CIs -- real, measurable, and
about a quarter of the naive estimate.

The figure is drawn with the inflated bar included and greyed rather than dropped, because the
reader's default mental comparison is the naive one and the figure has to intercept it.

Reconciles with the rest of the paper: profile-correlation metrics are memorisation-sensitive, while
thresholded ORF calling largely is not (cross-species drop-in F1 holds at 0.929-0.931). That is why
the ORF-call results carry the generalization claim and the profile-Pearson numbers do not.

Values are read from results/memorisation_control/hepatocytes_{attn,mamba4}.json at build time.

cas12a env. Regenerate: `python3 make_memorisation.py`.
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
HERE = NEW / "figures/C13_memorisation"
SRCDIR = NEW / "results/memorisation_control"
MODELS = [("attn", "attn (transformer)"), ("mamba4", "mamba4")]

# (json key for the mean, json key for the CI or None, label, colour, is_the_misleading_one)
TIERS = [
    ("mean_r_all_seen", None, "all seen\n(unmatched depth)", "#cbd5e0", True),
    ("mean_r_matched_seen", "ci_matched_seen", "seen,\ndepth-matched", "#2b6cb0", False),
    ("mean_r_unseen", "ci_unseen", "never seen", "#c05621", False),
]


def main():
    data = {}
    for m, _ in MODELS:
        p = SRCDIR / f"hepatocytes_{m}.json"
        if not p.exists():
            raise SystemExit(f"missing {p} -- run scripts/memorisation_control.py first")
        data[m] = json.loads(p.read_text())

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.7), constrained_layout=True, sharey=True)
    for ax, (m, title) in zip(axes, MODELS):
        J = data[m]
        for i, (key, cikey, lab, col, faded) in enumerate(TIERS):
            v = J[key]
            ax.bar(i, v, width=0.62, color=col, edgecolor="none",
                   alpha=0.55 if faded else 1.0, zorder=3)
            if cikey:
                lo, hi = J[cikey]
                ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], color="black", lw=1.2,
                            capsize=3, zorder=4)
            ax.text(i, v + 0.022, f"{v:.3f}", ha="center", fontsize=8.5, fontweight="bold")
        n = (J["n_all_seen"], J["n_matched_seen"], J["n_unseen"])
        for i, nn in enumerate(n):
            ax.text(i, -0.045, f"n={nn:,}", ha="center", fontsize=7, color="0.35")

        # The measured effect, drawn between the two comparable bars only.
        d = J["mean_r_matched_seen"] - J["mean_r_unseen"]
        top = max(J["mean_r_matched_seen"], J["mean_r_unseen"]) + 0.085
        ax.plot([1, 1, 2, 2], [top - 0.012, top, top, top - 0.012], color="black", lw=0.9)
        ax.text(1.5, top + 0.012, f"memorisation\n$\\Delta$ = {d:.3f}", ha="center",
                fontsize=7.5, fontweight="bold")

        ax.set_xticks(range(len(TIERS)))
        ax.set_xticklabels([t[2] for t in TIERS], fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(-0.07, 0.78)
        ax.grid(axis="y", alpha=0.25, zorder=0)
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("per-transcript profile Pearson r\n(held-out Hepatocytes)", fontsize=8.5)

    fig.suptitle("Matching on depth shrinks the memorisation effect from ~0.47 to ~0.10",
                 fontsize=10)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "C13_memorisation.pdf")
    fig.savefig(HERE / "C13_memorisation.png", dpi=200)
    plt.close(fig)

    out = {"source_dir": str(SRCDIR), "holdout": data["attn"]["holdout"],
           "models": {m: {k: data[m][k] for k in data[m] if k != "npz"} for m, _ in MODELS},
           "naive_delta": {m: data[m]["mean_r_all_seen"] - data[m]["mean_r_unseen"] for m, _ in MODELS},
           "matched_delta": {m: data[m]["mean_r_matched_seen"] - data[m]["mean_r_unseen"]
                             for m, _ in MODELS}}
    (HERE / "C13_values.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/'C13_memorisation.pdf'}")
    for m, _ in MODELS:
        J = data[m]
        print(f"  {m}: naive delta {J['mean_r_all_seen']-J['mean_r_unseen']:.3f} -> "
              f"depth-matched {J['mean_r_matched_seen']-J['mean_r_unseen']:.3f} "
              f"(CIs {J['ci_unseen'][0]:.3f}-{J['ci_unseen'][1]:.3f} vs "
              f"{J['ci_matched_seen'][0]:.3f}-{J['ci_matched_seen'][1]:.3f})")


if __name__ == "__main__":
    main()
