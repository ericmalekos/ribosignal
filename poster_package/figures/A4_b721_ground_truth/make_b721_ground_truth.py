#!/usr/bin/env python3
"""A4 -- B721.221: predicted vs MEASURED ORF calls, against the assay's own split-half ceiling.

This is the project's only NON-NULL benchmark and it had no figure. Every other comparison scores
the model against a null, another selection rule, or a single-dataset reference. Here both arms sit
in one cell line: Sarkizova RNA-seq drives the prediction, Ouspenskaia Ribo-seq (327 M unique
footprints) gives the measured calls, and the model never sees the Ribo-seq. Genomic keying,
restricted to the model's 11,527-gene space (21,974 measured calls genome-wide -> 16,331 in scope).

THE CEILING IS THE POINT, AND WITHOUT IT PANEL A IS MISLEADING. RiboCode was run independently on
two DISJOINT halves of the same Ribo-seq (split by HLA allele, so each half holds complete
libraries). Those two halves reproduce each other at F1 0.889 -- canonical recall 0.973, but
non-canonical recall only 0.498. Two independent measurements of the SAME cells agree on barely half
of each other's non-canonical calls. So a model non-canonical recall of 0.24 is NOT a quarter of
perfect; it is roughly half of what a replicate EXPERIMENT would achieve. Panel C states every model
number as a fraction of that ceiling, which is the only denominator that means anything.

WHAT THE FIGURE SAYS.
  * The two-arm framework behaves as designed on real measured data: Poisson buys precision
    (0.78 -> 0.89 mamba4, 0.76 -> 0.92 attn) and pays recall. This is the strongest independent
    support for the calibration dial, because the reference is measured translation and not a null.
  * mamba4 and attn are indistinguishable AGAIN (F1 0.665 vs 0.666) -- a third replication after the
    macrophage cross-subtype run and the 4-dataset MS panel.
  * Canonical recall (0.70-0.76) far exceeds non-canonical recall (0.07-0.27). The class the
    proteogenomics application depends on is the one the model recovers worst. Panel B exists so that
    cannot be read off a single pooled F1.

The fair headline is the one panel C supports: from RNA-seq alone, against a reference built from
327 M measured footprints, the model reaches ~75% of the assay's own self-agreement on F1 and 77-78%
on canonical ORFs. The non-canonical gap is real and stays visible.

cas12a env. Regenerate: `python3 make_b721_ground_truth.py`.
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
HERE = NEW / "figures/A4_b721_ground_truth"
DC = NEW / "proteogenomics/data/B721_pilot/dropin_compare"
# Order: theta=1 then Poisson within each model, so the calibration effect reads left-to-right.
ARMS = [("mamba4/theta1", "mamba4\n$\\theta$=1", "#2b6cb0"),
        ("mamba4/poisson", "mamba4\nPoisson", "#1b4f72"),
        ("attn/theta1", "attn\n$\\theta$=1", "#c1934a"),
        ("attn/poisson", "attn\nPoisson", "#8e6b3f")]
CEIL = "#2f855a"


def main():
    f = DC / "b721_dropin_metrics.json"
    if not f.exists():
        raise SystemExit(f"missing {f}")
    J = json.loads(f.read_text())
    A = J["arms"]
    # The ceiling is directional (A-vs-B and B-vs-A swap precision and recall); F1 is identical, so
    # take the A-vs-B orientation, which uses the LARGER half as the reference and therefore the
    # n_ref most comparable to the model's (15,703 vs 16,331).
    cj = json.loads((DC / "ceiling_AvsB" / "b721_dropin_metrics.json").read_text())
    C = next(iter(cj["arms"].values()))
    missing = [k for k, _, _ in ARMS if k not in A]
    if missing:
        print(f"  MISSING arms (absent from the figure, not zeroed): {missing}")
    arms = [(k, lab, col) for k, lab, col in ARMS if k in A]

    x = np.arange(len(arms))
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14.2, 4.4), constrained_layout=True)
    w = 0.26

    # A: precision / recall / F1 against the measured calls
    for i, (met, hatch) in enumerate((("precision", None), ("recall", "//"), ("f1", ".."))):
        v = [A[k][met] for k, _, _ in arms]
        b = axA.bar(x + (i - 1) * w, v, w, color=[c for _, _, c in arms], hatch=hatch,
                    edgecolor="white", linewidth=0.5, zorder=3)
        for rect, val in zip(b, v):
            axA.annotate(f"{val:.2f}", (rect.get_x() + rect.get_width() / 2, val),
                         xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6.2)
    axA.axhline(C["f1"], color=CEIL, lw=1.4, ls="--", zorder=4)
    # BELOW the line: above it the label ran into the mamba4/Poisson precision bar (0.89), which
    # sits within 0.01 of the ceiling.
    axA.text(-0.42, C["f1"] - 0.055, f"split-half ceiling, F1 {C['f1']:.3f}",
             ha="left", fontsize=7.4, color=CEIL, fontweight="bold")
    axA.set_ylabel("vs 16,331 MEASURED RiboCode calls", fontsize=8.5)
    axA.set_title("A  Predicted vs measured (bars: precision / recall$^{//}$ / F1$^{..}$)",
                  fontsize=9, loc="left")
    axA.set_ylim(0, 1.10)
    axA.text(0.02, 0.995, "Poisson buys precision, pays recall\n"
                          "-- confirmed against MEASURED translation",
             transform=axA.transAxes, fontsize=7.2, va="top", style="italic", color="0.3")

    # B: the split that a pooled F1 hides
    for i, (met, lab) in enumerate((("recall_annotated", "annotated CDS"),
                                    ("recall_noncanonical", "non-canonical"))):
        v = [A[k][met] for k, _, _ in arms]
        b = axB.bar(x + (i - 0.5) * (w + 0.06), v, w + 0.06,
                    color=[c for _, _, c in arms], alpha=1.0 if i == 0 else 0.45,
                    hatch=None if i == 0 else "///", edgecolor="white", linewidth=0.5, zorder=3)
        # Proxy handle only -- the real bars are per-arm coloured, so labelling them would put one
        # legend entry per ARM per class instead of one per class.
        axB.bar(0, 0, 0, color="0.6", alpha=1.0 if i == 0 else 0.45,
                hatch=None if i == 0 else "///", label=lab)
        for rect, val in zip(b, v):
            axB.annotate(f"{val:.2f}", (rect.get_x() + rect.get_width() / 2, val),
                         xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6.5)
    axB.axhline(C["recall_annotated"], color=CEIL, lw=1.4, ls="--", zorder=4)
    axB.axhline(C["recall_noncanonical"], color=CEIL, lw=1.4, ls=":", zorder=4)
    axB.text(len(arms) - 0.45, C["recall_annotated"] + 0.012,
             f"ceiling, annotated {C['recall_annotated']:.3f}", ha="right", fontsize=7.2,
             color=CEIL, fontweight="bold")
    axB.text(len(arms) - 0.45, C["recall_noncanonical"] + 0.012,
             f"ceiling, NON-CANONICAL {C['recall_noncanonical']:.3f}", ha="right", fontsize=7.2,
             color=CEIL, fontweight="bold")
    axB.set_ylabel("recall", fontsize=8.5)
    axB.set_title("B  Recall by class, vs the ceiling for each", fontsize=9, loc="left")
    # Lower-right: the only region of this panel with no bars (the Poisson non-canonical bars are
    # ~0.07-0.10). Centre-left ran into the mamba4/Poisson annotated bar.
    axB.text(0.98, 0.36, "two independent Ribo-seq halves agree on\nonly HALF the non-canonical calls",
             transform=axB.transAxes, ha="right", fontsize=7.4, style="italic", color=CEIL)
    axB.set_ylim(0, 1.10)
    axB.legend(fontsize=7, frameon=False, loc="upper left", ncol=2)

    # C: the only denominator that means anything
    mets = [("f1", "overall F1"), ("recall_annotated", "annotated CDS recall"),
            ("recall_noncanonical", "non-canonical recall")]
    for i, (met, lab) in enumerate(mets):
        v = [100 * A[k][met] / C[met] for k, _, _ in arms]
        b = axC.bar(x + (i - 1) * w, v, w, color=[c for _, _, c in arms],
                    hatch=[None, "//", ".."][i], edgecolor="white", linewidth=0.5, zorder=3)
        for rect, val in zip(b, v):
            axC.annotate(f"{val:.0f}%", (rect.get_x() + rect.get_width() / 2, val),
                         xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6.4)
    axC.axhline(100, color=CEIL, lw=1.4, ls="--", zorder=4)
    axC.text(len(arms) - 0.45, 102, "the assay's own self-agreement", ha="right", fontsize=7.4,
             color=CEIL, fontweight="bold")
    axC.set_ylabel("% of the split-half ceiling", fontsize=8.5)
    axC.set_title("C  Against the RIGHT denominator (bars as in A)", fontsize=9, loc="left")
    axC.set_ylim(0, 122)

    for ax in (axA, axB, axC):
        ax.set_xticks(x)
        ax.set_xticklabels([lab for _, lab, _ in arms], fontsize=7.6)
        ax.grid(axis="y", alpha=0.25, zorder=0)
        ax.tick_params(labelsize=7.5)

    best = max(arms, key=lambda a: A[a[0]]["f1"])[0]
    fig.suptitle("B721.221: ORF calls predicted from RNA-seq alone vs 327 M MEASURED Ribo-seq "
                 f"footprints -- {100*A[best]['f1']/C['f1']:.0f}% of the assay's own self-agreement "
                 "on F1, and the non-canonical gap is real but half as large as it looks",
                 fontsize=9.4)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "A4_b721_ground_truth.pdf")
    fig.savefig(HERE / "A4_b721_ground_truth.png", dpi=200)
    plt.close(fig)

    out = {"source": str(f), "ceiling_source": str(DC / "ceiling_AvsB"),
           "n_ref_genomewide": J["n_ref_genomewide"], "n_ref_restricted": J["n_ref_restricted"],
           "n_gene_universe": J["n_gene_universe"], "n_tx_universe": J["n_tx_universe"],
           "ceiling": {k: C[k] for k in ("precision", "recall", "f1", "recall_annotated",
                                         "recall_noncanonical", "n_real")
                       if k in C},
           "ceiling_n_ref": cj["n_ref_restricted"],
           "arms": {k: A[k] for k, _, _ in arms},
           "pct_of_ceiling": {k: {m: 100 * A[k][m] / C[m] for m, _ in mets} for k, _, _ in arms}}
    (HERE / "A4_values.json").write_text(json.dumps(out, indent=2))

    print(f"  wrote {HERE/'A4_b721_ground_truth.pdf'}")
    print(f"  ceiling (split-half): F1 {C['f1']:.3f}  annot {C['recall_annotated']:.3f}  "
          f"non-canon {C['recall_noncanonical']:.3f}  (n_ref {cj['n_ref_restricted']:,})")
    print(f"  {'arm':<16}{'prec':>7}{'rec':>7}{'F1':>7}{'annot':>8}{'nonC':>7}   "
          f"{'%ceil F1':>9}{'%annot':>8}{'%nonC':>7}")
    for k, _, _ in arms:
        d = A[k]
        print(f"  {k:<16}{d['precision']:>7.3f}{d['recall']:>7.3f}{d['f1']:>7.3f}"
              f"{d['recall_annotated']:>8.3f}{d['recall_noncanonical']:>7.3f}   "
              f"{100*d['f1']/C['f1']:>8.0f}%{100*d['recall_annotated']/C['recall_annotated']:>7.0f}%"
              f"{100*d['recall_noncanonical']/C['recall_noncanonical']:>6.0f}%")


if __name__ == "__main__":
    main()
