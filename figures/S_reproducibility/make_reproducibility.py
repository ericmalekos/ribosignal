#!/usr/bin/env python3
"""Figure S-REPRO: can the published pipeline regenerate the model's own training data?

For a methods claim to mean anything, someone must be able to run the pipeline on the public FASTQs
and get back the arrays the model trained on. This figure measures that directly: all 74 Chothani
Ribo-seq runs re-aligned under the final recipe, re-pooled per tissue, and compared against the
stored `target_counts.npy` the model actually consumed.

The figure exists because the claim was previously made on ONE tissue. HUVEC regenerates at 3.04%,
which sounded reassuring -- and turns out to be the BEST of the eight. Median is 3.86% and the worst
is 11.66%, so the single-tissue number understated typical divergence and understated the worst case
by nearly 4x. Panel A is drawn sorted with HUVEC marked so that point is unmissable.

Panel B exists because the two metrics are partly decoupled and a reader who sees only one will draw
the wrong conclusion. |divergence| is how much signal MOVED; per-nt r is whether it moved coherently.
Fat moved the most signal (11.66%) but coherently (r=0.844). HA_EC moved less (8.01%) but
incoherently (r=0.432) and is the one genuine outlier. Plotting them against each other shows there
is no single summary statistic that would do.

Deliberately NOT drawn: any causal explanation for HA_EC. Four hypotheses were tested and rejected
(uniform P-site shift, defective old pack, changed P-site calling, metaplots offset disagreement --
the last recorded as a prediction and falsified by ES). The figure reports the outlier and stops
there.

Also not drawn: a depth axis. There is no depth relationship -- Hepatocytes (1.16B P-sites) and
HUVEC (96M) both diverge at ~3.1% across a 12x range -- so a depth panel would imply a trend the data
does not support.

Values are read from results/chothani_regeneration/divergence_by_tissue.tsv at build time.

cas12a env. Regenerate: `python3 make_reproducibility.py`.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/S_reproducibility"
SRC = NEW / "results/chothani_regeneration/divergence_by_tissue.tsv"
HILITE = "HUVEC"   # the tissue the original single-tissue claim rested on


def main():
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} -- run scripts/summarise_chothani_divergence.py first")
    rows = list(csv.DictReader(open(SRC), delimiter="\t"))
    if len(rows) < 8:
        print(f"  WARNING: only {len(rows)} tissues in {SRC.name}; figure covers a subset")
    for r in rows:
        r["pct"] = float(r["pct_of_signal"])
        r["r"] = float(r["pearson"])
        r["n"] = int(r["old_psites"])

    med = float(np.median([r["pct"] for r in rows]))
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 3.8), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.15, 1.0]})

    # --- A: divergence per tissue, sorted, HUVEC called out ---
    srt = sorted(rows, key=lambda r: r["pct"])
    y = np.arange(len(srt))
    cols = ["#c05621" if r["tissue"] == HILITE else "#4a5568" for r in srt]
    axA.barh(y, [r["pct"] for r in srt], color=cols, height=0.62, zorder=3)
    axA.set_yticks(y)
    axA.set_yticklabels([r["tissue"] for r in srt], fontsize=8.5)
    for i, r in enumerate(srt):
        axA.text(r["pct"] + 0.16, i, f"{r['pct']:.2f}%", va="center", fontsize=7.5)
    axA.axvline(med, color="#2b6cb0", lw=1.4, ls="--", zorder=4)
    # Sit the median label inside the axes below the top bar, not above it: at the top it collided
    # with the panel title.
    axA.text(med + 0.14, len(srt) - 1.45, f"median\n{med:.2f}%", color="#2b6cb0",
             fontsize=7.5, fontweight="bold", va="top")
    axA.set_xlabel("P-site signal that moved on regeneration\n"
                   "$\\Sigma|new-old|\\,/\\,\\Sigma old$", fontsize=8.5)
    axA.set_xlim(0, max(r["pct"] for r in srt) * 1.22)
    axA.set_title("A  Regeneration divergence, per tissue", fontsize=9, loc="left")
    axA.grid(axis="x", alpha=0.25, zorder=0)
    axA.tick_params(labelsize=8)
    hv = next(r for r in srt if r["tissue"] == HILITE)
    axA.annotate("the single tissue the\nclaim used to rest on",
                 xy=(hv["pct"], srt.index(hv)), xytext=(hv["pct"] + 3.4, srt.index(hv) - 0.15),
                 fontsize=7.5, color="#c05621", va="center",
                 arrowprops=dict(arrowstyle="->", color="#c05621", lw=1.0))

    # --- B: the two metrics are not interchangeable ---
    for r in rows:
        out = r["r"] < 0.6
        axB.scatter(r["pct"], r["r"], s=52, zorder=3,
                    color="#c05621" if out else "#2b6cb0",
                    edgecolor="white", lw=0.6)
    # Deterministic label de-overlap: five tissues cluster in the top-left, so hand-placed offsets
    # would silently go wrong the moment a value changes. Walk points top-down and push each label
    # below the previous one if it would collide, with a leader line when it has been moved.
    MINGAP = 0.030
    last_y = None
    for r in sorted(rows, key=lambda r: -r["r"]):
        ly = r["r"] if last_y is None else min(r["r"], last_y - MINGAP)
        last_y = ly
        moved = abs(ly - r["r"]) > 1e-9
        axB.annotate(r["tissue"], xy=(r["pct"], r["r"]),
                     xytext=(r["pct"] + 0.40, ly), va="center", ha="left", fontsize=7.5,
                     color="#c05621" if r["r"] < 0.6 else "0.25",
                     arrowprops=dict(arrowstyle="-", color="0.65", lw=0.6,
                                     shrinkA=0, shrinkB=3) if moved else None)
    axB.set_xlabel("signal moved (%)", fontsize=8.5)
    axB.set_ylabel("per-nt Pearson r (new vs training)", fontsize=8.5)
    axB.set_title("B  How much moved $\\neq$ whether it moved coherently", fontsize=9, loc="left")
    axB.set_ylim(0.30, 1.05)
    axB.set_xlim(0, max(r["pct"] for r in rows) * 1.30)
    axB.grid(alpha=0.25, zorder=0)
    axB.tick_params(labelsize=8)

    # Computed, never typed: the first draft hardcoded "median 4.0%" while the code plotted 3.86%,
    # and that stale literal was briefly copied into results.md and the README.
    lo, hi = min(r["pct"] for r in rows), max(r["pct"] for r in rows)
    fig.suptitle(f"The published pipeline regenerates the training data to {lo:.1f}-{hi:.1f}% "
                 f"(median {med:.1f}%) across {len(rows)} tissues", fontsize=10)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "S_reproducibility.pdf")
    fig.savefig(HERE / "S_reproducibility.png", dpi=200)
    plt.close(fig)

    out = {"source_tsv": str(SRC), "n_tissues": len(rows), "median_pct_of_signal": med,
           "min_pct": min(r["pct"] for r in rows), "max_pct": max(r["pct"] for r in rows),
           "min_r": min(r["r"] for r in rows), "max_r": max(r["r"] for r in rows),
           "highlighted": HILITE,
           "tissues": [{k: r[k] for k in ("tissue", "pack", "n_bams", "old_psites", "new_psites",
                                          "pct_total_change", "pct_of_signal", "pearson")}
                       for r in rows]}
    (HERE / "S_reproducibility_values.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/'S_reproducibility.pdf'}")
    print(f"  {len(rows)} tissues, median {med:.2f}%, range "
          f"{min(r['pct'] for r in rows):.2f}-{max(r['pct'] for r in rows):.2f}%, "
          f"r {min(r['r'] for r in rows):.3f}-{max(r['r'] for r in rows):.3f}")


if __name__ == "__main__":
    main()
