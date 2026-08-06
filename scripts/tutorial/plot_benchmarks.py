#!/usr/bin/env python3
"""Living benchmark figure for the tutorial, built from results/orf_call_metrics/*.json.

Reads whatever per-model JSONs exist (written by scripts/heldout/orf_call_metrics.py) and plots the
standalone theta=1 ORF-call quality per model -- non-canonical F1 and novel-ORF precision -- grouped by
evaluation universe. Regenerated on every `make html`, so the figure always tracks the latest evals.
Run with a python that has matplotlib (e.g. conda_envs/cas12a); no-ops cleanly if nothing is present.
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
JDIR = os.path.join(REPO, "results", "orf_call_metrics")
OUT = os.path.join(REPO, "figures", "benchmarks")

# universe -> (display, colour); models sorted within universe by non-canon F1
UNI = {
    "fibroblast": ("Fibroblast universe (in-distribution)", "#3a6ea5"),
    "union": ("Union universe (broad)", "#c1666b"),
    "mouse-tcell": ("Mouse T-cell (cross-species held-out)", "#5b8c5a"),
}


def short(label):
    return (label.replace(" (baseline)", "").replace(" onehot", "").replace("onehot ", "")
            .replace(" (mouse held-out)", "").replace("GSE155087", "").strip())


def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(JDIR, "*.json"))):
        try:
            d = json.load(open(f))
        except (OSError, ValueError):
            continue
        p = d.get("preddepth_theta1_vs_real")
        if not p:
            continue
        pz = (d.get("poisson_cds_anchored") or {}).get("score", {})
        rows.append({
            "label": short(d.get("label", os.path.basename(f))),
            "uni": d.get("universe", "?"),
            "nonc_f1": p["non-canonical"]["f1"],
            "novel_p": p["novel"]["precision"],
            "pois_novel_p": (pz.get("novel") or {}).get("precision"),
        })
    # de-dup by (label, uni), keep first
    seen, out = set(), []
    for r in rows:
        k = (r["label"], r["uni"])
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out


def main():
    rows = load()
    if not rows:
        print("plot_benchmarks: no orf_call_metrics JSONs yet; skipping")
        return
    order = []
    for uni in UNI:
        grp = sorted([r for r in rows if r["uni"] == uni], key=lambda r: -r["nonc_f1"])
        order += grp
    order += [r for r in rows if r["uni"] not in UNI]

    labels = [r["label"] for r in order]
    colours = [UNI.get(r["uni"], ("", "#888"))[1] for r in order]
    y = list(range(len(order)))[::-1]

    fig, ax = plt.subplots(1, 2, figsize=(11, 0.5 * len(order) + 1.6), sharey=True)
    ax[0].barh(y, [r["nonc_f1"] for r in order], color=colours, height=0.7)
    ax[0].set_title("Non-canonical ORF-call F1\n(standalone $\\theta$=1)", fontsize=11)
    ax[0].set_xlim(0, 0.75)
    for yi, r in zip(y, order):
        ax[0].text(r["nonc_f1"] + 0.008, yi, f"{r['nonc_f1']:.3f}", va="center", fontsize=8)

    ax[1].barh(y, [r["novel_p"] for r in order], color=colours, height=0.7, alpha=0.85)
    for yi, r in zip(y, order):
        if r["pois_novel_p"] is not None:
            ax[1].plot(r["pois_novel_p"], yi, "D", color="#222", ms=5,
                       label="Poisson arm" if yi == y[0] else None)
    ax[1].set_title("Novel-ORF precision\n(bar = $\\theta$=1, diamond = Poisson)", fontsize=11)
    ax[1].set_xlim(0, 0.9)
    ax[1].legend(loc="lower right", fontsize=8, frameon=False)

    ax[0].set_yticks(y)
    ax[0].set_yticklabels(labels, fontsize=9)
    for a in ax:
        a.grid(axis="x", ls=":", alpha=0.4)
        a.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in UNI.values()]
    fig.legend(handles, [d for d, _ in UNI.values()], loc="upper center",
               ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    os.makedirs(OUT, exist_ok=True)
    dst = os.path.join(OUT, "benchmark_comparison.png")
    fig.savefig(dst, dpi=150, bbox_inches="tight")
    print(f"plot_benchmarks: wrote {dst} ({len(order)} models)")


if __name__ == "__main__":
    main()
