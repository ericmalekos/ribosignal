#!/usr/bin/env python3
"""Poster panel P3: what the model was trained on, and what it is asked to predict.

A poster reader's first question is "trained on what?", and the honest answer has three parts the
panel has to carry at once: which cell types, how unevenly they contribute, and what the model
actually sees at inference.

The depth imbalance is the part usually left out and the part that matters. Fibroblast contributes
32 of the 69 TRAINING libraries and HUVEC contributes 3, an 11x spread in P-sites. That is why the
training sampler caps at `--max_tx_per_tissue 6000` rather than pooling everything, and why LOTO
holds out a TISSUE: a random transcript split would let the deepest tissue dominate both sides.

Hepatocytes is drawn in a different colour because it is the LOTO holdout behind every human number
on this poster -- it is here to show its depth is unremarkable, not to suggest it was trained on. Its
5 libraries are counted separately from the 69, since "74 libraries" would read as if all 74 were
trained on.

Brain is NOT plotted. It was dropped from the source study's 9-tissue panel for low periodicity
(period_obs 0.044) before any pack was built, so there is no final-recipe pack to read and nothing to
draw; the 7 + 1 shown here are all that exist.

Panel B states the inputs, because "predicts ribosome density from sequence" is ambiguous about
whether Ribo-seq is used at inference. It is not, in the standalone arm: the model sees one-hot
sequence, a 5-channel ORF track derived from sequence + annotation, and RNA-seq coverage.

All counts are read from the final-recipe packs at build time.

cas12a env. Regenerate: `python3 make_training_data.py`.
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
HERE = NEW / "figures/P3_training_data"
TRAIN = ["Fibroblast", "VSMC", "ES", "Fat", "HA_EC", "HCAEC", "HUVEC"]
HOLDOUT = "Hepatocytes"


def pack_dir(t):
    return NEW / "data/packed_canon" if t == "Fibroblast" else NEW / f"data/packed_canon_{t}"


def main():
    rows = []
    for t in TRAIN + [HOLDOUT]:
        d = pack_dir(t)
        if not d.exists():
            print(f"  MISSING pack for {t}")
            continue
        norm = json.loads((d / "coverage_norm.json").read_text())
        ps = int(np.load(d / "target_counts.npy").sum())
        rows.append(dict(tissue=t, pack=d.name, n_ribo=norm.get("n_ribo_samples", 0),
                         psites=ps, holdout=(t == HOLDOUT)))
    if not rows:
        raise SystemExit("no final-recipe packs found")
    n_tx = sum(1 for _ in open(pack_dir("Fibroblast") / "tx_order.txt"))
    tr = [r for r in rows if not r["holdout"]]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 3.6), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.5, 1.0]})

    srt = sorted(rows, key=lambda r: -r["psites"])
    y = np.arange(len(srt))
    cols = ["#c05621" if r["holdout"] else "#2b6cb0" for r in srt]
    axA.barh(y, [r["psites"] / 1e6 for r in srt], color=cols, height=0.66, zorder=3)
    axA.set_yticks(y)
    axA.set_yticklabels([f"{r['tissue']}  ({r['n_ribo']})" for r in srt], fontsize=8)
    axA.invert_yaxis()
    for i, r in enumerate(srt):
        axA.text(r["psites"] / 1e6 * 1.02, i, f"{r['psites']/1e6:,.0f}M", va="center", fontsize=7.5)
    axA.set_xlabel("P-sites in the training pack (millions)", fontsize=8.5)
    axA.set_title("A  Training composition -- cell type (Ribo-seq libraries)", fontsize=9, loc="left")
    axA.grid(axis="x", alpha=0.25, zorder=0)
    axA.tick_params(labelsize=8)
    axA.set_xlim(0, max(r["psites"] for r in srt) / 1e6 * 1.20)
    hi, lo = max(r["psites"] for r in tr), min(r["psites"] for r in tr)
    axA.text(0.98, 0.04, f"{hi/lo:.0f}x depth spread across training tissues",
             transform=axA.transAxes, ha="right", fontsize=7.5, style="italic", color="0.3")
    axA.scatter([], [], marker="s", s=42, color="#c05621", label=f"{HOLDOUT} -- LOTO holdout")
    axA.scatter([], [], marker="s", s=42, color="#2b6cb0", label=f"training ({len(tr)} cell types)")
    axA.legend(fontsize=7, loc="lower right", bbox_to_anchor=(1.0, 0.10), framealpha=0.9)

    axB.axis("off")
    txt = [
        ("Universe", f"{n_tx:,} transcripts\nprotein-coding + lncRNA, salmon TPM $\\geq$ 1,\n"
                     "$\\leq$ 10,000 nt, chrM excluded"),
        ("Inputs at inference", "one-hot sequence (4 ch)\nORF track (5 ch: frame 0/1/2, start, stop)\n"
                                "RNA-seq coverage (1 ch)"),
        ("Target", "per-nucleotide Ribo-seq P-sites"),
        ("Standalone arm", "NO Ribo-seq at inference --\nthe arm the application uses"),
    ]
    ypos = 0.97
    for head, body in txt:
        axB.text(0.0, ypos, head, fontsize=8.5, fontweight="bold", va="top", transform=axB.transAxes)
        axB.text(0.0, ypos - 0.075, body, fontsize=7.8, va="top", transform=axB.transAxes,
                 color="0.25", linespacing=1.45)
        ypos -= 0.075 + 0.055 * (body.count("\n") + 1) + 0.055
    axB.set_title("B  What the model sees", fontsize=9, loc="left")

    # Count training and holdout libraries SEPARATELY: "74 libraries" includes the held-out tissue
    # and would read as if all 74 were trained on.
    n_tr = sum(r["n_ribo"] for r in tr)
    n_ho = sum(r["n_ribo"] for r in rows if r["holdout"])
    fig.suptitle(f"Trained on {len(tr)} human cell types ({n_tr} Ribo-seq libraries); "
                 f"{HOLDOUT} ({n_ho} libraries) held out entirely", fontsize=9.5)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "P3_training_data.pdf")
    fig.savefig(HERE / "P3_training_data.png", dpi=200)
    plt.close(fig)

    out = {"universe_tx": n_tx, "holdout": HOLDOUT,
           "n_training_tissues": len(tr),
           "n_ribo_libraries_total": sum(r["n_ribo"] for r in rows),
           "depth_spread_train": hi / lo, "tissues": rows}
    (HERE / "P3_values.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/'P3_training_data.pdf'}")
    print(f"  {len(tr)} training tissues + {HOLDOUT} holdout, "
          f"{sum(r['n_ribo'] for r in rows)} Ribo libraries, {n_tx:,} tx, "
          f"{hi/lo:.0f}x depth spread")


if __name__ == "__main__":
    main()
