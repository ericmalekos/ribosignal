#!/usr/bin/env python3
"""Figure C11: what each ORF-track channel buys, per ORF class.

The ORF track is the model's only structured prior about where reading frames lie. It has five
channels -- frame-0/1/2 occupancy (ch0/ch1/ch2), start-codon propensity (ch3), and is_stop (ch4).
This figure zeroes each in turn, re-runs prediction, and re-calls ORFs, so the readout is the change
in ORF-calling F1 rather than a change in loss. A sixth arm (`f012`) zeroes all three frame channels
together, because they are partly redundant with each other and removing one alone understates them.

The result is a double dissociation, and half of it is the opposite of what the track was designed
to do:

  * The FRAME channels SUPPRESS internal ORFs. Removing any one of them roughly triples internal F1.
    Those channels encode where the ANNOTATED reading frame is, which is exactly the prior that
    argues against a plausible ORF sitting inside a known CDS in a different frame. The model locates
    internal ORFs (saliency ranks their starts at 74) and is then told by its own input not to call
    them.
  * The START channel (ch3) is worth almost nothing on annotated CDS (delta ~0.000) and is by far the
    most valuable channel for uORFs, whose defining evidence is a start codon upstream of the
    annotated one.

This is why the figure is drawn as a diverging heatmap centred on zero rather than as bars: the sign
IS the finding, and a bar chart of magnitudes would hide it. Red means the model calls that class
BETTER without the channel.

It also settled a planned experiment. In-silico mutagenesis around start codons was scheduled to
probe the internal-ORF failure; removing the start channel entirely moves internal F1 by -0.0002, so
a gentler perturbation could only have confirmed a null. The sweep was cancelled before it ran.

All 12 arms (6 x 2 models) were re-run on canonically-aligned data; sign agreement with the earlier
off-recipe run is 57/60, with the three disagreements all within 0.005 of zero.

Values are read from results/orf_channel_ablation_canon/ablation_scores.tsv at build time, which is
written by scripts/score_channel_ablation.py using the same compare_dropin_calls.build_loader
filtering as every other drop-in number in this project. Nothing here is hardcoded.

cas12a env. Regenerate: `python3 make_channel_ablation.py`.
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
HERE = NEW / "figures/C11_channel_ablation"
SRC = NEW / "results/orf_channel_ablation_canon/ablation_scores.tsv"

# Channel order follows the track's own channel index, not effect size, so the figure can be read
# against build_orf_track.py without a mental remap.
ARMS = [("ch0", "frame-0 occupancy"), ("ch1", "frame-1 occupancy"), ("ch2", "frame-2 occupancy"),
        ("ch3", "start propensity"), ("ch4", "is_stop"), ("f012", "all frame channels")]
# Classes ordered by how much of the calling problem they represent, canonical first.
CLASSES = ["annotated", "uORF", "novel", "internal", "dORF"]
MODELS = [("attn", "attn (transformer)"), ("mamba4", "mamba4")]
DROPIN = "pred_obsdepth"   # the arm with an observed depth reference; preddepth is in the JSON too


def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        raise SystemExit(f"no rows in {path}")
    d = {}
    for r in rows:
        d[(r["model"], r["arm"], r["dropin"], r["orf_class"])] = r
    return d


def main():
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} -- run scripts/score_channel_ablation.py first")
    D = load(SRC)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.5), constrained_layout=True)
    vmax = 0.0
    grids = {}
    for m, _ in MODELS:
        g = np.full((len(ARMS), len(CLASSES)), np.nan)
        for i, (arm, _lab) in enumerate(ARMS):
            for j, c in enumerate(CLASSES):
                r = D.get((m, arm, DROPIN, c))
                if r:
                    g[i, j] = float(r["delta"])
        grids[m] = g
        vmax = max(vmax, np.nanmax(np.abs(g)))
    # Symmetric limits so the colour scale means the same thing on both panels and zero is white.
    vmax = float(np.ceil(vmax * 100) / 100)

    for ax, (m, title) in zip(axes, MODELS):
        g = grids[m]
        im = ax.imshow(g, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(CLASSES)))
        ax.set_xticklabels(CLASSES, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(ARMS)))
        ax.set_yticklabels([f"{a}  {lab}" for a, lab in ARMS], fontsize=8)
        ax.set_title(title, fontsize=9)
        for i in range(len(ARMS)):
            for j in range(len(CLASSES)):
                if np.isnan(g[i, j]):
                    continue
                v = g[i, j]
                ax.text(j, i, f"{v:+.3f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(v) > 0.6 * vmax else "black")
        ax.set_xlabel("ORF class", fontsize=8)
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label("$\\Delta$ F1 when the channel is REMOVED\n(red = model calls the class better without it)",
                 fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle("ORF-track channel ablation: the frame channels suppress internal ORFs",
                 fontsize=10)

    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "C11_channel_ablation.pdf")
    fig.savefig(HERE / "C11_channel_ablation.png", dpi=200)
    plt.close(fig)

    out = {"source_tsv": str(SRC), "dropin_arm": DROPIN, "classes": CLASSES,
           "arms": [a for a, _ in ARMS], "models": [m for m, _ in MODELS],
           "delta_f1": {m: {a: {c: (None if np.isnan(grids[m][i, j]) else float(grids[m][i, j]))
                                for j, c in enumerate(CLASSES)}
                            for i, (a, _l) in enumerate(ARMS)} for m, _t in MODELS},
           "base_f1": {m: {c: float(D[(m, "ch0", DROPIN, c)]["base_f1"])
                           for c in CLASSES if (m, "ch0", DROPIN, c) in D} for m, _t in MODELS},
           "n_ref": {m: {c: int(D[(m, "ch0", DROPIN, c)]["n_ref"])
                         for c in CLASSES if (m, "ch0", DROPIN, c) in D} for m, _t in MODELS}}
    (HERE / "C11_values.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/'C11_channel_ablation.pdf'}")
    for m, _t in MODELS:
        g = grids[m]
        i_int = CLASSES.index("internal")
        fr = [g[k, i_int] for k in range(3)]
        print(f"  {m}: internal delta for frame channels {['%+.3f' % x for x in fr]}, "
              f"ch3(start) on uORF {g[3, CLASSES.index('uORF')]:+.3f}")


if __name__ == "__main__":
    main()
