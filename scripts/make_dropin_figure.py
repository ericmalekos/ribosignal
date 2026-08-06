#!/usr/bin/env python3
"""Task 14 figure: do RiboCode's ORF calls on the model's PREDICTED density track its calls on the
real Ribo-seq? Two panels from dropin_metrics.json (compare_dropin_calls.py output):
  A  precision / recall / F1 for the harness check (real vs official) and the two predicted-density
     drop-ins (predicted shape at real depth; fully standalone at predicted depth), all vs the real
     RiboCode calls. Recall = fraction of real calls recovered by the predicted density.
  B  recall by ORF type (predicted shape vs real; standalone overlaid) -- canonical CDS and uORFs
     transfer strongly, dORF and internal are the blind spots (matching the Task 13 eval).
cas12a matplotlib. Usage: make_dropin_figure.py [--metrics <json>] [--out <png>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
DEFAULT = NEW / "results/loto/orf_v2_attn_holdout_Hepatocytes/dropin/dropin_metrics.json"
TYPE_ORDER = ["annotated", "uORF", "novel", "Overlap_uORF", "Overlap_dORF", "dORF", "internal"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default=str(DEFAULT))
    ap.add_argument("--out", default=str(NEW / "results/figures/dropin_ribocode.png"))
    args = ap.parse_args()
    m = json.loads(Path(args.metrics).read_text())

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 4.8),
                                   gridspec_kw={"width_ratios": [1.05, 1.25]})

    # Panel A: precision / recall / F1 for the three comparisons
    comps = [("real_vs_official", "real ->\nofficial\n(harness)", "#7f8c8d"),
             ("pred_obsdepth_vs_real", "pred shape\n-> real", "#d1782f"),
             ("pred_preddepth_vs_real", "pred standalone\n-> real", "#2b6cb0")]
    metrics = ["precision", "recall", "f1"]
    x = np.arange(len(comps))
    w = 0.26
    for j, met in enumerate(metrics):
        vals = [m[c[0]][met] for c in comps]
        bars = axA.bar(x + (j - 1) * w, vals, w, label=met,
                       color=["#c9ced2", "#e8a877", "#7ea6d4"][j] if met != "recall"
                       else [c[2] for c in comps],
                       edgecolor="#333", linewidth=0.6, alpha=0.95 if met == "recall" else 0.8)
        for rect, v in zip(bars, vals, strict=False):
            axA.text(rect.get_x() + rect.get_width() / 2, v + 0.012, f"{v:.2f}",
                     ha="center", fontsize=7.2)
    axA.axhline(0.5, ls=":", lw=1, color="#aaa")
    axA.set_xticks(x)
    axA.set_xticklabels([c[1] for c in comps], fontsize=8.5)
    axA.set_ylabel("agreement with real-density RiboCode calls")
    axA.set_ylim(0, 1.0)
    axA.legend(fontsize=8, loc="upper right", ncol=3, framealpha=0.9,
               columnspacing=0.9, handlelength=1.1)
    axA.set_title("A  Predicted density recovers most real ORF calls", fontsize=9.5, loc="left")

    # Panel B: recall by ORF type, pred shape vs standalone
    rt_obs = m["pred_obsdepth_vs_real"]["recall_per_type"]
    rt_pred = m["pred_preddepth_vs_real"]["recall_per_type"]
    types = [t for t in TYPE_ORDER if t in rt_obs]
    y = np.arange(len(types))
    h = 0.38
    r_obs = [rt_obs[t]["recall"] for t in types]
    r_pred = [rt_pred.get(t, {"recall": np.nan})["recall"] for t in types]
    ns = [rt_obs[t]["n_real"] for t in types]
    axB.barh(y + h / 2, r_obs, h, color="#d1782f", edgecolor="#333", linewidth=0.5,
             label="pred shape (real depth)")
    axB.barh(y - h / 2, r_pred, h, color="#2b6cb0", edgecolor="#333", linewidth=0.5,
             label="pred standalone")
    for yi, (ro, n) in enumerate(zip(r_obs, ns, strict=False)):
        axB.text(ro + 0.012, yi + h / 2, f"{ro:.2f}", va="center", fontsize=7)
        axB.text(1.0, yi, f"n={n:,}", va="center", ha="right", fontsize=6.8, color="#888")
    axB.axvline(0.5, ls=":", lw=1, color="#aaa")
    axB.set_yticks(y)
    axB.set_yticklabels([f"{t}\n({'canonical' if t == 'annotated' else 'non-canonical'})"
                         for t in types], fontsize=7.8)
    axB.invert_yaxis()
    axB.set_xlim(0, 1.05)
    axB.set_xlabel("recall (fraction of real ORF calls recovered)")
    axB.legend(fontsize=8, loc="lower right", framealpha=0.9)
    axB.set_title("B  Recall by ORF type (canonical + uORF transfer; dORF/internal blind spots)",
                  fontsize=9.0, loc="left")

    ml = m.get("min_len", 0)
    key = m.get("key", "genomic")
    me = m.get("min_enrichment", 0)
    keylab = "genomic-locus matched" if key == "genomic" else "transcript matched"
    flt = f", ORFs >= {ml} nt" if ml else ""
    enr = f", pred enrichment >= {me:g}" if me else ""
    fig.suptitle("RiboCode drop-in: ORF calls on the predicted density vs the real Ribo-seq "
                 f"(Hepatocytes, orf_v2_attn{flt}{enr}, {keylab})", fontsize=10.5, y=1.02)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
