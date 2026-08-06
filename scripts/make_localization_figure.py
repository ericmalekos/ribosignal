#!/usr/bin/env python3
"""Localization eval figure: does the model concentrate predicted P-site signal on translated ORFs?

Three panels from one run's localization outputs (default: within-tissue orf_v2_attn fold-0):
  A  Predicted frame0 (periodicity) distribution for annotated-translated vs non-canonical-
     translated vs untranslated candidate AUG ORFs, with the 1/3 no-translation null. The
     separation IS the localization signal; non-canonical positives sit between canonical CDS
     and the untranslated sea.
  B  Discrimination AUROC on the NON-CANONICAL stratum (where length is uninformative): the
     ORF-length floor vs the model's length-controlled predicted-frame0 vs the observed-profile
     ceiling, baseline vs orf_v2_attn. Shows the model beats length and the ORF track + attention
     adds to it.
  C  Fidelity: median predicted vs observed frame0 per ORF type -- the model reproduces the 3-nt
     periodicity it was never explicitly trained on (loss is multinomial NLL on raw counts).

Panels A/C read <run>/localization_orfs.tsv; panel B reads localization_metrics.json for the two
compared runs. cas12a (matplotlib) on CPU. Writes <out> (default figures/localization.png) + a
sibling FIGURE_DATA_INPUTS note is maintained separately.

Usage: make_localization_figure.py [--run <dir>] [--baseline <dir>] [--out <png>]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
NONCANON = {"uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "novel", "internal"}
TYPE_ORDER = ["annotated", "uORF", "novel", "dORF", "Overlap_uORF", "Overlap_dORF", "internal"]


def load_orfs(run):
    """Panel A/C inputs, matching the collapsed metric set: candidates that are INTERNAL in-frame
    alt starts (cds_rel == cds_inframe: an AUG re-initiating the same protein truncated inside the
    CDS) are dropped from both positives and negatives. dORFs (in-frame and off-frame), out-of-frame
    overlapping ORFs, upstream extensions, uORFs and lncRNA candidates all stay. (The untranslated-
    decomposition figure keeps the internal alt starts to show why they are collapsed.)"""
    pos_an, pos_nc, neg = [], [], []
    by_type = {}
    with (run / "localization_orfs.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            f0 = r["pred_frame0"]
            of0 = r["obs_frame0"]
            f0 = float(f0) if f0 != "nan" else np.nan
            of0 = float(of0) if of0 != "nan" else np.nan
            is_pos = r["is_translated"] == "1"
            # drop internal in-frame alt starts (except annotated positives) -- the collapsed class
            if r.get("cds_rel") == "cds_inframe" and not (is_pos and r["orf_type"] == "annotated"):
                continue
            if is_pos:
                (pos_nc if r["orf_type"] in NONCANON else pos_an).append(f0)
                by_type.setdefault(r["orf_type"], []).append((f0, of0))
            else:
                neg.append(f0)
    return pos_an, pos_nc, neg, by_type


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(NEW / "results" / "improve" / "orf_v2_attn_f0"))
    ap.add_argument("--baseline", default=str(NEW / "results" / "improve" / "baseline_f0"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    run = Path(args.run)
    base = Path(args.baseline)
    pos_an, pos_nc, neg, by_type = load_orfs(run)
    m_run = json.loads((run / "localization_metrics.json").read_text())
    m_base = json.loads((base / "localization_metrics.json").read_text())

    fig = plt.figure(figsize=(14, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.15], wspace=0.32)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    # Panel A: predicted frame0 distributions
    groups = [("untranslated\ncandidates", neg, "#9aa4ad"),
              ("non-canonical\ntranslated", pos_nc, "#d1782f"),
              ("annotated CDS\ntranslated", pos_an, "#2b6cb0")]
    parts = axA.violinplot([np.array([v for v in g if v == v]) for _, g, _ in groups],
                           showmedians=True, showextrema=False, widths=0.85)
    for pc_, (_, _, col) in zip(parts["bodies"], groups, strict=False):
        pc_.set_facecolor(col)
        pc_.set_alpha(0.75)
    parts["cmedians"].set_color("#222")
    axA.axhline(1 / 3, ls="--", lw=1, color="#444")
    axA.text(0.55, 1 / 3 + 0.015, "1/3 (no translation)", fontsize=8, color="#444")
    axA.set_xticks([1, 2, 3])
    axA.set_xticklabels([g[0] for g in groups], fontsize=8.5)
    axA.set_ylabel("predicted in-ORF frame-0 fraction")
    axA.set_ylim(0, 1)
    for i, (_, g, _) in enumerate(groups):
        gg = [v for v in g if v == v]
        axA.text(i + 1, 0.03, f"n={len(gg):,}\nmed {np.median(gg):.2f}",
                 ha="center", fontsize=7.5, color="#333")
    axA.set_title("A  Periodicity localizes to translated ORFs", fontsize=9.5, loc="left")

    # Panel B: non-canonical discrimination AUROC (length-controlled)
    def nc(m, key):
        d = m["discrimination"]["noncanonical"]
        return d.get(key, float("nan"))
    bars = [
        ("ORF length\n(floor)", nc(m_run, "auroc_orf_length_floor"),
         nc(m_run, "auroc_orf_length_floor"), "#9aa4ad"),
        ("pred frame0\n(len-ctrl)", nc(m_base, "auroc_pred_frame0_lengthctrl"),
         nc(m_run, "auroc_pred_frame0_lengthctrl"), "#d1782f"),
        ("observed\nceiling", nc(m_base, "auroc_obs_frame0_ceiling_lengthctrl"),
         nc(m_run, "auroc_obs_frame0_ceiling_lengthctrl"), "#2b6cb0"),
    ]
    x = np.arange(len(bars))
    w = 0.38
    b0 = axB.bar(x - w / 2, [b[1] for b in bars], w, color=[b[3] for b in bars],
                 alpha=0.5, label="baseline", edgecolor="#555", linewidth=0.6)
    b1 = axB.bar(x + w / 2, [b[2] for b in bars], w, color=[b[3] for b in bars],
                 alpha=1.0, label="orf_v2_attn", edgecolor="#222", linewidth=0.6)
    for bars_ in (b0, b1):
        for rect in bars_:
            h = rect.get_height()
            if h == h:
                axB.text(rect.get_x() + rect.get_width() / 2, h + 0.008, f"{h:.2f}",
                         ha="center", fontsize=7)
    axB.axhline(0.5, ls=":", lw=1, color="#777")
    axB.set_xticks(x)
    axB.set_xticklabels([b[0] for b in bars], fontsize=8)
    axB.set_ylabel("AUROC (translated vs untranslated)")
    axB.set_ylim(0.4, max(0.8, np.nanmax([b[2] for b in bars]) + 0.06))
    axB.legend(fontsize=8, loc="upper left", framealpha=0.9)
    axB.set_title("B  Non-canonical: beats length floor", fontsize=9.5, loc="left")

    # Panel C: predicted vs observed frame0 per ORF type (fidelity)
    types = [t for t in TYPE_ORDER if t in by_type and len(by_type[t]) >= 20]
    pred_med, obs_med, ns = [], [], []
    for t in types:
        pv = [p for p, o in by_type[t] if p == p]
        ov = [o for p, o in by_type[t] if o == o]
        pred_med.append(np.median(pv) if pv else np.nan)
        obs_med.append(np.median(ov) if ov else np.nan)
        ns.append(len(by_type[t]))
    y = np.arange(len(types))
    axC.plot(obs_med, y, "o", color="#2b6cb0", ms=7, label="observed")
    axC.plot(pred_med, y, "D", color="#d1782f", ms=6, label="predicted")
    for i in range(len(types)):
        axC.plot([obs_med[i], pred_med[i]], [i, i], "-", color="#bbb", lw=1, zorder=0)
    axC.axvline(1 / 3, ls="--", lw=1, color="#444")
    axC.set_yticks(y)
    axC.set_yticklabels([f"{t}\n(n={n})" for t, n in zip(types, ns, strict=False)], fontsize=8)
    axC.invert_yaxis()
    axC.set_xlabel("median in-ORF frame-0 fraction")
    axC.set_xlim(0, 1)
    axC.legend(fontsize=8, loc="lower right", framealpha=0.9)
    axC.set_title("C  Reproduces observed periodicity per type", fontsize=9.5, loc="left")

    fig.suptitle(f"Translated-ORF localization on held-out {m_run['tissue']} "
                 f"({m_run['run']} vs {m_base['run']})", fontsize=11, y=1.02)
    out = Path(args.out) if args.out else (run / "figures" / "localization.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
