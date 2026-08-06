#!/usr/bin/env python3
"""Supplement to the localization figure: WHY is the untranslated-candidate frame0 distribution
(main-figure panel A) bimodal, with a large cluster above 0.7?

Decomposes the untranslated candidate AUG ORFs by their relationship to the annotated CDS (the
`cds_rel` column that eval_localization.py now writes into localization_orfs.tsv): in-frame CDS
overlap / off-frame CDS overlap / 5'UTR / 3'UTR / no-CDS(lncRNA). The high-periodicity cluster is
almost entirely `cds_inframe` -- candidate ORFs that share the reading frame of a genuinely
translated CDS (RiboCode collapses each locus to one call, so these in-frame fragments are labeled
'untranslated' but sit inside real periodic translation). The model predicts high frame0 for them
correctly; observed frame0 (overlaid) confirms they really are periodic.

Panels: (left) predicted in-ORF frame0 violins per cds_rel with observed-frame0 medians overlaid;
(right) cds_rel composition of the high (pred frame0 > 0.7) vs low (< 0.45) clusters.

Reads <run>/localization_orfs.tsv (default within-tissue orf_v2_attn fold-0). cas12a matplotlib.
Usage: make_untranslated_decomp_figure.py [--run <dir>] [--out <png>]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ORDER = ["cds_inframe", "dorf_inframe", "cds_offframe", "ext_inframe", "utr5", "utr3", "no_cds"]
LABELS = {"cds_inframe": "internal\nin-frame\nalt start",
          "dorf_inframe": "in-frame\ndORF\n(past stop)",
          "cds_offframe": "off-frame\nCDS overlap", "ext_inframe": "in-frame\nN-term\nextension",
          "utr5": "5'UTR", "utr3": "3'UTR", "no_cds": "no CDS\n(lncRNA)"}
# red = the ONLY collapsed class (internal in-frame alt start); all others kept as candidates
COLORS = {"cds_inframe": "#c0392b", "dorf_inframe": "#7f8c8d", "cds_offframe": "#95a5a6",
          "ext_inframe": "#c39bd3", "utr5": "#e08e0b", "utr3": "#aab7b8", "no_cds": "#8e6fb0"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(NEW / "results" / "improve" / "orf_v2_attn_f0"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    run = Path(args.run)
    rows = [r for r in csv.DictReader((run / "localization_orfs.tsv").open(), delimiter="\t")
            if r["is_translated"] == "0"]
    if not rows or "cds_rel" not in rows[0]:
        raise SystemExit("localization_orfs.tsv lacks cds_rel; re-run eval_localization.py")

    def col(rs, k):
        return np.array([float(r[k]) if r[k] not in ("nan", "") else np.nan for r in rs])

    rel = np.array([r["cds_rel"] for r in rows], dtype=object)
    pf0 = col(rows, "pred_frame0")
    of0 = col(rows, "obs_frame0")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})

    # left: predicted frame0 violins per cds_rel, observed median overlaid
    present = [c for c in ORDER if (rel == c).sum() >= 20]
    data = [pf0[(rel == c)][~np.isnan(pf0[(rel == c)])] for c in present]
    parts = axL.violinplot(data, showmedians=True, showextrema=False, widths=0.85)
    for body, c in zip(parts["bodies"], present, strict=False):
        body.set_facecolor(COLORS[c])
        body.set_alpha(0.75)
    parts["cmedians"].set_color("#111")
    for i, c in enumerate(present):
        m = rel == c
        obs = np.nanmedian(of0[m])
        axL.plot(i + 1, obs, "D", color="#1b4f72", ms=7, zorder=5)
        axL.text(i + 1, 0.03, f"n={int(m.sum()):,}\n>{0.7:.0%}: {np.mean(pf0[m] > 0.7):.0%}",
                 ha="center", fontsize=7.5, color="#333")
    axL.axhline(1 / 3, ls="--", lw=1, color="#444")
    axL.text(len(present) - 0.4, 1 / 3 + 0.015, "1/3 null", fontsize=8, color="#444")
    axL.plot([], [], "D", color="#1b4f72", ms=7, label="observed frame0 (median)")
    axL.set_xticks(range(1, len(present) + 1))
    axL.set_xticklabels([LABELS[c] for c in present], fontsize=8.5)
    axL.set_ylabel("predicted in-ORF frame-0 fraction")
    axL.set_ylim(0, 1)
    axL.legend(fontsize=8, loc="upper right", framealpha=0.9)
    axL.set_title("Untranslated candidates by CDS relationship\n(the high cluster IS in-frame CDS "
                  "overlap)", fontsize=10, loc="left")

    # right: cds_rel composition of high vs low predicted-frame0 clusters
    hi, lo = pf0 > 0.7, pf0 < 0.45
    comp = {}
    for c in ORDER:
        comp[c] = ((rel == c) & hi).sum(), ((rel == c) & lo).sum()
    xs = ["pred frame0\n> 0.7\n(high cluster)", "pred frame0\n< 0.45\n(low cluster)"]
    bottoms = np.zeros(2)
    for c in ORDER:
        vals = np.array([comp[c][0] / max(hi.sum(), 1), comp[c][1] / max(lo.sum(), 1)])
        axR.bar(xs, vals, bottom=bottoms, color=COLORS[c], label=LABELS[c].replace("\n", " "),
                edgecolor="white", linewidth=0.5)
        for j, v in enumerate(vals):
            if v > 0.06:
                axR.text(j, bottoms[j] + v / 2, f"{v:.0%}", ha="center", va="center",
                         fontsize=8, color="white", fontweight="bold")
        bottoms += vals
    axR.set_ylabel("fraction of cluster")
    axR.set_ylim(0, 1)
    axR.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(1.0, 0.5), framealpha=0.9)
    axR.set_title("Cluster composition", fontsize=10, loc="left")

    fig.suptitle("Why the untranslated-candidate periodicity is bimodal "
                 "(within-tissue Fibroblast, orf_v2_attn)", fontsize=11, y=1.02)
    fig.tight_layout()
    out = Path(args.out) if args.out else (NEW / "results" / "figures" / "untranslated_decomp.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
