#!/usr/bin/env python3
"""GTF2I RNA coverage under five STAR multimap postures, against the Ribo-seq the model must predict.

THE QUESTION. GTF2I ENST00000901263.1 has a 1,617 nt stretch of its CDS (1682-3299) where observed
Ribo-seq is ~0 and the model also predicts ~0. That reads as memorisation ONLY if the model had no
input telling it to. But the pack's RNA coverage -- a MODEL INPUT -- is 11x depleted in the same
window (12.6 vs 141.5 per nt over the CDS), and that channel is `--outFilterMultimapNmax 1` too.

If the hole fills in when multimappers are kept, the model was reading a corrupted feature
correctly and the fix is the coverage channel, not the training set.

  make_gtf2i_rna_postures.py
"""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402
import numpy as np                      # noqa: E402

NEW = Path(__file__).resolve().parents[2]
D = NEW / "data/mm25_diagnostic/hep_rna_mm/gtf2i"
OUT = Path(__file__).resolve().parent
TX, L, CS, CE = "ENST00000901263.1", 5178, 451, 3324
GAP = (1682, 3299)
POST = [("mm1", "mm1  (production: multimappers DISCARDED)", "#B4574A"),
        ("mm10", "mm10 (all alignments kept)", "#C8A45B"),
        ("mm25", "mm25 (all alignments kept)", "#7BA05B"),
        ("best_Old_2_4", "best-only, Old_2.4 tie order", "#2C6FBB"),
        ("best_Random", "best-only, Random tie order", "#6B4C9A")]


def load(tag):
    tot = np.zeros(L)
    n = 0
    for srr in ("SRR15513268", "SRR15513269"):
        f = D / f"{srr}.{tag}.cov.tsv"
        if not f.exists():
            continue
        v = np.loadtxt(f, dtype=float, usecols=1)
        tot[:len(v)] += v[:L]; n += 1
    return tot / max(n, 1)


def main():
    cov = {t: load(t) for t, _, _ in POST}
    # the Ribo-seq the model is asked to predict, same transcript, same tissue
    P = NEW / "data/packed_union_Hepatocytes"
    order = (P / "tx_order.txt").read_text().split(); i = order.index(TX)
    off = np.load(P / "offsets.npy")
    obs = np.load(P / "target_counts.npy")[off[i]:off[i+1]].astype(float)

    fig, axes = plt.subplots(6, 1, figsize=(13.5, 10), sharex=True)
    g0, g1 = GAP
    for ax, (tag, lab, c) in zip(axes, POST):
        y = cov[tag]
        ax.axvspan(CS-1, CE-1, color="#000000", alpha=0.045, lw=0)
        ax.axvspan(g0, g1, color="#B4574A", alpha=0.10, lw=0)
        ax.fill_between(np.arange(L), 0, y, color=c, lw=0)
        ing = y[g0:g1].mean(); incds = y[CS-1:CE-1].mean()
        ax.set_ylabel("RNA cov", fontsize=8)
        ax.set_title(f"{lab}     gap mean {ing:,.0f}/nt   vs CDS mean {incds:,.0f}/nt"
                     f"   ({100*ing/max(incds,1e-9):.0f}% of CDS)", fontsize=9, loc="left", pad=3)
        ax.tick_params(labelsize=7)
    ax = axes[-1]
    ax.axvspan(CS-1, CE-1, color="#000000", alpha=0.045, lw=0)
    ax.axvspan(g0, g1, color="#B4574A", alpha=0.10, lw=0)
    ax.vlines(np.arange(L), 0, obs, color="#333333", lw=0.5)
    ax.set_ylabel("Ribo P-sites", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title(f"observed Ribo-seq (the target)     gap total {int(obs[g0:g1].sum())} of "
                 f"{int(obs.sum()):,}", fontsize=9, loc="left", pad=3)
    ax.set_xlabel("position along transcript (nt)   |   grey = CDS 451-3324   |   "
                  "red = the 1,617 nt gap", fontsize=9)
    ax.set_xlim(0, L)
    fig.suptitle(f"GTF2I {TX}, held-out hepatocytes: does the RNA-input hole survive keeping "
                 f"multimappers?", fontsize=11.5, y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"X_gtf2i_rna_postures.{e}", dpi=175, bbox_inches="tight")

    vals = {"model": None,
            "source": "data/mm25_diagnostic/hep_rna_mm/gtf2i/*.cov.tsv (STAR postures on the 2 "
                      "Hepatocytes RNA libs); observed from data/packed_union_Hepatocytes",
            "provenance_note": "MODEL-FREE. RNA coverage is a model INPUT; mm10/mm25 keep all "
                               "alignments and are count-inflating, best-only commits each read to "
                               "one alignment.",
            "transcript": TX, "cds": [CS, CE], "gap": list(GAP),
            "postures": {tag: {"gap_mean_per_nt": round(float(cov[tag][g0:g1].mean()), 2),
                               "cds_mean_per_nt": round(float(cov[tag][CS-1:CE-1].mean()), 2),
                               "gap_pct_of_cds": round(100*float(cov[tag][g0:g1].mean())
                                                       / max(float(cov[tag][CS-1:CE-1].mean()), 1e-9), 1),
                               "total": int(cov[tag].sum())} for tag, _, _ in POST},
            "observed_ribo_in_gap": int(obs[g0:g1].sum()), "observed_ribo_total": int(obs.sum())}
    (OUT / "X_gtf2i_rna_postures_values.json").write_text(json.dumps(vals, indent=2))
    print(f"  {'posture':<16} {'gap/nt':>9} {'CDS/nt':>9} {'gap as % of CDS':>17}")
    for tag, _, _ in POST:
        v = vals["postures"][tag]
        print(f"  {tag:<16} {v['gap_mean_per_nt']:>9,.0f} {v['cds_mean_per_nt']:>9,.0f} "
              f"{v['gap_pct_of_cds']:>16.0f}%")
    print(f"\n  observed Ribo-seq in the gap: {vals['observed_ribo_in_gap']} of "
          f"{vals['observed_ribo_total']:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
