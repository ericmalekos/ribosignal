#!/usr/bin/env python3
"""GTF2I: the corrupted input, the corrected input, and what each model does with them.

The whole argument in one panel. Rows top to bottom:

  1-2  RNA coverage as the model actually receives it (the pack's own coverage.npy),
       mm1 vs mm10. This is a MODEL INPUT, and mm1 leaves a 1,617 nt hole in the CDS.
  3-6  Each model's predicted profile given each coverage. IDENTICAL WEIGHTS within a
       model -- nothing is retrained, only which pack supplies the coverage channel.
  7    The observed Ribo-seq target, which is mm1 and therefore holed in the same window.

Row 7 is why this is not simply "the model got better": the measurement itself is blind
here (1 P-site of 13,258 in the gap), so there is no unholed ground truth for this window.
What the panel shows is that the models READ their coverage input rather than reciting a
memorised output -- the hole in rows 3 and 5 is inherited from row 1, not from training.

  make_gtf2i_inputswap.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NEW = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
G = NEW / "data/mm25_diagnostic/gtf2i_inference"
TX, CS, CE = "ENST00000901263.1", 451, 3324
GAP = (1682, 3299)


def pack_arrays(pack, name):
    order = (NEW / "data" / pack / "tx_order.txt").read_text().split()
    i = order.index(TX)
    off = np.load(NEW / "data" / pack / "offsets.npy")
    a, b = int(off[i]), int(off[i + 1])
    return np.load(NEW / "data" / pack / name, mmap_mode="r")[a:b].astype(float)


def pred(dump):
    z = np.load(G / dump / "pred_profiles.npz")
    ids = [str(x) for x in z["tx_ids"]]
    i = ids.index(TX)
    off = np.concatenate([[0], np.cumsum(z["lengths"])])
    return z["pred_flat"][int(off[i]):int(off[i + 1])].astype(float)


def main():
    g0, g1 = GAP
    cov1 = pack_arrays("packed_union_Hepatocytes", "coverage.npy")
    cov10 = pack_arrays("packed_mm10cov_Hepatocytes", "coverage.npy")
    obs = pack_arrays("packed_union_Hepatocytes", "target_counts.npy")
    rows = [
        ("RNA coverage input, mm1  (multimappers DISCARDED)", cov1, "#B4574A", "input"),
        ("RNA coverage input, mm10 (multimappers KEPT)", cov10, "#C8A45B", "input"),
        ("attn prediction  GIVEN mm1 coverage", pred("attn_dump_union"), "#B4574A", "pred"),
        ("attn prediction  GIVEN mm10 coverage", pred("attn_dump_mm10cov"), "#2C6FBB", "pred"),
        ("mamba4 prediction  GIVEN mm1 coverage", pred("mamba4_dump_union"), "#B4574A", "pred"),
        ("mamba4 prediction  GIVEN mm10 coverage", pred("mamba4_dump_mm10cov"), "#2C6FBB", "pred"),
        ("observed Ribo-seq (the TARGET, also mm1 -- blind here)", obs, "#555555", "obs"),
    ]
    fig, axes = plt.subplots(len(rows), 1, figsize=(13.5, 12), sharex=True)
    vals = {}
    for ax, (lab, y, c, kind) in zip(axes, rows, strict=False):
        x = np.arange(len(y))
        ax.fill_between(x, 0, y, color=c, lw=0, zorder=3)
        ax.axvspan(g0, g1, color="#D94A3D", alpha=0.10, zorder=1)
        for v in (CS, CE):
            ax.axvline(v, color="#333", lw=0.8, ls=":", zorder=4)
        gm = float(y[g0:g1].mean())
        cm = float(y[CS:CE].mean())
        pctv = 100 * gm / cm if cm else float("nan")
        # predictions live at ~1e-3, so a fixed 1-decimal mean prints an uninformative "0.0/nt"
        gs = f"{gm:,.1f}" if gm >= 1 else f"{gm:.2e}"
        ax.set_title(f"{lab}        gap mean {gs}/nt   =  {pctv:.1f}% of CDS mean",
                     fontsize=9, loc="left")
        ax.set_ylabel({"input": "depth", "pred": "pred", "obs": "P-sites"}[kind], fontsize=8)
        ax.set_xlim(0, len(y))
        ax.tick_params(labelsize=8)
        vals[lab] = {"gap_mean_per_nt": round(gm, 4), "cds_mean_per_nt": round(cm, 4),
                     "gap_pct_of_cds": round(pctv, 1)}
    axes[-1].set_xlabel(f"{TX} transcript position (nt).  dotted = CDS {CS}-{CE};  "
                        f"red band = the {g1-g0:,} nt multimapped window", fontsize=9)
    fig.suptitle("GTF2I: the model reads its RNA input. Same weights, "
                 "only the coverage channel differs.", fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"X_gtf2i_inputswap.{ext}", dpi=200)
    (OUT / "X_gtf2i_inputswap_values.json").write_text(json.dumps(
        {"model": ["attn", "mamba4"],
         "source": "RNA + observed from data/packed_{union,mm10cov}_Hepatocytes; predictions "
                   "from data/mm25_diagnostic/gtf2i_inference/{attn,mamba4}_dump_"
                   "{union,mm10cov}/pred_profiles.npz (dump_pred_profiles.py, input-swap: "
                   "identical weights, only RIBO_PACK_SUFFIX differs -- NOT a retrain)",
         "transcript": TX, "cds": [CS, CE], "gap": list(GAP), "rows": vals}, indent=2) + "\n")
    print(f"  wrote {OUT/'X_gtf2i_inputswap.png'}")
    for k, v in vals.items():
        print(f"    {v['gap_pct_of_cds']:7.1f}% of CDS mean   {k}")


if __name__ == "__main__":
    main()
