#!/usr/bin/env python3
"""Example predictions, OBSERVED and PREDICTED drawn as separate frame-coloured bar panels.

Companion to plot_prediction_examples.py, which overlays the two (predicted as bars, observed as a
black step line). Overlaying makes the comparison compact but asymmetric: the observed trace is a
single colour, so the reader can only judge the PREDICTED profile's reading frame, and has to take
on trust that the observed signal is periodic in the same frame.

Here each series gets its own panel, coloured by codon frame on the same rule, stacked so the
observed panel sits directly above its prediction and the two share an x and y axis. Frame-0
dominance is then legible independently in each, and the comparison is between two like objects.

Layout: 3 columns (uORF / CDS / lncRNA), 2 rows (observed on top, predicted beneath). Selection,
zoom-window choice and all quality gates are imported from plot_prediction_examples, so the two
figures always showcase the same transcripts for a given dump.

  python -m scripts.plot_prediction_examples_split --npz <pred_profiles.npz> --auto --out <png>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_prediction_examples import (  # noqa: E402
    CALLS, FRAME_C, NEW, NPZ, ZOOM_W, autoselect, best_zoom, load_tsv, orf_windows, pearson,
    region_quality,
)

OBS_EDGE = "#111111"


def frame_bars(ax, lo, hi, y, ref_start, title, ylabel, edge=None):
    """Per-nucleotide bars over [lo, hi), coloured by codon frame relative to `ref_start`."""
    x = np.arange(lo, hi)
    frames = (x - ref_start) % 3
    ax.bar(x, y[lo:hi], width=0.9, linewidth=(0.4 if edge else 0),
           edgecolor=(edge or "none"), color=[FRAME_C[f] for f in frames])
    ax.margins(x=0)
    ax.set_title(title, fontsize=8.5)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", default=str(NPZ))
    ap.add_argument("--auto", action="store_true", help="re-select exemplars for this dump")
    ap.add_argument("--uorf_tx", default="ENST00000399697.7")
    ap.add_argument("--cds_tx", default="ENST00000429711.7")
    ap.add_argument("--lnc_tx", default="ENST00000647872.1")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default=str(NEW / "figures/prediction_examples/"
                                             "prediction_examples_split.png"))
    a = ap.parse_args()

    z = np.load(a.npz, allow_pickle=False)
    tx_ids = [str(t) for t in z["tx_ids"]]
    idx = {t: i for i, t in enumerate(tx_ids)}
    lengths = z["lengths"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    pred, obs = z["pred_flat"], z["obs_flat"]
    bio = load_tsv(NEW / "data/tx2biotype.tsv", "tx_id", ("gene_name", "transcript_type"))
    cds = load_tsv(NEW / "data/tx2cds.tsv", "tx_id", ("utr5_len", "cds_len", "utr3_len"))
    calls = orf_windows(CALLS)

    if a.auto:
        sel = autoselect(tx_ids, idx, off, pred, obs, bio, cds, calls)
        missing = [k for k, (t, _r) in sel.items() if t is None]
        if missing:
            raise SystemExit(f"auto-selection found no candidate for: {missing}")
        a.uorf_tx, a.cds_tx, a.lnc_tx = sel["uorf"][0], sel["cds"][0], sel["lnc"][0]
        for k in ("uorf", "cds", "lnc"):
            t, r = sel[k]
            print(f"auto-selected {k:5s}: {t} {bio[t][0]:<14s} r = {r:.3f}")

    def get(tx):
        j = idx[tx]
        p = pred[off[j]:off[j + 1]].astype(np.float64)
        c = obs[off[j]:off[j + 1]].astype(np.float64)
        return p, c, c / c.sum()

    # (label, tx, ORF window, reference start for the frame colouring)
    cols = []
    tx = a.uorf_tx
    p, c, q = get(tx)
    uo = next(w for w in calls[tx] if w[0] in ("uORF", "Overlap_uORF"))
    ua, ue = uo[1] - 1, uo[2]
    cols.append(("uORF", tx, p, c, q, ua, ue, pearson(p[:int(cds[tx][0])], q[:int(cds[tx][0])])))

    tx = a.cds_tx
    p, c, q = get(tx)
    u5, cl, _ = (int(v) for v in cds[tx])
    cols.append(("CDS", tx, p, c, q, u5, u5 + cl, pearson(p, q)))

    tx = a.lnc_tx
    p, c, q = get(tx)
    lo_ = next(w for w in calls[tx] if w[0] in ("novel", "novel_NonPCGs"))
    cols.append(("lncRNA ORF", tx, p, c, q, lo_[1] - 1, lo_[2], pearson(p, q)))

    fig, axes = plt.subplots(2, 3, figsize=(15, 6.2))
    for k, (lab, tx, p, c, q, oa, oe, r) in enumerate(cols):
        w = oe - oa
        if w <= ZOOM_W:
            pad = max(6, (ZOOM_W - w) // 2)
            zlo, zhi = max(0, oa - pad), min(len(c), oe + pad)
        else:
            zlo, zhi = best_zoom(c, oa, oe)
        gn = bio[tx][0]
        f0_obs = region_quality(c, zlo, zhi, oa)[2]
        f0_pred = region_quality(p, zlo, zhi, oa)[2]

        frame_bars(axes[0, k], zlo, zhi, q, oa,
                   f"{lab}  |  {gn}   OBSERVED   frame-0 = {f0_obs:.0%}",
                   "P-site fraction", edge=OBS_EDGE)
        frame_bars(axes[1, k], zlo, zhi, p, oa,
                   f"PREDICTED   frame-0 = {f0_pred:.0%}   (profile r = {r:.3f})",
                   "P-site fraction")
        axes[1, k].set_xlabel(f"transcript position (nt)   {tx}", fontsize=8)

        # same y on both panels so the two series are compared on one scale, not auto-rescaled
        top = max(q[zlo:zhi].max(), p[zlo:zhi].max()) * 1.12
        for row in (0, 1):
            axes[row, k].set_ylim(0, top)
            axes[row, k].set_xlim(zlo, zhi)
            for b, style in ((oa, "--"), (oe, ":")):
                if zlo <= b <= zhi:
                    axes[row, k].axvline(b, color="#555", ls=style, lw=0.9)

    handles = [Patch(color=FRAME_C[f], label=f"frame {f}" + (" (in-frame)" if f == 0 else ""))
               for f in range(3)]
    handles.append(Line2D([0], [0], color=OBS_EDGE, lw=0.8, label="observed bars are outlined"))
    fig.legend(handles=handles, fontsize=8, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(a.title or ("Observed vs predicted Ribo-seq P-site profile, coloured by codon "
                             "frame -- held-out Hepatocytes"), fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0.03, 1, 0.98))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    pdf = out.with_suffix(".pdf")
    if pdf != out:
        fig.savefig(pdf, bbox_inches="tight")
        print(f"wrote {pdf}")


if __name__ == "__main__":
    main()
