#!/usr/bin/env python3
"""Architecture / training-pipeline figure for the one-hot attention model (orf_v2_attn, union universe),
in THREE alternative presentation styles so one can be picked and iterated on:

  style1_panel : BPNet-style multi-panel journal figure (a-d). Dense, compact, muted. Input tracks ->
                 dilated-residual body + dual heads -> real observed-vs-predicted profile -> eval bars.
  style2_flow  : RiNALMo-style wide left-to-right pipeline. Training data -> central model tower
                 (exploded block detail) -> dual heads -> downstream evaluations fanning out right.
  style3_tiers : Orthrus-style A/B/C tiers. (A) data construction + LOTO split, (B) training pipeline
                 (10-track encoding -> body -> heads -> loss), (C) evaluation breadth matrix.

Panel (c)/profile traces use REAL data: predicted + observed per-nt P-site profiles are read straight from
the union model's held-out Hepatocytes dump (pred_profiles.npz carries both pred_flat and obs_flat), so the
periodicity shown is measured, not drawn. Every number printed on the figures is sourced in
FIGURE_DATA_INPUTS.md. cas12a env (matplotlib + numpy). Run: python make_arch_attn.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/arch_attn"
RUN = NEW / "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
NPZ = RUN / "dropin/pred_profiles.npz"
TX2CDS = NEW / "data/tx2cds.tsv"

# ---- model facts (args.json / history.json / test_metrics.json; see FIGURE_DATA_INPUTS.md) ----
M = dict(params="5,071,106", channels=256, n_blocks=10, n_attn=2, n_heads=8,
         epochs=28, hours=14.2, gpu="1x A5500", best_epoch=21, val_pearson=0.610,
         test_pearson=0.659, budget=16000, lr="3e-4", count_w=0.1)
TRAIN_TISSUES = ["Fibroblast", "VSMC", "ES", "Fat", "HA_EC", "HCAEC", "HUVEC"]
DILATIONS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
EVALS = [
    ("In-distribution", "held-out tissue\n(Hepatocytes)", "profile r", 0.659),
    ("ORF calling", "RiboCode drop-in\nCDS F1", "F1", 0.953),
    ("Cross-study", "Ruiz-Orera\n(human)", "CDS F1", 0.931),
    ("Cross-species", "5 mouse sets\n(median)", "CDS F1", 0.971),
    ("Proteogenomics", "12 macrophage\npopulations", "MS validated", None),
]
XSPECIES = [("Janich liver", 0.995), ("T-cell", 0.980), ("GSE120762 LPS", 0.971),
            ("Wang liver", 0.962), ("GSE120762 NT", 0.960)]


# ------------------------------------------------------------------ real data
def load_example(min_len=900, max_len=2400, min_counts=400):
    """Pick a deterministic well-covered held-out transcript; return (tx, obs, pred, cds, r)."""
    d = np.load(NPZ, allow_pickle=False)
    tx_ids = [str(t) for t in d["tx_ids"]]
    lengths = d["lengths"].astype(np.int64)
    off = np.concatenate([[0], np.cumsum(lengths)])
    obs_flat, pred_flat = d["obs_flat"], d["pred_flat"]
    totals = np.add.reduceat(obs_flat, off[:-1])                     # cheap per-tx sums
    cds = {}
    with open(TX2CDS) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len")}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            cl = int(f[ix["cds_len"]])
            if cl > 0:
                u5 = int(f[ix["utr5_len"]])
                cds[f[ix["tx_id"]]] = (u5, u5 + cl)
    ok = np.where((lengths >= min_len) & (lengths <= max_len) & (totals >= min_counts))[0]
    best = None
    for i in sorted(ok, key=lambda j: -totals[j])[:400]:             # deterministic scan
        tx = tx_ids[i]
        if tx not in cds:
            continue
        a, b = off[i], off[i + 1]
        o, p = obs_flat[a:b].astype(float), pred_flat[a:b].astype(float)
        if o.std() == 0 or p.std() == 0:
            continue
        r = float(np.corrcoef(o, p)[0, 1])
        if best is None or r > best[-1]:
            best = (tx, o, p, cds[tx], r)
        if best[-1] > 0.90:                                          # good enough, stop early
            break
    return best


# ------------------------------------------------------------------ drawing helpers
def box(ax, x, y, w, h, title, sub=None, fc="#ffffff", ec="#333333", tc="#111111",
        fs=8, subfs=6.5, lw=1.0, rad=0.012, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rad}",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ty = y + h / 2 + (0.012 if sub else 0)
    ax.text(x + w / 2, ty, title, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", zorder=3)
    if sub:
        ax.text(x + w / 2, y + h / 2 - 0.022, sub, ha="center", va="center",
                fontsize=subfs, color=tc, alpha=0.75, zorder=3)


def arrow(ax, x1, y1, x2, y2, c="#555555", lw=1.1, ms=7, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, color=c, lw=lw,
                                 mutation_scale=ms, shrinkA=0, shrinkB=0, zorder=4))


def dilated_motif(ax, x, y, w, h, rows=5, n=17, c="#2C6FBB", dot=2.2, lw=0.55):
    """Compact WaveNet/BPNet dilated-connectivity motif: each row doubles the stride."""
    xs = np.linspace(x, x + w, n)
    ys = np.linspace(y, y + h, rows)
    for r in range(rows):
        ax.plot(xs, [ys[r]] * n, "o", ms=dot, color=c, alpha=0.35 + 0.13 * r,
                mec="none", zorder=3)
    for r in range(rows - 1):
        d = 2 ** r
        for i in range(n):
            if i - d >= 0:
                ax.plot([xs[i], xs[i - d]], [ys[r], ys[r + 1]], "-", color=c,
                        lw=lw, alpha=0.5, zorder=2)
    return ys


def profile_axes(ax, obs, pred, cds, r, tx, c_obs="#9aa0a6", c_pred="#2C6FBB",
                 title=True, zoom=None):
    """Real observed vs predicted P-site profile; predicted rescaled to the observed total."""
    L = len(obs)
    pred_s = pred / pred.sum() * obs.sum() if pred.sum() > 0 else pred
    xs = np.arange(L)
    if zoom:
        s, e = zoom
        ax.fill_between(xs[s:e], obs[s:e], color=c_obs, lw=0, alpha=0.85, step="mid")
        ax.plot(xs[s:e], pred_s[s:e], color=c_pred, lw=0.9)
        ax.set_xlim(s, e)
    else:
        ax.fill_between(xs, obs, color=c_obs, lw=0, alpha=0.85, step="mid")
        ax.plot(xs, pred_s, color=c_pred, lw=0.7)
        ax.set_xlim(0, L)
        cs, ce = cds
        for v in (cs, ce):
            ax.axvline(v, color="#B5714F", lw=0.8, ls=(0, (3, 2)), zorder=5)
    ax.set_ylim(bottom=0)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=6, length=2, pad=1)
    ax.set_yticks([])
    if title:
        ax.set_title(f"{tx}   observed vs predicted   r = {r:.2f}", fontsize=7, pad=3)


def save(fig, name):
    HERE.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"{name}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


# ================================================================== STYLE 1
def style1(ex):
    """BPNet-style multi-panel journal figure: compact, muted, panel letters a-d."""
    INK, GREY, BLUE, ORANGE = "#1a1a1a", "#9aa0a6", "#2C6FBB", "#E08214"
    LIGHT, MID = "#eef2f6", "#dbe4ee"
    tx, obs, pred, cds, r = ex
    fig = plt.figure(figsize=(7.4, 7.0), facecolor="white")

    def lbl(x, y, s):
        fig.text(x, y, s, fontsize=10.5, fontweight="bold", color=INK, va="top")

    # ---- (a) inputs -------------------------------------------------------
    lbl(0.015, 0.985, "a")
    ax = fig.add_axes([0.07, 0.775, 0.88, 0.185]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0, 0.93, "Input: one transcript, per nucleotide (L up to 10,000 nt)",
            fontsize=8.5, fontweight="bold", color=INK)
    # transcript schematic
    ax.add_patch(Rectangle((0.06, 0.66), 0.16, 0.10, fc=MID, ec="#7f8c99", lw=0.7))
    ax.add_patch(Rectangle((0.22, 0.63), 0.44, 0.16, fc="#c9d8e8", ec="#5b7a99", lw=0.7))
    ax.add_patch(Rectangle((0.66, 0.66), 0.24, 0.10, fc=MID, ec="#7f8c99", lw=0.7))
    ax.text(0.14, 0.71, "5' UTR", ha="center", va="center", fontsize=6.5)
    ax.text(0.44, 0.71, "CDS", ha="center", va="center", fontsize=7, fontweight="bold")
    ax.text(0.78, 0.71, "3' UTR", ha="center", va="center", fontsize=6.5)
    rows = [("one-hot sequence", "A C G U", 4, BLUE),
            ("reading-frame track", "frame 0/1/2, start, stop", 5, ORANGE),
            ("RNA-seq coverage", "log1p depth", 1, "#6E8B74")]
    y0 = 0.50
    for i, (nm, sub, nch, c) in enumerate(rows):
        y = y0 - i * 0.175
        ax.add_patch(Rectangle((0.06, y), 0.84, 0.115, fc=c, ec="none", alpha=0.18))
        ax.text(0.075, y + 0.058, nm, fontsize=7.2, va="center", fontweight="bold", color=INK)
        ax.text(0.40, y + 0.058, sub, fontsize=6.4, va="center", color="#444444")
        ax.text(0.895, y + 0.058, f"{nch} ch", fontsize=6.6, va="center", ha="right", color=c,
                fontweight="bold")
    ax.annotate("", xy=(0.955, 0.16), xytext=(0.955, 0.615),
                arrowprops=dict(arrowstyle="-", color="#999999", lw=0.8))
    ax.text(0.965, 0.39, "10 channels", fontsize=6.6, rotation=90, va="center", color="#555555")

    # ---- (b) architecture -------------------------------------------------
    lbl(0.015, 0.755, "b")
    ax = fig.add_axes([0.07, 0.435, 0.88, 0.30]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    box(ax, 0.0, 0.40, 0.115, 0.20, "1x1 conv", "10 -> 256", fc=LIGHT, ec="#7f8c99", fs=7, subfs=6)
    arrow(ax, 0.118, 0.50, 0.155, 0.50)
    # dilated body
    ax.add_patch(FancyBboxPatch((0.16, 0.14), 0.36, 0.72,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                fc="#f7fafd", ec=BLUE, lw=1.0, zorder=1))
    ax.text(0.34, 0.805, "10 residual dilated conv blocks", ha="center", fontsize=7.6,
            fontweight="bold", color=INK)
    ax.text(0.34, 0.755, "Conv(k=3) - GroupNorm - GELU, residual", ha="center", fontsize=6.2,
            color="#555555")
    dilated_motif(ax, 0.19, 0.30, 0.30, 0.38, rows=5, n=17, c=BLUE)
    ax.text(0.34, 0.235, "dilation 1, 2, 4, ... 512   ->   receptive field ~4 kb",
            ha="center", fontsize=6.4, color="#444444", style="italic")
    ax.text(0.34, 0.175, "O(L) cost", ha="center", fontsize=6.2, color=BLUE, fontweight="bold")
    arrow(ax, 0.525, 0.50, 0.565, 0.50)
    # transformer
    ax.add_patch(FancyBboxPatch((0.57, 0.26), 0.155, 0.48,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                fc="#fdf6ec", ec=ORANGE, lw=1.0, zorder=1))
    ax.text(0.6475, 0.685, "Transformer", ha="center", fontsize=7.4, fontweight="bold", color=INK)
    ax.text(0.6475, 0.635, "2 layers, 8 heads", ha="center", fontsize=6.2, color="#555555")
    for i in range(6):
        for j in range(6):
            ax.plot([0.60 + i * 0.018], [0.55 - j * 0.018], "s", ms=2.4,
                    color=ORANGE, alpha=0.25 + 0.1 * ((i + j) % 3), mec="none")
    ax.text(0.6475, 0.435, "all-to-all", ha="center", fontsize=6.2, color=ORANGE,
            fontweight="bold")
    ax.text(0.6475, 0.335, "whole-transcript\nframe context", ha="center", fontsize=5.9,
            color="#555555")
    arrow(ax, 0.73, 0.50, 0.775, 0.615)
    arrow(ax, 0.73, 0.50, 0.775, 0.385)
    box(ax, 0.78, 0.545, 0.215, 0.155, "profile head", "1x1 conv -> softmax over L",
        fc="#eaf1fa", ec=BLUE, fs=7, subfs=5.8)
    box(ax, 0.78, 0.30, 0.215, 0.155, "count head", "masked mean -> MLP -> scalar",
        fc="#fdf0e3", ec=ORANGE, fs=7, subfs=5.8)
    ax.text(0.8875, 0.735, "SHAPE  (B, L)", ha="center", fontsize=6.2, color=BLUE,
            fontweight="bold")
    ax.text(0.8875, 0.245, "DEPTH  (B,)", ha="center", fontsize=6.2, color=ORANGE,
            fontweight="bold")
    ax.text(0.0, 0.075, f"{M['params']} trainable parameters   |   PyTorch   |   "
                        f"loss = multinomial NLL(profile) + {M['count_w']} x MSE(log count)",
            fontsize=6.6, color="#333333")
    ax.text(0.0, 0.015, f"{M['epochs']} epochs, {M['hours']} h on {M['gpu']}; early stop "
                        f"(patience 6) on median per-transcript Pearson; best epoch {M['best_epoch']}",
            fontsize=6.6, color="#666666")

    # ---- (c) real profile -------------------------------------------------
    lbl(0.015, 0.415, "c")
    ax = fig.add_axes([0.07, 0.245, 0.55, 0.145])
    profile_axes(ax, obs, pred, cds, r, tx, c_obs=GREY, c_pred=BLUE)
    ax.set_xlabel("position in transcript (nt)", fontsize=6.5, labelpad=1)
    ax.text(0.015, 0.90, "observed Ribo-seq", transform=ax.transAxes, fontsize=6, color="#6b7075")
    ax.text(0.015, 0.75, "predicted", transform=ax.transAxes, fontsize=6, color=BLUE,
            fontweight="bold")
    cs, ce = cds
    z0 = cs + int(0.28 * (ce - cs)); z1 = z0 + 90
    ax2 = fig.add_axes([0.685, 0.245, 0.265, 0.145])
    profile_axes(ax2, obs, pred, cds, r, tx, c_obs=GREY, c_pred=BLUE, title=False, zoom=(z0, z1))
    ax2.set_title("90 nt zoom: 3-nt periodicity", fontsize=6.6, pad=3)
    ax2.set_xlabel("position (nt)", fontsize=6.5, labelpad=1)

    # ---- (d) evals --------------------------------------------------------
    lbl(0.015, 0.205, "d")
    ax = fig.add_axes([0.07, 0.055, 0.40, 0.125])
    names = [e[0] for e in EVALS[:4]]
    vals = [e[3] for e in EVALS[:4]]
    cols = [BLUE, BLUE, "#6E8B74", "#6E8B74"]
    b = ax.barh(range(len(names))[::-1], vals, color=cols, height=0.62, alpha=0.9)
    for i, (rect, v) in enumerate(zip(b, vals)):
        ax.text(v + 0.012, rect.get_y() + rect.get_height() / 2, f"{v:.3f}",
                va="center", fontsize=6.2, fontweight="bold", color="#333333")
    ax.set_yticks(range(len(names))[::-1])
    ax.set_yticklabels([f"{e[0]}\n{e[2]}" for e in EVALS[:4]], fontsize=6)
    ax.set_xlim(0, 1.12); ax.set_xticks([0, 0.5, 1.0])
    ax.tick_params(labelsize=6, length=2, pad=1)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_title("Held-out performance", fontsize=7, pad=3, loc="left")

    ax = fig.add_axes([0.575, 0.055, 0.375, 0.125])
    xn = [x[0] for x in XSPECIES]; xv = [x[1] for x in XSPECIES]
    ax.bar(range(len(xn)), xv, color="#6E8B74", width=0.62, alpha=0.9)
    ax.set_ylim(0.9, 1.0)
    ax.set_xticks(range(len(xn)))
    ax.set_xticklabels(xn, fontsize=5.6, rotation=28, ha="right")
    ax.set_ylabel("CDS F1", fontsize=6.2, labelpad=2)
    ax.tick_params(labelsize=6, length=2, pad=1)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_title("Cross-species transfer (human-trained, mouse held-outs)",
                 fontsize=7, pad=3, loc="left")
    save(fig, "style1_panel")


# ================================================================== STYLE 2
def style2(ex):
    """RiNALMo-style wide horizontal pipeline: data -> model tower -> heads -> eval fan-out."""
    TEAL, SLATE, AMBER, PLUM = "#12726B", "#33415C", "#D98A0B", "#7D5A7B"
    BG, SOFT = "#fbfcfc", "#eef4f3"
    tx, obs, pred, cds, r = ex
    fig = plt.figure(figsize=(13.0, 5.6), facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ax.text(0.012, 0.955, "Predicting ribosome P-site density from transcript sequence",
            fontsize=13, fontweight="bold", color=SLATE)
    ax.text(0.012, 0.905, "dilated convolutional body + transformer, dual-head output",
            fontsize=8.5, color="#6b7684")

    for x, t in [(0.135, "TRAINING DATA"), (0.375, "MODEL"), (0.635, "OUTPUT"),
                 (0.855, "EVALUATION")]:
        ax.text(x, 0.845, "   ".join(t), fontsize=7.2, fontweight="bold", color="#8a94a2",
                ha="center")
        ax.plot([x - 0.10, x + 0.10], [0.825, 0.825], color="#d5dde2", lw=1.0)

    # ---- training data ----
    ax.add_patch(FancyBboxPatch((0.035, 0.30), 0.20, 0.49,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc="white", ec="#d5dde2", lw=1.0))
    ax.text(0.135, 0.745, "8 human tissues", ha="center", fontsize=9, fontweight="bold",
            color=SLATE)
    ax.text(0.135, 0.705, "Chothani Ribo-seq + RNA-seq", ha="center", fontsize=6.6,
            color="#6b7684")
    for i, t in enumerate(TRAIN_TISSUES):
        cx = 0.055 + (i % 2) * 0.098
        cy = 0.645 - (i // 2) * 0.052
        ax.add_patch(FancyBboxPatch((cx, cy - 0.019), 0.088, 0.038,
                                    boxstyle="round,pad=0,rounding_size=0.008",
                                    fc=SOFT, ec=TEAL, lw=0.7))
        ax.text(cx + 0.044, cy, t, ha="center", va="center", fontsize=6.2, color=SLATE)
    ax.add_patch(FancyBboxPatch((0.055, 0.437 - 0.019), 0.186, 0.038,
                                boxstyle="round,pad=0,rounding_size=0.008",
                                fc="#fdeeee", ec="#C0504D", lw=1.0, ls="--"))
    ax.text(0.148, 0.437, "Hepatocytes  =  HELD OUT", ha="center", va="center",
            fontsize=6.4, color="#9c3f3c", fontweight="bold")
    ax.text(0.135, 0.385, "84,472 transcripts (union universe)", ha="center", fontsize=6.4,
            color="#6b7684")
    ax.text(0.135, 0.348, "validation genes chromosome-disjoint", ha="center", fontsize=6.4,
            color="#6b7684")

    # inputs strip
    ax.text(0.135, 0.275, "per-nucleotide input, 10 channels", ha="center", fontsize=6.8,
            fontweight="bold", color=SLATE)
    for i, (nm, nch, c) in enumerate([("one-hot seq", 4, TEAL), ("frame track", 5, AMBER),
                                      ("RNA coverage", 1, PLUM)]):
        y = 0.225 - i * 0.048
        ax.add_patch(Rectangle((0.045, y - 0.016), 0.18, 0.032, fc=c, alpha=0.16, ec="none"))
        ax.text(0.052, y, nm, fontsize=6.2, va="center", color=SLATE)
        ax.text(0.219, y, f"{nch}", fontsize=6.2, va="center", ha="right", color=c,
                fontweight="bold")
    arrow(ax, 0.243, 0.50, 0.283, 0.50, c="#9aa6b2", lw=1.6, ms=11)

    # ---- model tower ----
    ax.add_patch(FancyBboxPatch((0.29, 0.135), 0.175, 0.66,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc="white", ec=TEAL, lw=1.4))
    ax.text(0.3775, 0.755, "1x1 conv stem", ha="center", fontsize=7, color="#6b7684")
    ax.add_patch(FancyBboxPatch((0.305, 0.715), 0.145, 0.028,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                fc="#e7f1f0", ec=TEAL, lw=0.7))
    ax.text(0.3775, 0.729, "10 -> 256 channels", ha="center", va="center", fontsize=6, color=SLATE)
    ys = np.linspace(0.655, 0.345, 10)
    for i, (yy, dl) in enumerate(zip(ys, DILATIONS)):
        w = 0.055 + 0.0085 * i
        ax.add_patch(FancyBboxPatch((0.3775 - w / 2, yy - 0.0125), w, 0.025,
                                    boxstyle="round,pad=0,rounding_size=0.005",
                                    fc=TEAL, ec="none", alpha=0.22 + 0.06 * i))
        ax.text(0.3775, yy, f"dilation {dl}", ha="center", va="center", fontsize=5.5,
                color=SLATE)
    ax.text(0.3775, 0.685, "10 residual dilated blocks", ha="center", fontsize=7,
            fontweight="bold", color=SLATE)
    ax.annotate("", xy=(0.452, 0.345), xytext=(0.452, 0.655),
                arrowprops=dict(arrowstyle="<->", color="#b6c2cc", lw=0.9))
    ax.text(0.470, 0.655, "receptive\nfield ~4 kb", fontsize=5.9, va="top", color="#6b7684")
    ax.add_patch(FancyBboxPatch((0.305, 0.225), 0.145, 0.075,
                                boxstyle="round,pad=0,rounding_size=0.008",
                                fc="#fdf4e6", ec=AMBER, lw=1.0))
    ax.text(0.3775, 0.277, "Transformer x2", ha="center", fontsize=7.4, fontweight="bold",
            color=SLATE)
    ax.text(0.3775, 0.247, "8 heads, all-to-all", ha="center", fontsize=6.1, color="#6b7684")
    ax.text(0.3775, 0.185, f"{M['params']} parameters", ha="center", fontsize=6.6,
            fontweight="bold", color=TEAL)
    ax.text(0.3775, 0.157, f"{M['epochs']} epochs, {M['hours']} h, {M['gpu']}", ha="center",
            fontsize=6, color="#6b7684")
    arrow(ax, 0.473, 0.50, 0.513, 0.50, c="#9aa6b2", lw=1.6, ms=11)

    # ---- outputs ----
    ax.text(0.6325, 0.795, "loss = multinomial NLL  +  0.1 x MSE", ha="center", fontsize=6.6,
            fontweight="bold", color=SLATE)
    ax.text(0.6325, 0.750, "the SHAPE: where translation happens", ha="center", fontsize=6.2,
            color=TEAL, style="italic")
    box(ax, 0.52, 0.595, 0.225, 0.115, "PROFILE HEAD", "softmax over positions -> (B, L)",
        fc="#e7f1f0", ec=TEAL, fs=8, subfs=6.2)
    # real held-out output sits directly under the head that produces it
    axp = fig.add_axes([0.523, 0.395, 0.219, 0.135])
    profile_axes(axp, obs, pred, cds, r, tx, c_obs="#c3c9cf", c_pred=TEAL, title=False)
    axp.set_title(f"real held-out prediction, r = {r:.2f}", fontsize=6.4, pad=2, color=SLATE)
    axp.set_xlabel("position (nt)", fontsize=6, labelpad=1)
    box(ax, 0.52, 0.195, 0.225, 0.115, "COUNT HEAD", "log total P-sites -> (B,)",
        fc="#fdf4e6", ec=AMBER, fs=8, subfs=6.2)
    ax.text(0.6325, 0.155, "the DEPTH: how much, dataset-specific", ha="center", fontsize=6.2,
            color=AMBER, style="italic")
    arrow(ax, 0.752, 0.50, 0.792, 0.50, c="#9aa6b2", lw=1.6, ms=11)

    # ---- eval fan ----
    ylab = np.linspace(0.735, 0.185, len(EVALS))
    for (grp, sub, met, val), yy in zip(EVALS, ylab):
        ax.add_patch(FancyBboxPatch((0.80, yy - 0.048), 0.185, 0.096,
                                    boxstyle="round,pad=0,rounding_size=0.01",
                                    fc="white", ec="#d5dde2", lw=0.9))
        ax.text(0.812, yy + 0.022, grp, fontsize=7.2, fontweight="bold", color=SLATE, va="center")
        ax.text(0.812, yy - 0.004, sub.replace("\n", " "), fontsize=5.9, color="#6b7684",
                va="center")
        if val is not None:
            ax.text(0.812, yy - 0.030, f"{met} {val:.3f}", fontsize=6.4, color=TEAL,
                    fontweight="bold", va="center")
        else:
            ax.text(0.812, yy - 0.030, "440 novel peptides at 1% class FDR", fontsize=5.8,
                    color=PLUM, fontweight="bold", va="center")
        arrow(ax, 0.792, 0.50, 0.798, yy, c="#c9d3da", lw=0.9, ms=7)
    save(fig, "style2_flow")


# ================================================================== STYLE 3
def style3(ex):
    """Orthrus-style A/B/C tiers: data construction, training pipeline, evaluation breadth."""
    SLATE, SAGE, CLAY, STONE = "#40506B", "#6E8B74", "#B5714F", "#8C8681"
    PAPER, PANEL = "#ffffff", "#f6f7f5"
    tx, obs, pred, cds, r = ex
    fig = plt.figure(figsize=(9.0, 9.2), facecolor=PAPER)

    def tier(y, h, letter, title):
        a = fig.add_axes([0.055, y, 0.90, h]); a.axis("off")
        a.set_xlim(0, 1); a.set_ylim(0, 1)
        a.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=0.008",
                                   fc=PANEL, ec="#e3e5e0", lw=1.0, zorder=0))
        fig.text(0.022, y + h - 0.004, letter, fontsize=12, fontweight="bold", color=SLATE,
                 va="top")
        a.text(0.018, 0.93, title, fontsize=9.2, fontweight="bold", color=SLATE, va="center")
        return a

    fig.text(0.055, 0.982, "Ribo-seq signal model: data, training, evaluation",
             fontsize=12.5, fontweight="bold", color=SLATE, va="top")
    fig.text(0.055, 0.960, "one-hot attention architecture (orf_v2_attn), union universe",
             fontsize=8, color=STONE, va="top")

    # ---------- A: data construction ----------
    a = tier(0.655, 0.275, "A", "Data construction and leave-one-tissue-out split")
    for i, t in enumerate(TRAIN_TISSUES):
        cx = 0.035 + i * 0.088
        a.add_patch(FancyBboxPatch((cx, 0.60), 0.078, 0.16,
                                   boxstyle="round,pad=0,rounding_size=0.02",
                                   fc=SAGE, ec="none", alpha=0.30))
        a.text(cx + 0.039, 0.68, t, ha="center", va="center", fontsize=6.1, color=SLATE)
    a.add_patch(FancyBboxPatch((0.035 + 7 * 0.088, 0.60), 0.078, 0.16,
                               boxstyle="round,pad=0,rounding_size=0.02",
                               fc="none", ec=CLAY, lw=1.3, ls="--"))
    a.text(0.035 + 7 * 0.088 + 0.039, 0.68, "Hepato-\ncytes", ha="center", va="center",
           fontsize=6.1, color=CLAY, fontweight="bold")
    a.text(0.035, 0.545, "7 tissues -> TRAIN", fontsize=6.8, color=SAGE, fontweight="bold")
    a.text(0.655, 0.545, "1 tissue -> HELD OUT (never seen)", fontsize=6.8, color=CLAY,
           fontweight="bold")
    a.plot([0.03, 0.63], [0.575, 0.575], color=SAGE, lw=1.6)
    a.plot([0.648, 0.965], [0.575, 0.575], color=CLAY, lw=1.6)
    for i, (nm, sub) in enumerate([
            ("Ribo-seq", "per-nt P-site counts = the LABEL"),
            ("RNA-seq", "per-nt coverage = an INPUT"),
            ("GENCODE v49", "sequence + reading-frame track = INPUTS")]):
        y = 0.40 - i * 0.115
        a.add_patch(Rectangle((0.035, y - 0.035), 0.012, 0.07,
                              fc=[CLAY, SAGE, SLATE][i], ec="none"))
        a.text(0.058, y, nm, fontsize=7.4, fontweight="bold", va="center", color=SLATE)
        a.text(0.185, y, sub, fontsize=6.6, va="center", color="#6b6b6b")
    a.text(0.60, 0.40, "84,472 transcripts", fontsize=7.4, fontweight="bold", color=SLATE,
           va="center")
    a.text(0.60, 0.34, "protein-coding + lncRNA, up to 10 kb, chrM excluded",
           fontsize=6.4, color="#6b6b6b", va="center")
    a.text(0.60, 0.245, "validation genes are chromosome-disjoint from training",
           fontsize=6.4, color="#6b6b6b", va="center")
    a.text(0.60, 0.185, "so early stopping cannot leak into the held-out result",
           fontsize=6.4, color=CLAY, va="center", style="italic")

    # ---------- B: training pipeline ----------
    b = tier(0.325, 0.305, "B", "Training pipeline")
    b.text(0.018, 0.83, "10-track per-nucleotide encoding", fontsize=7, color="#6b6b6b")
    for i, (nm, nch, c) in enumerate([("one-hot A C G U", 4, SLATE),
                                      ("reading frame + start/stop", 5, CLAY),
                                      ("log RNA-seq coverage", 1, SAGE)]):
        y = 0.70 - i * 0.10
        b.add_patch(Rectangle((0.02, y - 0.033), 0.20, 0.066, fc=c, alpha=0.16, ec="none"))
        b.text(0.03, y, nm, fontsize=6.3, va="center", color=SLATE)
        b.text(0.214, y, f"{nch}", fontsize=6.3, va="center", ha="right", color=c,
               fontweight="bold")
    arrow(b, 0.232, 0.60, 0.263, 0.60, c="#a8aca6", lw=1.3, ms=9)
    b.add_patch(FancyBboxPatch((0.27, 0.30), 0.20, 0.60,
                               boxstyle="round,pad=0,rounding_size=0.02",
                               fc="white", ec=SLATE, lw=1.1))
    b.text(0.37, 0.855, "encoder", ha="center", fontsize=7.6, fontweight="bold", color=SLATE)
    dilated_motif(b, 0.295, 0.545, 0.15, 0.24, rows=5, n=13, c=SLATE, dot=2.0, lw=0.5)
    b.text(0.37, 0.505, "10 dilated blocks, 1 to 512", ha="center", fontsize=6.0,
           color="#6b6b6b")
    b.add_patch(FancyBboxPatch((0.295, 0.365), 0.15, 0.10,
                               boxstyle="round,pad=0,rounding_size=0.012",
                               fc=CLAY, ec="none", alpha=0.22))
    b.text(0.37, 0.415, "Transformer x2, 8 heads", ha="center", va="center", fontsize=6.2,
           color=SLATE)
    b.text(0.37, 0.325, f"{M['params']} params", ha="center", fontsize=6.2, color=SLATE,
           fontweight="bold")
    arrow(b, 0.478, 0.60, 0.508, 0.68, c="#a8aca6", lw=1.3, ms=9)
    arrow(b, 0.478, 0.60, 0.508, 0.45, c="#a8aca6", lw=1.3, ms=9)
    box(b, 0.515, 0.62, 0.19, 0.135, "profile head", "softmax over L", fc="white", ec=SLATE,
        fs=7, subfs=5.9)
    box(b, 0.515, 0.385, 0.19, 0.135, "count head", "scalar log total", fc="white", ec=CLAY,
        fs=7, subfs=5.9)
    b.add_patch(FancyBboxPatch((0.735, 0.385), 0.245, 0.37,
                               boxstyle="round,pad=0,rounding_size=0.015",
                               fc="white", ec="#dcdfda", lw=1.0))
    b.text(0.8575, 0.705, "objective", ha="center", fontsize=7.2, fontweight="bold", color=SLATE)
    b.text(0.8575, 0.635, "multinomial NLL (profile)", ha="center", fontsize=6.3, color=SLATE)
    b.text(0.8575, 0.585, "+  0.1 x MSE (log count)", ha="center", fontsize=6.3, color=CLAY)
    b.plot([0.762, 0.953], [0.545, 0.545], color="#e3e5e0", lw=1.0)
    b.text(0.8575, 0.495, "select on median per-transcript", ha="center", fontsize=6.0,
           color="#6b6b6b")
    b.text(0.8575, 0.452, "Pearson (validation)", ha="center", fontsize=6.0, color="#6b6b6b")
    b.text(0.8575, 0.410, f"best epoch {M['best_epoch']}  |  {M['hours']} h  |  {M['gpu']}",
           ha="center", fontsize=5.9, color=STONE)
    # real held-out output, tucked into the free space under the objective box (tier B)
    axp = fig.add_axes([0.615, 0.352, 0.265, 0.060])
    profile_axes(axp, obs, pred, cds, r, tx, c_obs="#cfd2cc", c_pred=SLATE, title=False)
    axp.set_title(f"real held-out prediction, r = {r:.2f}", fontsize=5.8, pad=2, color=SLATE)
    axp.set_xticks([])

    # ---------- C: evaluation ----------
    c = tier(0.035, 0.265, "C", "Evaluation: four axes of generalization")
    cols = [
        ("In-distribution", ["held-out tissue", "profile r = 0.659", "CDS F1 = 0.953"], SLATE),
        ("Cross-study", ["Ruiz-Orera human", "CDS F1 = 0.931", "beats obs ceiling"], SAGE),
        ("Cross-species", ["5 mouse datasets", "CDS F1 0.960 to 0.995", "human-trained only"], CLAY),
        ("Proteogenomics", ["12 macrophage pops", "440 novel peptides", "1% class FDR"], STONE),
    ]
    for i, (title, rows_, col) in enumerate(cols):
        x = 0.02 + i * 0.246
        c.add_patch(FancyBboxPatch((x, 0.30), 0.226, 0.53,
                                   boxstyle="round,pad=0,rounding_size=0.015",
                                   fc="white", ec="#dcdfda", lw=1.0))
        c.add_patch(Rectangle((x, 0.795), 0.226, 0.035, fc=col, ec="none", alpha=0.75))
        c.text(x + 0.113, 0.745, title, ha="center", fontsize=7.6, fontweight="bold", color=col)
        for j, t in enumerate(rows_):
            c.text(x + 0.113, 0.655 - j * 0.095, t, ha="center", fontsize=6.3,
                   color="#4a4a4a" if j else "#6b6b6b",
                   fontweight="bold" if j == 1 else "normal")
    c.text(0.02, 0.185, "Trained once on human tissue Ribo-seq; every axis above is data the model "
                        "never saw during training.", fontsize=6.8, color="#6b6b6b")
    c.text(0.02, 0.105, "No ribosome profiling is required at inference: the model predicts the "
                        "profile from sequence + RNA-seq alone.", fontsize=6.8, color=CLAY,
           style="italic")
    save(fig, "style3_tiers")


def main():
    print("loading a real held-out example from the union model dump ...")
    ex = load_example()
    if ex is None:
        raise SystemExit("no suitable example transcript found")
    print(f"  example: {ex[0]}  L={len(ex[1])}  obs_total={ex[1].sum():.0f}  r={ex[4]:.3f}")
    (HERE).mkdir(parents=True, exist_ok=True)
    json.dump({"example_tx": ex[0], "length": int(len(ex[1])),
               "obs_total": float(ex[1].sum()), "pearson": float(ex[4]),
               "cds": [int(ex[3][0]), int(ex[3][1])], "source_npz": str(NPZ)},
              open(HERE / "arch_attn_values.json", "w"), indent=1)
    print("style 1 (BPNet-style panel) ..."); style1(ex)
    print("style 2 (RiNALMo-style flow) ..."); style2(ex)
    print("style 3 (Orthrus-style tiers) ..."); style3(ex)
    print(f"\nall three styles -> {HERE}")


if __name__ == "__main__":
    main()
