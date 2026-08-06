#!/usr/bin/env python3
"""Model-architecture panel ONLY (no data / training / evaluation panels), for the one-hot attention model
(orf_v2_attn). Three presentation styles to pick from, each a standalone figure sized for a manuscript
architecture panel. Everything drawn here is read off scripts/model.py (RiboSignalModel +
ResidualDilatedBlock) and the union run's args.json -- layer types, kernel sizes, dilation schedule,
normalisation, head structure and tensor shapes are all literal, so the panel is reimplementable.

  arch_style1 : horizontal journal flow, residual-block internals as an inset  (BPNet-flavoured)
  arch_style2 : vertical tower with the dilation pyramid, exploded block detail (RiNALMo-flavoured)
  arch_style3 : modular blocks with explicit tensor shapes between stages       (Orthrus-flavoured)

cas12a env (matplotlib). Run: python make_arch_panel.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

HERE = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
            "riboseq_signal_model/figures/arch_attn")
DIL = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
PARAMS = "5,071,106"


def box(ax, x, y, w, h, title, sub=None, fc="#ffffff", ec="#333333", tc="#111111",
        fs=8, subfs=6.3, lw=1.0, rad=0.012, bold=True, dy=0.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rad}",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2 + (0.055 * h if sub else 0) + dy, title, ha="center",
            va="center", fontsize=fs, color=tc, fontweight="bold" if bold else "normal", zorder=3)
    if sub:
        ax.text(x + w / 2, y + h / 2 - 0.20 * h + dy, sub, ha="center", va="center",
                fontsize=subfs, color=tc, alpha=0.78, zorder=3)


def arrow(ax, x1, y1, x2, y2, c="#666666", lw=1.2, ms=8, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, color=c, lw=lw,
                                 mutation_scale=ms, shrinkA=0, shrinkB=0, zorder=4))


def shape_tag(ax, x, y, txt, c="#7a8590", fs=6.2, ha="center"):
    ax.text(x, y, txt, ha=ha, va="center", fontsize=fs, color=c, family="monospace", zorder=5)


def dilated_motif(ax, x, y, w, h, rows=5, n=17, c="#2C6FBB", dot=2.3, lw=0.6):
    xs = np.linspace(x, x + w, n)
    ys = np.linspace(y, y + h, rows)
    for r in range(rows):
        ax.plot(xs, [ys[r]] * n, "o", ms=dot, color=c, alpha=0.35 + 0.13 * r, mec="none", zorder=3)
    for r in range(rows - 1):
        d = 2 ** r
        for i in range(n):
            if i - d >= 0:
                ax.plot([xs[i], xs[i - d]], [ys[r], ys[r + 1]], "-", color=c, lw=lw,
                        alpha=0.5, zorder=2)


def save(fig, name):
    HERE.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"{name}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


# ================================================================= STYLE 1
def style1():
    """Horizontal journal flow; residual-block internals broken out as an inset below."""
    INK, GREY, BLUE, ORANGE = "#1a1a1a", "#8a9199", "#2C6FBB", "#E08214"
    fig = plt.figure(figsize=(9.2, 4.5), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ax.text(0.012, 0.955, "Dual-head dilated CNN + Transformer", fontsize=11,
            fontweight="bold", color=INK)
    ax.text(0.012, 0.905, f"per-nucleotide Ribo-seq P-site prediction   |   PyTorch   |   "
                          f"{PARAMS} parameters", fontsize=7.2, color="#666666")

    yb, hb = 0.50, 0.30                                   # main row geometry
    # input
    box(ax, 0.012, yb, 0.105, hb, "input", "one-hot 4 + ORF 5\n+ coverage 1", fc="#f4f6f8",
        ec="#9aa3ab", fs=7.4, subfs=5.9)
    shape_tag(ax, 0.064, yb - 0.055, "(B, L, 10)")
    arrow(ax, 0.120, yb + hb / 2, 0.152, yb + hb / 2)
    # stem
    box(ax, 0.155, yb, 0.088, hb, "Conv1d", "kernel 1\n10 -> 256", fc="#f4f6f8", ec="#9aa3ab",
        fs=7.4, subfs=5.9)
    shape_tag(ax, 0.199, yb - 0.055, "(B, 256, L)")
    arrow(ax, 0.246, yb + hb / 2, 0.278, yb + hb / 2)
    # dilated body
    ax.add_patch(FancyBboxPatch((0.281, yb - 0.035), 0.30, hb + 0.07,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc="#f5f9fd", ec=BLUE, lw=1.2, zorder=1))
    ax.text(0.431, yb + hb + 0.005, "10 x residual dilated block", ha="center", fontsize=8,
            fontweight="bold", color=INK)
    dilated_motif(ax, 0.303, yb + 0.045, 0.256, 0.175, rows=5, n=17, c=BLUE)
    ax.text(0.431, yb + 0.012, "dilation 1, 2, 4, 8 ... 512   |   receptive field ~4 kb",
            ha="center", fontsize=6.4, color="#4a5560", style="italic")
    shape_tag(ax, 0.431, yb - 0.055, "(B, 256, L)")
    arrow(ax, 0.584, yb + hb / 2, 0.616, yb + hb / 2)
    # transformer
    ax.add_patch(FancyBboxPatch((0.619, yb - 0.035), 0.135, hb + 0.07,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc="#fdf7ee", ec=ORANGE, lw=1.2, zorder=1))
    ax.text(0.6865, yb + hb + 0.005, "2 x Transformer", ha="center", fontsize=8,
            fontweight="bold", color=INK)
    for i in range(7):
        for j in range(5):
            ax.plot([0.638 + i * 0.0165], [yb + 0.20 - j * 0.032], "s", ms=3.0, color=ORANGE,
                    alpha=0.22 + 0.10 * ((i * j) % 4), mec="none")
    ax.text(0.6865, yb + 0.012, "8 heads, pre-norm\nd_ff 512, all-to-all", ha="center",
            fontsize=6.2, color="#4a5560")
    shape_tag(ax, 0.6865, yb - 0.055, "(B, 256, L)")
    arrow(ax, 0.757, yb + hb / 2, 0.788, yb + hb * 0.80)
    arrow(ax, 0.757, yb + hb / 2, 0.788, yb + hb * 0.20)
    # heads
    box(ax, 0.791, yb + 0.155, 0.197, 0.145, "profile head", "Conv1d(256->1, k=1)\nsoftmax over L",
        fc="#eaf2fb", ec=BLUE, fs=7.4, subfs=5.9)
    box(ax, 0.791, yb - 0.035, 0.197, 0.145, "count head", "masked mean + log cov\nMLP 257->256->1",
        fc="#fdf1e4", ec=ORANGE, fs=7.4, subfs=5.9)
    shape_tag(ax, 0.889, yb + 0.325, "(B, L)   the SHAPE", c=BLUE, fs=6.4)
    shape_tag(ax, 0.889, yb - 0.075, "(B,)   the DEPTH", c=ORANGE, fs=6.4)

    # ---- inset: what one residual block contains ----
    ax.plot([0.281, 0.281], [yb - 0.045, 0.345], color=BLUE, lw=0.7, ls=(0, (2, 2)))
    ax.plot([0.581, 0.581], [yb - 0.045, 0.345], color=BLUE, lw=0.7, ls=(0, (2, 2)))
    ax.add_patch(FancyBboxPatch((0.281, 0.115), 0.30, 0.225,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec=BLUE, lw=0.9, zorder=1))
    ax.text(0.296, 0.305, "one block, dilation d", fontsize=6.8, fontweight="bold", color=BLUE)
    steps = ["Conv1d\nk=3, dil d", "GroupNorm\n8 groups", "GELU", "Conv1d\nk=3, dil d",
             "GroupNorm"]
    xs0, wstep = 0.293, 0.0535
    for i, s in enumerate(steps):
        x = xs0 + i * wstep
        ax.add_patch(FancyBboxPatch((x, 0.175), 0.0455, 0.075,
                                    boxstyle="round,pad=0,rounding_size=0.006",
                                    fc="#f5f9fd", ec="#a9c4e0", lw=0.7))
        ax.text(x + 0.0228, 0.2125, s, ha="center", va="center", fontsize=5.2, color="#33414d")
        if i < len(steps) - 1:
            arrow(ax, x + 0.0455, 0.2125, x + wstep, 0.2125, c="#a9c4e0", lw=0.8, ms=5)
    ax.annotate("", xy=(0.5645, 0.155), xytext=(0.2925, 0.155),
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8,
                                connectionstyle="arc3,rad=0.0"))
    arrow(ax, 0.5645, 0.155, 0.5645, 0.172, c=GREY, lw=0.8, ms=6)
    ax.text(0.428, 0.138, "residual add, then GELU", ha="center", fontsize=5.6, color="#5b6670")

    ax.text(0.012, 0.055, "loss = multinomial NLL(profile)  +  0.1 x MSE(log count).  The profile head is "
                          "scale-free (shape only); depth is isolated in the count head.",
            fontsize=6.5, color="#555555")
    save(fig, "arch_style1")


# ================================================================= STYLE 2
def style2():
    """Vertical tower: dilation pyramid centre, exploded block detail right, heads at the bottom."""
    TEAL, SLATE, AMBER = "#12726B", "#33415C", "#D98A0B"
    fig = plt.figure(figsize=(6.6, 8.6), facecolor="#fbfcfc")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ax.text(0.5, 0.972, "Dual-head dilated CNN + Transformer", fontsize=12.5,
            fontweight="bold", color=SLATE, ha="center")
    ax.text(0.5, 0.947, f"{PARAMS} trainable parameters   |   PyTorch", fontsize=7.6,
            color="#6b7684", ha="center")

    cx, cw = 0.30, 0.34                                    # tower centre / width
    # input
    box(ax, cx - cw / 2, 0.875, cw, 0.052, "per-nucleotide input", None, fc="white",
        ec="#c7d0d6", fs=8)
    for i, (nm, nch, c) in enumerate([("one-hot A C G U", 4, TEAL),
                                      ("reading frame + start/stop", 5, AMBER),
                                      ("log RNA-seq coverage", 1, "#7D5A7B")]):
        y = 0.845 - i * 0.030
        ax.add_patch(Rectangle((cx - cw / 2, y - 0.012), cw, 0.024, fc=c, alpha=0.15, ec="none"))
        ax.text(cx - cw / 2 + 0.012, y, nm, fontsize=6.1, va="center", color=SLATE)
        ax.text(cx + cw / 2 - 0.012, y, str(nch), fontsize=6.1, va="center", ha="right",
                color=c, fontweight="bold")
    shape_tag(ax, cx + cw / 2 + 0.045, 0.845, "(B, L, 10)", c="#8f9aa3", ha="left")
    arrow(ax, cx, 0.748, cx, 0.722, c="#9aa6b2", lw=1.4, ms=10)
    box(ax, cx - cw / 2, 0.678, cw, 0.042, "Conv1d  kernel 1,  10 -> 256", None, fc="white",
        ec=TEAL, fs=7.2)
    arrow(ax, cx, 0.676, cx, 0.652, c="#9aa6b2", lw=1.4, ms=10)

    # dilation pyramid
    ax.text(cx, 0.632, "10 x residual dilated block", ha="center", fontsize=8.4,
            fontweight="bold", color=SLATE)
    ys = np.linspace(0.598, 0.335, 10)
    for i, (yy, d) in enumerate(zip(ys, DIL)):
        w = 0.115 + 0.0225 * i
        ax.add_patch(FancyBboxPatch((cx - w / 2, yy - 0.0115), w, 0.023,
                                    boxstyle="round,pad=0,rounding_size=0.005",
                                    fc=TEAL, ec="none", alpha=0.20 + 0.062 * i))
        ax.text(cx, yy, f"dilation {d}", ha="center", va="center", fontsize=5.8, color=SLATE)
    ax.annotate("", xy=(cx + 0.185, 0.335), xytext=(cx + 0.185, 0.598),
                arrowprops=dict(arrowstyle="<->", color="#b6c2cc", lw=1.0))
    ax.text(cx + 0.196, 0.467, "receptive field\n~4 kb, O(L)", fontsize=6.2, va="center",
            color="#6b7684")
    shape_tag(ax, cx - 0.20, 0.467, "(B, 256, L)", c="#8f9aa3")
    arrow(ax, cx, 0.325, cx, 0.300, c="#9aa6b2", lw=1.4, ms=10)

    # transformer
    ax.add_patch(FancyBboxPatch((cx - cw / 2, 0.228), cw, 0.070,
                                boxstyle="round,pad=0,rounding_size=0.010",
                                fc="#fdf4e6", ec=AMBER, lw=1.2))
    ax.text(cx, 0.279, "2 x Transformer encoder", ha="center", fontsize=8,
            fontweight="bold", color=SLATE)
    ax.text(cx, 0.252, "8 heads, pre-norm, d_ff 512", ha="center", fontsize=6.3, color="#6b7684")
    ax.text(cx, 0.236, "all-to-all: whole-transcript frame context", ha="center", fontsize=5.9,
            color=AMBER, style="italic")
    arrow(ax, cx - 0.06, 0.226, cx - 0.115, 0.185, c="#9aa6b2", lw=1.3, ms=9)
    arrow(ax, cx + 0.06, 0.226, cx + 0.115, 0.185, c="#9aa6b2", lw=1.3, ms=9)

    # heads
    box(ax, cx - 0.245, 0.088, 0.235, 0.092, "profile head", "Conv1d(256->1, k=1)\nsoftmax over L",
        fc="#e7f1f0", ec=TEAL, fs=7.6, subfs=5.9)
    box(ax, cx + 0.010, 0.088, 0.235, 0.092, "count head", "masked mean + log cov\nMLP 257->256->1",
        fc="#fdf4e6", ec=AMBER, fs=7.6, subfs=5.9)
    shape_tag(ax, cx - 0.1275, 0.068, "(B, L)  SHAPE", c=TEAL, fs=6.4)
    shape_tag(ax, cx + 0.1275, 0.068, "(B,)  DEPTH", c=AMBER, fs=6.4)
    ax.text(cx, 0.032, "loss = multinomial NLL(profile)  +  0.1 x MSE(log count)",
            ha="center", fontsize=6.6, color=SLATE, fontweight="bold")

    # exploded block detail (right column)
    bx = 0.755
    ax.add_patch(FancyBboxPatch((bx - 0.115, 0.335), 0.23, 0.30,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec=TEAL, lw=1.0))
    ax.text(bx, 0.612, "inside one block", ha="center", fontsize=7.4, fontweight="bold",
            color=TEAL)
    steps = ["Conv1d  k=3, dilation d", "GroupNorm (8)", "GELU", "Dropout 0.1",
             "Conv1d  k=3, dilation d", "GroupNorm (8)"]
    for i, s in enumerate(steps):
        y = 0.578 - i * 0.036
        ax.add_patch(FancyBboxPatch((bx - 0.100, y - 0.0135), 0.20, 0.027,
                                    boxstyle="round,pad=0,rounding_size=0.005",
                                    fc="#f2f8f7", ec="#bcd6d3", lw=0.7))
        ax.text(bx, y, s, ha="center", va="center", fontsize=5.7, color=SLATE)
        if i < len(steps) - 1:
            arrow(ax, bx, y - 0.0145, bx, y - 0.0215, c="#bcd6d3", lw=0.7, ms=5)
    ax.add_patch(FancyArrowPatch((bx - 0.104, 0.585), (bx - 0.104, 0.368),
                                 connectionstyle="arc3,rad=0.55", arrowstyle="-",
                                 color="#8fb3ae", lw=0.9))
    arrow(ax, bx - 0.104, 0.368, bx - 0.035, 0.362, c="#8fb3ae", lw=0.9, ms=6)
    ax.text(bx + 0.005, 0.356, "residual add -> GELU", fontsize=5.7, color="#5c6f6d")
    ax.plot([cx + 0.185, bx - 0.115], [0.52, 0.52], color="#d5dde2", lw=0.8, ls=(0, (3, 2)))
    save(fig, "arch_style2")


# ================================================================= STYLE 3
def style3():
    """Modular stages with explicit tensor shapes on every edge; block internals as a sub-panel."""
    SLATE, SAGE, CLAY, STONE = "#40506B", "#6E8B74", "#B5714F", "#8C8681"
    fig = plt.figure(figsize=(9.6, 5.2), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ax.text(0.012, 0.955, "Architecture", fontsize=12, fontweight="bold", color=SLATE)
    ax.text(0.012, 0.912, f"dual-head dilated CNN + Transformer  .  {PARAMS} parameters  .  PyTorch",
            fontsize=7.4, color=STONE)

    ytop, hgt = 0.545, 0.255
    stages = [
        (0.012, 0.148, "INPUT", "one-hot 4  +  ORF track 5\n+  log coverage 1", "(B, L, 10)", SLATE),
        (0.196, 0.128, "STEM", "Conv1d\nkernel 1,  10 -> 256", "(B, 256, L)", SLATE),
        (0.360, 0.230, "BODY", None, "(B, 256, L)", SAGE),
        (0.626, 0.150, "MIXER", None, "(B, 256, L)", CLAY),
    ]
    for x, w, tag, sub, shp, col in stages:
        ax.add_patch(FancyBboxPatch((x, ytop), w, hgt, boxstyle="round,pad=0,rounding_size=0.013",
                                    fc="white", ec=col, lw=1.2, zorder=2))
        ax.add_patch(Rectangle((x, ytop + hgt - 0.030), w, 0.030, fc=col, ec="none",
                               alpha=0.80, zorder=3))
        ax.text(x + w / 2, ytop + hgt - 0.015, tag, ha="center", va="center", fontsize=6.8,
                color="white", fontweight="bold", zorder=4)
        if sub:
            ax.text(x + w / 2, ytop + hgt * 0.40, sub, ha="center", va="center", fontsize=6.6,
                    color="#3f3f3f", zorder=4)
        shape_tag(ax, x + w / 2, ytop - 0.042, shp, c=STONE)
    for x1, x2 in [(0.160, 0.194), (0.324, 0.358), (0.590, 0.624)]:
        arrow(ax, x1, ytop + hgt / 2, x2, ytop + hgt / 2, c="#b9bdb6", lw=1.2, ms=8)
    # BODY: graphic on top, labels beneath it (no overlap)
    dilated_motif(ax, 0.378, ytop + 0.108, 0.194, 0.098, rows=5, n=15, c=SAGE, dot=2.1, lw=0.55)
    ax.text(0.475, ytop + 0.078, "10 x residual dilated block,  dilation 1 ... 512", ha="center",
            fontsize=6.4, color="#3f3f3f", zorder=4)
    ax.text(0.475, ytop + 0.040, "receptive field ~4 kb, linear in L", ha="center", fontsize=5.9,
            color="#5f6b60", style="italic", zorder=4)
    # MIXER: same treatment
    for i in range(6):
        for j in range(3):
            ax.plot([0.652 + i * 0.016], [ytop + 0.190 - j * 0.026], "s", ms=2.8, color=CLAY,
                    alpha=0.20 + 0.10 * ((i + j) % 3), mec="none")
    ax.text(0.701, ytop + 0.078, "2 x Transformer, 8 heads", ha="center", fontsize=6.4,
            color="#3f3f3f", zorder=4)
    ax.text(0.701, ytop + 0.040, "pre-norm, all-to-all", ha="center", fontsize=5.9, color=CLAY,
            style="italic", zorder=4)

    arrow(ax, 0.778, ytop + hgt / 2, 0.812, ytop + hgt * 0.86, c="#b9bdb6", lw=1.2, ms=8)
    arrow(ax, 0.778, ytop + hgt / 2, 0.812, ytop + hgt * 0.14, c="#b9bdb6", lw=1.2, ms=8)
    box(ax, 0.815, ytop + hgt * 0.55, 0.173, 0.115, "profile head",
        "Conv1d(256->1)\nsoftmax over L", fc="white", ec=SLATE, fs=7.2, subfs=5.8)
    box(ax, 0.815, ytop - 0.010, 0.173, 0.115, "count head",
        "masked mean + log cov\nMLP 257->256->1", fc="white", ec=CLAY, fs=7.2, subfs=5.8)
    shape_tag(ax, 0.9015, ytop + hgt * 0.55 + 0.132, "(B, L)   SHAPE", c=SLATE, fs=6.3)
    shape_tag(ax, 0.9015, ytop - 0.042, "(B,)   DEPTH", c=CLAY, fs=6.3)

    # ---- sub-panel: one residual block ----
    ax.add_patch(FancyBboxPatch((0.012, 0.075), 0.560, 0.335,
                                boxstyle="round,pad=0,rounding_size=0.013",
                                fc="#f7f8f6", ec="#dcdfda", lw=1.0, zorder=1))
    ax.text(0.030, 0.368, "one residual dilated block", fontsize=8, fontweight="bold", color=SAGE)
    ax.text(0.030, 0.335, "identical for all 10 blocks; only the dilation d changes",
            fontsize=6.2, color="#6b6b6b")
    steps = ["Conv1d\nk=3, dil d", "GroupNorm\n8 groups", "GELU", "Dropout\n0.1",
             "Conv1d\nk=3, dil d", "GroupNorm\n8 groups"]
    x0, w, gap = 0.032, 0.0765, 0.0115
    for i, s in enumerate(steps):
        x = x0 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 0.185), w, 0.090,
                                    boxstyle="round,pad=0,rounding_size=0.007",
                                    fc="white", ec="#c3ceC3", lw=0.8))
        ax.text(x + w / 2, 0.230, s, ha="center", va="center", fontsize=5.6, color="#3f3f3f")
        if i < len(steps) - 1:
            arrow(ax, x + w, 0.230, x + w + gap, 0.230, c="#c3ceC3", lw=0.8, ms=5)
    # skip connection routed BELOW the block boxes so it never crosses them
    ax.plot([0.032, 0.032], [0.185, 0.163], color=SAGE, lw=0.9)
    ax.plot([0.032, 0.5245], [0.163, 0.163], color=SAGE, lw=0.9)
    arrow(ax, 0.5245, 0.163, 0.5245, 0.184, c=SAGE, lw=0.9, ms=6)
    ax.text(0.278, 0.128, "residual add, then GELU   (padding keeps length L exactly)",
            ha="center", fontsize=6.0, color="#5f6b60")

    ax.text(0.600, 0.330, "objective", fontsize=8, fontweight="bold", color=SLATE)
    ax.text(0.600, 0.288, "multinomial NLL over the profile", fontsize=6.6, color="#3f3f3f")
    ax.text(0.600, 0.252, "+  0.1 x MSE on log total count", fontsize=6.6, color=CLAY)
    ax.plot([0.600, 0.985], [0.222, 0.222], color="#dcdfda", lw=1.0)
    ax.text(0.600, 0.186, "The profile head is scale-free: it predicts only the", fontsize=6.2,
            color="#6b6b6b")
    ax.text(0.600, 0.152, "SHAPE of translation. Sequencing depth is isolated", fontsize=6.2,
            color="#6b6b6b")
    ax.text(0.600, 0.118, "in the count head, so the transferable signal is not", fontsize=6.2,
            color="#6b6b6b")
    ax.text(0.600, 0.084, "entangled with dataset-specific magnitude.", fontsize=6.2,
            color="#6b6b6b")
    save(fig, "arch_style3")


# ================================================================= STYLE 2H
def style2h():
    """CHOSEN DIRECTION: style-2 visual language, laid out horizontally.

    Same palette / bold-block treatment / dilation-pyramid motif / exploded block detail as the vertical
    style2, but the flow runs left to right and the pyramid is rotated so dilation grows ALONG the reading
    direction (bar height encodes depth), which is what makes the receptive-field growth legible in a
    landscape panel. Block internals move from a right-hand column to a bottom strip.
    """
    TEAL, SLATE, AMBER, PLUM = "#12726B", "#33415C", "#D98A0B", "#7D5A7B"
    fig = plt.figure(figsize=(13.6, 5.6), facecolor="#fbfcfc")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ax.text(0.012, 0.955, "Dual-head dilated CNN + Transformer", fontsize=13,
            fontweight="bold", color=SLATE)
    ax.text(0.012, 0.908, f"per-nucleotide Ribo-seq P-site prediction   |   PyTorch   |   "
                          f"{PARAMS} trainable parameters", fontsize=7.8, color="#6b7684")

    ROW = 0.615                                            # flow centre line
    SHP = 0.395                                            # shape-tag baseline

    # ---- input ----
    ax.add_patch(FancyBboxPatch((0.012, 0.487), 0.150, 0.272,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec="#c7d0d6", lw=1.1))
    ax.text(0.087, 0.727, "per-nucleotide input", ha="center", fontsize=8.2,
            fontweight="bold", color=SLATE)
    for i, (nm, nch, c) in enumerate([("one-hot A C G U", 4, TEAL),
                                      ("reading frame + start/stop", 5, AMBER),
                                      ("log RNA-seq coverage", 1, PLUM)]):
        y = 0.667 - i * 0.058
        ax.add_patch(Rectangle((0.024, y - 0.024), 0.126, 0.048, fc=c, alpha=0.15, ec="none"))
        ax.text(0.032, y, nm, fontsize=6.1, va="center", color=SLATE)
        ax.text(0.144, y, str(nch), fontsize=6.4, va="center", ha="right", color=c,
                fontweight="bold")
    shape_tag(ax, 0.087, SHP, "(B, L, 10)", c="#8f9aa3")
    arrow(ax, 0.166, ROW, 0.196, ROW, c="#9aa6b2", lw=1.5, ms=11)

    # ---- stem ----
    ax.add_patch(FancyBboxPatch((0.199, 0.520), 0.108, 0.190,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec=TEAL, lw=1.2))
    ax.text(0.253, 0.652, "Conv1d", ha="center", fontsize=8.4, fontweight="bold", color=SLATE)
    ax.text(0.253, 0.606, "kernel 1", ha="center", fontsize=6.6, color="#6b7684")
    ax.text(0.253, 0.567, "10 -> 256", ha="center", fontsize=6.6, color=TEAL, fontweight="bold")
    shape_tag(ax, 0.253, SHP, "(B, 256, L)", c="#8f9aa3")
    arrow(ax, 0.311, ROW, 0.341, ROW, c="#9aa6b2", lw=1.5, ms=11)

    # ---- body: dilation pyramid rotated to grow left -> right ----
    bx0, bw = 0.344, 0.300
    ax.add_patch(FancyBboxPatch((bx0, 0.455), bw, 0.325,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec=TEAL, lw=1.4))
    ax.text(bx0 + bw / 2, 0.748, "10 x residual dilated block", ha="center", fontsize=8.6,
            fontweight="bold", color=SLATE)
    base, hmax = 0.512, 0.195
    bwid = (bw - 0.036) / len(DIL)
    for i, d in enumerate(DIL):
        h = 0.038 + (hmax - 0.038) * (i / (len(DIL) - 1))
        x = bx0 + 0.018 + i * bwid
        ax.add_patch(FancyBboxPatch((x + 0.0022, base), bwid - 0.0044, h,
                                    boxstyle="round,pad=0,rounding_size=0.004",
                                    fc=TEAL, ec="none", alpha=0.22 + 0.062 * i))
        ax.text(x + bwid / 2, 0.494, str(d), ha="center", va="center", fontsize=5.3,
                color="#5b6b70")
    ax.text(bx0 + bw / 2, 0.472, "dilation", ha="center", fontsize=5.8, color="#8f9aa3")
    ax.annotate("", xy=(bx0 + bw - 0.020, 0.722), xytext=(bx0 + 0.020, 0.722),
                arrowprops=dict(arrowstyle="-|>", color="#b6c2cc", lw=1.0))
    ax.text(bx0 + bw / 2, 0.700, "receptive field ~4 kb, linear in L", ha="center",
            fontsize=6.2, color="#6b7684", style="italic")
    shape_tag(ax, bx0 + bw / 2, SHP, "(B, 256, L)", c="#8f9aa3")
    arrow(ax, bx0 + bw + 0.004, ROW, bx0 + bw + 0.034, ROW, c="#9aa6b2", lw=1.5, ms=11)

    # ---- mixer ----
    mx0, mw = 0.682, 0.140
    ax.add_patch(FancyBboxPatch((mx0, 0.478), mw, 0.278,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="#fdf4e6", ec=AMBER, lw=1.4))
    ax.text(mx0 + mw / 2, 0.724, "2 x Transformer", ha="center", fontsize=8.6,
            fontweight="bold", color=SLATE)
    for i in range(7):
        for j in range(4):
            ax.plot([mx0 + 0.022 + i * 0.0158], [0.678 - j * 0.030], "s", ms=3.2, color=AMBER,
                    alpha=0.22 + 0.10 * ((i + j) % 3), mec="none")
    ax.text(mx0 + mw / 2, 0.536, "8 heads, pre-norm, d_ff 512", ha="center", fontsize=6.3,
            color="#6b7684")
    ax.text(mx0 + mw / 2, 0.502, "all-to-all: whole-transcript context", ha="center",
            fontsize=5.9, color=AMBER, style="italic")
    shape_tag(ax, mx0 + mw / 2, SHP, "(B, 256, L)", c="#8f9aa3")
    arrow(ax, mx0 + mw + 0.004, ROW, mx0 + mw + 0.030, 0.705, c="#9aa6b2", lw=1.4, ms=10)
    arrow(ax, mx0 + mw + 0.004, ROW, mx0 + mw + 0.030, 0.525, c="#9aa6b2", lw=1.4, ms=10)

    # ---- heads ----
    hx = 0.855
    ax.add_patch(FancyBboxPatch((hx, 0.648), 0.133, 0.128,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="#e7f1f0", ec=TEAL, lw=1.3))
    ax.text(hx + 0.0665, 0.740, "profile head", ha="center", fontsize=8,
            fontweight="bold", color=SLATE)
    ax.text(hx + 0.0665, 0.706, "Conv1d(256 -> 1, k=1)", ha="center", fontsize=6.1,
            color="#5b6b70")
    ax.text(hx + 0.0665, 0.676, "softmax over L", ha="center", fontsize=6.1, color="#5b6b70")
    shape_tag(ax, hx + 0.0665, 0.628, "(B, L)   the SHAPE", c=TEAL, fs=6.4)
    ax.add_patch(FancyBboxPatch((hx, 0.462), 0.133, 0.128,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="#fdf4e6", ec=AMBER, lw=1.3))
    ax.text(hx + 0.0665, 0.554, "count head", ha="center", fontsize=8,
            fontweight="bold", color=SLATE)
    ax.text(hx + 0.0665, 0.520, "masked mean + log cov", ha="center", fontsize=6.1,
            color="#5b6b70")
    ax.text(hx + 0.0665, 0.490, "MLP 257 -> 256 -> 1", ha="center", fontsize=6.1, color="#5b6b70")
    shape_tag(ax, hx + 0.0665, 0.442, "(B,)   the DEPTH", c=AMBER, fs=6.4)

    # ---- bottom strip: inside one block ----
    ax.plot([bx0, bx0], [0.450, 0.318], color=TEAL, lw=0.8, ls=(0, (2, 2)))
    ax.plot([bx0 + bw, 0.612], [0.450, 0.318], color=TEAL, lw=0.8, ls=(0, (2, 2)))
    ax.add_patch(FancyBboxPatch((0.012, 0.070), 0.600, 0.248,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec=TEAL, lw=1.0))
    ax.text(0.030, 0.283, "inside one residual dilated block", fontsize=8,
            fontweight="bold", color=TEAL)
    ax.text(0.030, 0.252, "identical for all 10 blocks; only the dilation d changes",
            fontsize=6.2, color="#6b7684")
    steps = ["Conv1d\nk=3, dil d", "GroupNorm\n8 groups", "GELU", "Dropout\n0.1",
             "Conv1d\nk=3, dil d", "GroupNorm\n8 groups"]
    x0, w, gap = 0.032, 0.0815, 0.0125
    for i, s in enumerate(steps):
        x = x0 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 0.140), w, 0.082,
                                    boxstyle="round,pad=0,rounding_size=0.007",
                                    fc="#f2f8f7", ec="#bcd6d3", lw=0.8))
        ax.text(x + w / 2, 0.181, s, ha="center", va="center", fontsize=5.7, color=SLATE)
        if i < len(steps) - 1:
            arrow(ax, x + w, 0.181, x + w + gap, 0.181, c="#bcd6d3", lw=0.8, ms=5)
    ax.plot([0.032, 0.032], [0.140, 0.118], color=TEAL, lw=0.9)
    ax.plot([0.032, 0.5645], [0.118, 0.118], color=TEAL, lw=0.9)
    arrow(ax, 0.5645, 0.118, 0.5645, 0.139, c=TEAL, lw=0.9, ms=6)
    ax.text(0.298, 0.092, "residual add, then GELU   (padding keeps length L exactly)",
            ha="center", fontsize=6.1, color="#5c6f6d")

    # ---- objective ----
    ax.text(0.640, 0.283, "objective", fontsize=8, fontweight="bold", color=SLATE)
    ax.text(0.640, 0.243, "multinomial NLL over the profile", fontsize=6.8, color="#3f3f3f")
    ax.text(0.640, 0.208, "+   0.1 x MSE on the log total count", fontsize=6.8, color=AMBER)
    ax.plot([0.640, 0.988], [0.182, 0.182], color="#dfe5e4", lw=1.0)
    for i, t in enumerate([
            "The profile head is scale-free: it predicts only the SHAPE",
            "of translation. Sequencing depth is isolated in the count",
            "head, so the transferable signal is not entangled with",
            "dataset-specific magnitude."]):
        ax.text(0.640, 0.150 - i * 0.030, t, fontsize=6.3, color="#6b7684")
    save(fig, "arch_style2h")


# ================================================================= RINALMO STYLE
def vbox(ax, x, y, w, h, label, fs=6.4, fc="white", ec="#1a1a1a", lw=0.9, italic=False):
    """Narrow rectangle with 90-degree rotated text -- the RiNALMo Fig-1 idiom."""
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=3))
    ax.text(x + w / 2, y + h / 2, label, rotation=90, ha="center", va="center",
            fontsize=fs, color="#111111", zorder=4,
            style="italic" if italic else "normal")
    return x + w


def oplus(ax, x, y, r=0.0085, ec="#1a1a1a", lw=0.9):
    """Circled-plus residual junction (drawn round in display space)."""
    from matplotlib.patches import Ellipse
    asp = ax.figure.get_figwidth() / ax.figure.get_figheight()
    ax.add_patch(Ellipse((x, y), 2 * r, 2 * r * asp, fc="white", ec=ec, lw=lw, zorder=4))
    ax.text(x, y, "+", ha="center", va="center", fontsize=7.2, color="#111111", zorder=5)


def hline(ax, x1, x2, y, c="#1a1a1a", lw=0.9):
    ax.plot([x1, x2], [y, y], color=c, lw=lw, zorder=2, solid_capstyle="butt")


def track_grid(ax, x, y, cw, ch, rows, labels, lw=0.45, fill="#111111",
               ec="#6b7280", fs=5.0, lab_dx=0.005):
    """Orthrus-style literal encoding matrix: filled cell = 1, empty = 0, row labels at the left."""
    n = len(rows[0])
    for r, (row, lab) in enumerate(zip(rows, labels)):
        yy = y + (len(rows) - 1 - r) * ch
        ax.text(x - lab_dx, yy + ch / 2, lab, ha="right", va="center", fontsize=fs,
                family="monospace", color="#111111", zorder=5)
        for c, v in enumerate(row):
            v = float(v)
            # binary channels fill solid; graded channels (start_ext) fill proportionally grey
            fcc = "white" if v <= 0 else (fill if v >= 1 else str(round(1.0 - v, 3)))
            ax.add_patch(Rectangle((x + c * cw, yy), cw, ch,
                                   fc=fcc, ec=ec, lw=lw, zorder=4))
    return x + n * cw


def onehot_rows(seq, bases="ACGU"):
    return [[1 if ch == b else 0 for ch in seq] for b in bases]


def style_rinalmo(mixer="transformer"):
    """RiNALMo Figure-1 idiom for the network (rotated-label boxes, tinted container, dashed
    'N x' block regions, residual arcs into circled-plus) + ORTHRUS Figure-1 idiom for the inputs
    (literal black/white encoding matrices, one per track, separated by '+', split into two chunks
    with an ellipsis to denote an arbitrary-length transcript).

    mixer='transformer' -> the deployed attention model (2 transformer blocks, 5,071,106 params)
    mixer='mamba'       -> the mamba4 variant (4 bidirectional Mamba blocks, 7,521,026 params)
    Only the second dashed region differs, so the two panels are directly comparable.
    """
    LAV, GRN, AMB = "#e7e3f6", "#dcece4", "#fbf0dc"
    PAR = PARAMS if mixer == "transformer" else "7,521,026"
    fig = plt.figure(figsize=(14.0, 5.2), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ROW, BY, BH = 0.545, 0.355, 0.380
    TOP = BY + BH
    ARC = TOP + 0.050

    # ================= INPUT: Orthrus-style explicit encoding =================
    import sys
    sys.path.insert(0, str(HERE.parent.parent / "scripts"))
    from build_orf_track import START_W as _SW
    from build_orf_track import orf_track as _orf_track

    IX = 0.030                                   # left edge of the matrices
    CW, GAP = 0.0058, 0.013                      # cell width, gap between the two chunks
    # 24-nt display window chosen so the REAL track shows the thing that matters: two candidate
    # ORFs open in DIFFERENT frames at once (AUG at 0 -> f0, AUG at 4 -> f1), a near-cognate CUG
    # start, and stop codons in all three frames. Track values below are computed by the
    # production builder (scripts/build_orf_track.py), not drawn by hand.
    DSEQ = "AUGCAUGCUGAAGCCUAACUAGCC"
    NC = len(DSEQ) // 2
    S1, S2 = DSEQ[:NC], DSEQ[NC:]
    x2 = IX + NC * CW + GAP                      # start of the second chunk
    RIGHT = x2 + NC * CW

    ax.text((IX + RIGHT) / 2 - 0.006, 0.955, "Encode transcript", ha="center", fontsize=8.6,
            fontweight="bold", color="#111111")

    # transcript cartoon
    ty, th = 0.885, 0.030
    for xa, wa, fcc, lab in [(IX, 0.030, "#dfe6ee", "5' UTR"), (IX + 0.030, 0.086, "#b9cde3", "CDS"),
                             (IX + 0.116, RIGHT - IX - 0.116, "#dfe6ee", "3' UTR")]:
        ax.add_patch(Rectangle((xa, ty), wa, th, fc=fcc, ec="#5b7a99", lw=0.6, zorder=4))
        ax.text(xa + wa / 2, ty + th / 2, lab, ha="center", va="center", fontsize=4.9,
                color="#22384d", zorder=5)

    def chunked(rows, labels, y, ch, title, note=None):
        ax.text(IX - 0.019, y + len(rows) * ch + 0.016, title, fontsize=6.3,
                fontweight="bold", color="#111111", ha="left")
        if note:
            ax.text(RIGHT, y + len(rows) * ch + 0.016, note, fontsize=5.3, color="#6b7280",
                    ha="right")
        track_grid(ax, IX, y, CW, ch, [r[:NC] for r in rows], labels)
        track_grid(ax, x2, y, CW, ch, [r[NC:] for r in rows], [""] * len(rows))
        ax.text(IX + NC * CW + GAP / 2, y + len(rows) * ch / 2, "...", ha="center",
                va="center", fontsize=7, color="#111111", zorder=5)

    # 1. one-hot sequence
    oh = onehot_rows(S1 + S2)
    chunked(oh, list("ACGU"), 0.735, 0.0255, "One-hot sequence", "4 ch")
    ax.text((IX + RIGHT) / 2, 0.705, "+", ha="center", va="center", fontsize=10,
            color="#111111")

    # 2. ORF-candidate track, computed by the production builder (mode ext, --kozak none)
    trk = _orf_track(DSEQ.replace("U", "T"), mode="ext", kozak="none")
    chunked([list(trk[:, k]) for k in range(5)],
            ["f0", "f1", "f2", "start", "stop"], 0.545, 0.0255,
            "ORF-candidate track", "5 ch")
    ax.text(IX - 0.019, 0.527, f"occupancy = AUG..in-frame stop, per frame;  start graded over "
                               f"{len(_SW)} codons", fontsize=4.9, color="#6b7280")
    ax.text((IX + RIGHT) / 2, 0.515, "+", ha="center", va="center", fontsize=10, color="#111111")

    # 3. RNA-seq coverage (continuous, so drawn as a depth profile rather than a binary matrix)
    ax.text(IX - 0.019, 0.478, "RNA-seq coverage", fontsize=6.3, fontweight="bold",
            color="#111111", ha="left")
    ax.text(RIGHT, 0.478, "1 ch", fontsize=5.3, color="#6b7280", ha="right")
    cov = [0.30, 0.48, 0.66, 0.80, 0.90, 0.86, 0.74, 0.68, 0.77, 0.88, 0.83, 0.72,
           0.64, 0.71, 0.85, 0.79, 0.66, 0.57, 0.63, 0.75, 0.68, 0.52, 0.39, 0.27]
    cy, chh = 0.395, 0.058
    for i, v in enumerate(cov):
        xx = IX + i * CW if i < NC else x2 + (i - NC) * CW
        ax.add_patch(Rectangle((xx, cy), CW, chh * v, fc="#6E8B74", ec="none", alpha=0.85,
                               zorder=4))
    hline(ax, IX, IX + NC * CW, cy, c="#6b7280", lw=0.5)
    hline(ax, x2, RIGHT, cy, c="#6b7280", lw=0.5)
    ax.text(IX + NC * CW + GAP / 2, cy + chh / 2, "...", ha="center", va="center", fontsize=7,
            color="#111111", zorder=5)

    ax.add_patch(Rectangle((IX - 0.022, 0.352), RIGHT - IX + 0.030, 0.575, fc="none",
                           ec="#9aa3ab", lw=0.8, zorder=1))
    ax.text((IX + RIGHT) / 2 - 0.006, 0.325, "(B, L, 10)   per nucleotide", ha="center",
            fontsize=6.2, family="monospace", color="#6b7280")

    arrow(ax, RIGHT + 0.012, ROW, RIGHT + 0.038, ROW, c="#1a1a1a", lw=0.9, ms=8)

    # ================= ENCODER: RiNALMo idiom =================
    CX0 = RIGHT + 0.042
    CX1 = 0.884
    ax.add_patch(Rectangle((CX0, 0.275), CX1 - CX0, 0.545, fc=LAV, ec="#8b83b8", lw=1.0, zorder=1))
    ax.text((CX0 + CX1) / 2, 0.955, "RiboSignal encoder", ha="center", fontsize=9.6,
            fontweight="bold", color="#111111")
    ax.text((CX0 + CX1) / 2, 0.908, f"{PAR} parameters   |   PyTorch", ha="center",
            fontsize=6.4, color="#5b6270")

    BW, BG = 0.027, 0.012
    x = vbox(ax, CX0 + 0.014, BY, BW, BH, "1x1 Conv Projection", fs=6.0)
    ax.text(CX0 + 0.014 + BW / 2, BY - 0.042, "10 -> 256", ha="center", fontsize=5.6,
            family="monospace", color="#6b7280")
    arrow(ax, x, ROW, x + 0.020, ROW, c="#1a1a1a", lw=0.9, ms=8); x += 0.020

    # --- dashed region A: residual dilated blocks ---
    a0 = x + 0.004
    AW = 0.296
    ax.add_patch(Rectangle((a0, 0.300), AW, 0.494, fc=GRN, ec="none", zorder=2))
    ax.add_patch(Rectangle((a0, 0.300), AW, 0.494, fc="none", ec="#4b7a63", lw=0.9,
                           ls=(0, (4, 3)), zorder=6))
    bx = a0 + 0.015
    hline(ax, a0, bx, ROW)
    res0 = bx - 0.008
    lab = ["Conv1d\nk=3, dilation d", "GroupNorm", "GELU", "Dropout",
           "Conv1d\nk=3, dilation d", "GroupNorm"]
    for i, s in enumerate(lab):
        nx = vbox(ax, bx, BY, BW, BH, s, fs=5.6)
        if i < len(lab) - 1:
            hline(ax, nx, nx + BG, ROW); bx = nx + BG
        else:
            bx = nx
    hline(ax, bx, bx + 0.016, ROW)
    oplus(ax, bx + 0.026, ROW)
    hline(ax, bx + 0.034, a0 + AW, ROW)
    ax.plot([res0, res0], [ROW, ARC], color="#1a1a1a", lw=0.9, zorder=3)
    hline(ax, res0, bx + 0.026, ARC)
    ax.plot([bx + 0.026, bx + 0.026], [ARC, ROW + 0.016], color="#1a1a1a", lw=0.9, zorder=3)
    arrow(ax, bx + 0.026, ROW + 0.026, bx + 0.026, ROW + 0.015, c="#1a1a1a", lw=0.9, ms=6)
    ax.text(a0 + AW / 2, 0.322, "10 x  Residual Dilated Blocks    d = 1, 2, 4 ... 512",
            ha="center", fontsize=6.8, style="italic", color="#2f5a49", zorder=7)
    ax.text(a0 + AW / 2, 0.833, "receptive field ~4 kb, linear in L", ha="center", fontsize=6.1,
            color="#5b6270")
    x = a0 + AW
    arrow(ax, x, ROW, x + 0.020, ROW, c="#1a1a1a", lw=0.9, ms=8); x += 0.020

    # --- dashed region B: the MIXER (transformer or bidirectional Mamba) ---
    b0 = x + 0.004
    BWD = 0.232 if mixer == "transformer" else 0.250
    ax.add_patch(Rectangle((b0, 0.300), BWD, 0.494, fc=AMB, ec="none", zorder=2))
    ax.add_patch(Rectangle((b0, 0.300), BWD, 0.494, fc="none", ec="#b07d2a", lw=0.9,
                           ls=(0, (4, 3)), zorder=6))
    bx = b0 + 0.014
    hline(ax, b0, bx, ROW)

    if mixer == "transformer":
        for grp in (["LayerNorm", "Multi-Head Attention"], ["LayerNorm", "Feed Forward Network"]):
            g0 = bx - 0.007
            for j, s in enumerate(grp):
                nx = vbox(ax, bx, BY, BW, BH, s, fs=5.6)
                if j < len(grp) - 1:
                    hline(ax, nx, nx + BG, ROW); bx = nx + BG
                else:
                    bx = nx
            hline(ax, bx, bx + 0.012, ROW)
            oplus(ax, bx + 0.021, ROW)
            ax.plot([g0, g0], [ROW, ARC], color="#1a1a1a", lw=0.9, zorder=3)
            hline(ax, g0, bx + 0.021, ARC)
            ax.plot([bx + 0.021, bx + 0.021], [ARC, ROW + 0.016], color="#1a1a1a", lw=0.9, zorder=3)
            arrow(ax, bx + 0.021, ROW + 0.026, bx + 0.021, ROW + 0.015, c="#1a1a1a", lw=0.9, ms=6)
            bx += 0.029
            hline(ax, bx, bx + 0.010, ROW); bx += 0.010
        blocklab = "2 x  Transformer Blocks    8 heads, pre-norm"
        toplab = "all-to-all: whole-transcript frame context"
    else:
        # BiMambaBlock: residual + Dropout( Mamba_fwd(LN(x)) + Mamba_bwd(flip(LN(x))).flip() )
        g0 = bx - 0.007
        nx = vbox(ax, bx, BY, BW, BH, "LayerNorm", fs=5.6)
        hline(ax, nx, nx + 0.014, ROW)
        sx = nx + 0.014                                   # split point
        ytop, ybot = ROW + 0.102, ROW - 0.102             # the two scan lanes
        hh = 0.166
        ax.plot([sx, sx], [ybot, ytop], color="#1a1a1a", lw=0.9, zorder=3)
        for yy, lab in ((ytop, "Mamba  forward scan"), (ybot, "Mamba  reverse scan")):
            arrow(ax, sx, yy, sx + 0.010, yy, c="#1a1a1a", lw=0.9, ms=6)
            ax.add_patch(Rectangle((sx + 0.010, yy - hh / 2), BW, hh, fc="white", ec="#1a1a1a",
                                   lw=0.9, zorder=3))
            ax.text(sx + 0.010 + BW / 2, yy, lab, rotation=90, ha="center", va="center",
                    fontsize=5.2, zorder=4)
        mx = sx + 0.010 + BW + 0.012                      # merge point
        for yy in (ytop, ybot):
            hline(ax, sx + 0.010 + BW, mx, yy)
        ax.plot([mx, mx], [ybot, ytop], color="#1a1a1a", lw=0.9, zorder=3)
        oplus(ax, mx + 0.010, ROW)                        # sum of the two scans
        hline(ax, mx + 0.018, mx + 0.030, ROW)
        nx2 = vbox(ax, mx + 0.030, BY, BW, BH, "Dropout", fs=5.6)
        hline(ax, nx2, nx2 + 0.012, ROW)
        oplus(ax, nx2 + 0.021, ROW)                       # residual add
        ax.plot([g0, g0], [ROW, ARC], color="#1a1a1a", lw=0.9, zorder=3)
        hline(ax, g0, nx2 + 0.021, ARC)
        ax.plot([nx2 + 0.021, nx2 + 0.021], [ARC, ROW + 0.016], color="#1a1a1a", lw=0.9, zorder=3)
        arrow(ax, nx2 + 0.021, ROW + 0.026, nx2 + 0.021, ROW + 0.015, c="#1a1a1a", lw=0.9, ms=6)
        bx = nx2 + 0.030
        blocklab = "4 x  Bi-Mamba Blocks    d_state 16, d_conv 4, expand 2"
        toplab = "bidirectional SSM: full-transcript context at O(L)"
    hline(ax, bx, b0 + BWD, ROW)
    ax.text(b0 + BWD / 2, 0.322, blocklab, ha="center", fontsize=6.8, style="italic",
            color="#8a6118", zorder=7)
    ax.text(b0 + BWD / 2, 0.833, toplab, ha="center", fontsize=6.1, color="#5b6270")

    # ================= HEADS =================
    hx = CX1 + 0.020
    ax.plot([CX1, hx - 0.009], [ROW, ROW], color="#1a1a1a", lw=0.9, zorder=2)
    ax.plot([hx - 0.009, hx - 0.009], [0.420, 0.670], color="#1a1a1a", lw=0.9, zorder=2)
    for yy, name, shp, note, col in [
            (0.670, "Profile Head", "(B, L)", "softmax over positions    the SHAPE", "#2f5a49"),
            (0.420, "Count Head", "(B,)", "scalar log total    the DEPTH", "#8a6118")]:
        arrow(ax, hx - 0.009, yy, hx, yy, c="#1a1a1a", lw=0.9, ms=7)
        ax.add_patch(Rectangle((hx, yy - 0.072), 0.024, 0.144, fc="white", ec="#1a1a1a",
                               lw=0.9, zorder=3))
        ax.text(hx + 0.012, yy, name, rotation=90, ha="center", va="center", fontsize=6.0,
                zorder=4)
        arrow(ax, hx + 0.024, yy, hx + 0.040, yy, c="#1a1a1a", lw=0.9, ms=7)
        ax.text(hx + 0.046, yy + 0.019, shp, fontsize=6.3, family="monospace", va="center",
                color=col, fontweight="bold")
        ax.text(hx + 0.046, yy - 0.024, note, fontsize=5.7, va="center", color="#5b6270")

    ax.text(0.012, 0.208, "loss  =  multinomial NLL over the profile   +   0.1 x MSE on the log "
                          "total count", fontsize=7.0, color="#111111")
    ax.text(0.012, 0.157, "The profile head is scale-free (shape only); sequencing depth is isolated "
                          "in the count head, so the transferable signal is not entangled with "
                          "dataset-specific magnitude.", fontsize=6.2, color="#5b6270")
    ax.text(0.012, 0.108, "The ORF-candidate track is derived from the sequence alone (no annotated CDS). "
                          "It encodes EVERY candidate ORF, not one: f0/f1/f2 mark AUG..in-frame-stop "
                          "occupancy bucketed by the start codon's frame, so overlapping ORFs in "
                          "different frames are on simultaneously;", fontsize=6.2, color="#5b6270",
            style="italic")
    ax.text(0.012, 0.062, "'stop' marks every stop codon in any frame, and 'start' is graded over AUG "
                          "(1.0) plus 9 near-cognates (CUG 0.5, GUG/ACG 0.35 ... AAG/AGG 0.15). "
                          "Occupancy opens at AUG only.", fontsize=6.2, color="#5b6270", style="italic")
    save(fig, "arch_rinalmo" if mixer == "transformer" else "arch_rinalmo_mamba")


if __name__ == "__main__":
    print("architecture-only panels ...")
    style1(); style2(); style3(); style2h()
    style_rinalmo("transformer"); style_rinalmo("mamba")
    print(f"-> {HERE}")
