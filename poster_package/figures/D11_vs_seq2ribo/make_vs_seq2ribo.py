#!/usr/bin/env python3
"""Figure D11: this model vs the purpose-built sequence->Ribo-seq model seq2ribo, on WITHIN-transcript
profile shape (per-codon CDS Shape r, elongation only). This model (0.526) substantially outperforms
seq2ribo's published within-transcript Shape r (0.05-0.19), because this model is trained on the
within-transcript profile (multinomial NLL) while seq2ribo optimizes cross-transcript magnitude (its
advertised 0.92 is 'Elemwise r', pooled across transcripts ~ expression, NOT within-transcript shape). The
reproducible ceiling (replicate concordance, per-codon) is shown as the upper bound. Values from results.md
Task 22/per-codon section. cas12a env.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/D11_vs_seq2ribo"
REP = NEW / "results/replicate_concordance.json"

THIS_MODEL = 0.526                 # orf_v2_attn fold-0, per-codon CDS Shape r (results.md)
SEQ2RIBO = (0.05, 0.19)            # published within-transcript "Shape r" range (Kaynar & Kingsford 2026)


def sb(r):
    return 2 * r / (1 + r)


def main():
    rep = json.load(open(REP))
    ceiling = sum(sb(s["median"]["pc_cds_codon"]) for s in rep["per_seed"]) / len(rep["per_seed"])
    labels = ["seq2ribo\n(published)", "This model", "Reproducible\nceiling"]
    vals = [(SEQ2RIBO[0] + SEQ2RIBO[1]) / 2, THIS_MODEL, ceiling]
    errs = [[(SEQ2RIBO[1] - SEQ2RIBO[0]) / 2], [0], [0]]
    cols = ["#C8A45B", "#2C6FBB", "#B0B0B0"]
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    x = range(len(labels))
    ax.bar(x, vals, 0.6, color=cols, edgecolor="white", lw=0.5)
    # seq2ribo range whisker
    ax.errorbar(0, vals[0], yerr=errs[0], fmt="none", ecolor="#6b5628", capsize=4, lw=1.2)
    ax.text(0, SEQ2RIBO[1] + 0.02, f"{SEQ2RIBO[0]:.2f}-{SEQ2RIBO[1]:.2f}", ha="center", va="bottom", fontsize=7.5, color="#6b5628")
    ax.text(1, THIS_MODEL + 0.02, f"{THIS_MODEL:.3f}", ha="center", va="bottom", fontsize=8, color="#2C6FBB")
    ax.text(2, ceiling + 0.02, f"{ceiling:.2f}", ha="center", va="bottom", fontsize=7.5, color="#555")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("Within-transcript profile shape\n(per-codon CDS Shape r)", fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Beats the purpose-built seq2ribo on within-transcript shape",
                 fontsize=9.5)
    ann = "seq2ribo's advertised 0.92 is cross-transcript 'Elemwise r'\n(~ expression), not within-transcript shape"
    ax.text(0.98, 0.55, ann, transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color="#666", style="italic")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"D11_vs_seq2ribo.{ext}", dpi=300, bbox_inches="tight")
    (HERE / "D11_values.json").write_text(json.dumps(
        {"this_model": THIS_MODEL, "seq2ribo_range": SEQ2RIBO, "ceiling": ceiling}, indent=2) + "\n")
    print(f"  this model {THIS_MODEL} vs seq2ribo {SEQ2RIBO} vs ceiling {ceiling:.3f}")


if __name__ == "__main__":
    main()
