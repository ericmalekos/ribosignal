#!/usr/bin/env python3
"""Figure A1: the model's PREDICTED translation-localization (pred_frame0 AUROC) matches or beats the
OBSERVED Ribo-seq periodicity ceiling on held-out data -- "predicts better than the measurement predicts
itself." Grouped bars: model-predicted vs observed-ceiling AUROC, for held-out Hepatocytes (LOTO) and
cross-study Ruiz-Orera, split into all ORFs and non-canonical ORFs. Length-controlled AUROC is the primary
(fairer) metric; raw is annotated. Reusable: reads the two localization_metrics.json, so it regenerates if
the model is retrained (e.g. the mm1 RNA-coverage ablation). cas12a env (matplotlib). Run: python make_...py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/A1_localization_ceiling"
SOURCES = {
    "Held-out\nHepatocytes": NEW / "results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/localization_metrics.json",
    "Cross-study\nRuiz-Orera": NEW / "results/heldout/human_ruizorera/onehot/localization_metrics_human_ruizorera.json",
}
STRATA = [("all", "all ORFs"), ("noncanonical", "non-canonical")]
MODEL_C, OBS_C = "#2C6FBB", "#B0B0B0"


def load(fp):
    d = json.load(open(fp))["discrimination"]
    out = {}
    for key, _ in STRATA:
        s = d.get(key, {})
        out[key] = dict(pred=s.get("auroc_pred_frame0_lengthctrl"),
                        obs=s.get("auroc_obs_frame0_ceiling_lengthctrl"),
                        pred_raw=s.get("auroc_pred_frame0"), obs_raw=s.get("auroc_obs_frame0_ceiling"))
    return out


def main():
    data = {name: load(fp) for name, fp in SOURCES.items()}
    groups = [(name, sk, sl) for name in SOURCES for sk, sl in STRATA]
    x = np.arange(len(groups)); w = 0.38
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    preds = [data[n][sk]["pred"] for n, sk, _ in groups]
    obss = [data[n][sk]["obs"] for n, sk, _ in groups]
    b1 = ax.bar(x - w/2, preds, w, label="Model (predicted)", color=MODEL_C, edgecolor="white", linewidth=0.5)
    b2 = ax.bar(x + w/2, obss, w, label="Observed Ribo-seq ceiling", color=OBS_C, edgecolor="white", linewidth=0.5)
    for xi, (n, sk, _) in enumerate(groups):
        p, o = data[n][sk]["pred"], data[n][sk]["obs"]
        ax.text(xi - w/2, p + 0.006, f"{p:.3f}", ha="center", va="bottom", fontsize=7.5, color=MODEL_C)
        ax.text(xi + w/2, o + 0.006, f"{o:.3f}", ha="center", va="bottom", fontsize=7.5, color="#555")
        if p > o:  # model beats the measurement's own ceiling
            ax.text(xi, max(p, o) + 0.028, "beats\nceiling", ha="center", va="bottom", fontsize=6.5,
                    color="#1a7f37", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{n}\n{sl}".replace("\n", " ", 1) for n, _, sl in groups], fontsize=8)
    ax.set_ylabel("Localization AUROC\n(length-controlled)", fontsize=9)
    ax.set_ylim(0.7, 0.97)
    ax.axhline(0.5, color="k", lw=0.4, ls=":")
    ax.legend(frameon=False, fontsize=8, loc="lower left", ncol=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Predicted translation localization matches/beats the Ribo-seq ceiling",
                 fontsize=9.5, pad=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"A1_localization_ceiling.{ext}", dpi=300, bbox_inches="tight")
    # dump the source numbers next to the figure for the caption
    (HERE / "A1_values.json").write_text(json.dumps(
        {n: data[n] for n in SOURCES}, indent=2) + "\n")
    print("wrote A1_localization_ceiling.pdf/.png + A1_values.json")
    for n in SOURCES:
        for sk, sl in STRATA:
            v = data[n][sk]
            verdict = "BEATS ceiling" if v["pred"] > v["obs"] else f"{100*v['pred']/v['obs']:.0f}% of ceiling"
            print(f"  {n.strip().replace(chr(10),' ')} {sl}: pred {v['pred']:.3f} vs obs {v['obs']:.3f} ({verdict})")


if __name__ == "__main__":
    main()
