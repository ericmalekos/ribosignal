#!/usr/bin/env python3
"""GTF2I (ENST00000901263.1) in held-out hepatocytes: the POSTER's attn prediction beside mamba4,
against the same observed track.

Both profiles come from each run's own `dropin/pred_profiles.npz`, so the two arms are read the
same way rather than one from a prepared exemplar TSV and one recomputed. The transcript is the
poster's panel-B annotated-CDS exemplar (5,178 nt, called ORF 451..3324).

  make_gtf2i_attn_vs_mamba.py
"""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
from scipy.stats import pearsonr          # noqa: E402

NEW = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TX = "ENST00000901263.1"
CDS = (451, 3324)                          # 1-based inclusive, from human_ribocode_annot_primary
RUN = "results/loto/orf_v2_{m}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/pred_profiles.npz"
COL = {"attn": "#2C6FBB", "mamba4": "#B4574A"}


def grab(m):
    z = np.load(NEW / RUN.format(m=m), allow_pickle=False)
    tx = [str(t) for t in z["tx_ids"]]
    i = tx.index(TX)
    off = np.concatenate([[0], np.cumsum(z["lengths"])])
    a, b = off[i], off[i + 1]
    return (np.asarray(z["pred_flat"][a:b], dtype=float),
            np.asarray(z["obs_flat"][a:b], dtype=float),
            float(z["pred_total"][i]))


def main():
    pa, obs, ta = grab("attn")
    pm, obs2, tm = grab("mamba4")
    assert np.array_equal(obs, obs2), "observed track differs between runs"
    L = len(obs); cs, ce = CDS
    # profiles are scale-free shapes; put them on the observed scale for a like-for-like picture
    sa = pa / max(pa.sum(), 1e-12) * obs.sum()
    sm = pm / max(pm.sum(), 1e-12) * obs.sum()
    r = {"attn": pearsonr(obs, sa)[0], "mamba4": pearsonr(obs, sm)[0],
         "attn_vs_mamba": pearsonr(sa, sm)[0]}

    fig, axes = plt.subplots(3, 1, figsize=(13, 7.2), sharex=True)
    for ax, (y, lab, c) in zip(axes, [(obs, "observed", "#6B7684"),
                                      (sa, f'attn  (POSTER)   r={r["attn"]:.3f}', COL["attn"]),
                                      (sm, f'mamba4            r={r["mamba4"]:.3f}', COL["mamba4"])]):
        ax.axvspan(cs - 1, ce - 1, color="#000000", alpha=0.05, lw=0)
        ax.vlines(np.arange(L), 0, y, color=c, lw=0.5)
        ax.set_ylabel("P-sites", fontsize=8)
        ax.set_title(lab, fontsize=9.5, loc="left", pad=3)
        ax.tick_params(labelsize=7)
    axes[-1].set_xlabel("position along transcript (nt)   |   grey band = called ORF 451-3324",
                        fontsize=9)
    axes[0].set_xlim(0, L)
    fig.suptitle(f"GTF2I {TX}, held-out hepatocytes: poster's attn prediction vs mamba4 "
                 f"(attn-vs-mamba r={r['attn_vs_mamba']:.3f})", fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"X_gtf2i_attn_vs_mamba.{e}", dpi=180, bbox_inches="tight")

    def instats(y):
        inc = y[cs - 1:ce - 1]
        f0 = inc[::3].sum()
        return dict(frac_in_orf=round(float(inc.sum() / y.sum()), 4),
                    frame0_in_orf=round(float(f0 / max(inc.sum(), 1e-12)), 4),
                    max_pos=int(np.argmax(y)), top1_frac=round(float(y.max() / y.sum()), 4))
    vals = {"model": ["attn", "mamba4"],
            "source": "results/loto/orf_v2_{attn,mamba4}_onehot_union_noBrain_nokozak_mm1_"
                      "holdout_Hepatocytes/dropin/pred_profiles.npz",
            "provenance_note": "Both arms read from their own dropin npz so they are directly "
                               "comparable. Predicted SHAPES are rescaled to the observed total "
                               "for plotting; the profile head is scale-free.",
            "transcript": TX, "length": L, "cds": list(CDS),
            "obs_total": int(obs.sum()), "pred_total_attn": ta, "pred_total_mamba4": tm,
            "pearson": {k: round(float(v), 4) for k, v in r.items()},
            "observed": instats(obs), "attn": instats(sa), "mamba4": instats(sm)}
    (OUT / "X_gtf2i_attn_vs_mamba_values.json").write_text(json.dumps(vals, indent=2))
    print(f"  obs total {int(obs.sum()):,}   attn r={r['attn']:.4f}   mamba4 r={r['mamba4']:.4f}"
          f"   attn-vs-mamba r={r['attn_vs_mamba']:.4f}")
    for k in ("observed", "attn", "mamba4"):
        v = vals[k]
        print(f"  {k:<9} in-ORF {100*v['frac_in_orf']:>5.1f}%  frame0 {100*v['frame0_in_orf']:>5.1f}%"
              f"  top1 {100*v['top1_frac']:>4.1f}%  peak@{v['max_pos']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
