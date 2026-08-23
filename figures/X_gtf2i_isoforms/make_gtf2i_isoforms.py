#!/usr/bin/env python3
"""GTF2I in held-out human hepatocytes: observed Ribo-seq across all 35 annotated isoforms.

WHAT THIS SHOWS. A worked example of the posture-A multimapper convention (methods.md 2.3). A
footprint that is genomically unique still maps to every compatible isoform of its gene, and
posture A counts it at weight 1 on all of them. GTF2I has 35 isoforms in the union universe, so
one gene's ~13k footprints are reported 35 times over.

The consequence is visible in panel (b) and is the point of the figure: observed Ribo-seq carries
almost NO information about WHICH isoform is expressed. RNA-seq does.

  make_gtf2i_isoforms.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

NEW = Path(__file__).resolve().parents[2]   # figures/<name>/ -> project root
PACK = NEW / "data/packed_union_Hepatocytes"
UNIV = NEW / "data/union_universe.tsv"
OUT = Path(__file__).resolve().parent
GENE = "GTF2I"


def load():
    order = (PACK / "tx_order.txt").read_text().split()
    off = np.load(PACK / "offsets.npy")
    lens = np.load(PACK / "lengths.npy")
    tc = np.load(PACK / "target_counts.npy")
    idx = {t: i for i, t in enumerate(order)}
    tpm = {}
    with open(UNIV) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] == GENE:
                tpm[f[0]] = float(f[6])
    rows = []
    for t, v in tpm.items():
        i = idx[t]
        prof = tc[off[i]:off[i + 1]].astype(float)
        rows.append({"tx": t, "len": int(lens[i]), "tpm": v,
                     "total": int(prof.sum()), "profile": prof})
    rows.sort(key=lambda r: -r["tpm"])
    return rows


def main():
    rows = load()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7.2),
                                   gridspec_kw={"height_ratios": [1.35, 1]})

    # (a) per-nt observed profile, the 6 highest-RNA isoforms overlaid
    for r in rows[:6]:
        ax1.plot(np.arange(r["len"]), r["profile"], lw=0.8, alpha=0.75,
                 label=f'{r["tx"]}  TPM {r["tpm"]:.1f}')
    ax1.set_xlim(0, 5300)
    ax1.set_xlabel("position along transcript (nt)")
    ax1.set_ylabel("observed P-sites")
    ax1.set_title(f"(a) {GENE}, held-out hepatocytes: per-nt observed Ribo-seq, "
                  f"6 highest-RNA isoforms overlaid", fontsize=10, loc="left")
    ax1.legend(fontsize=6.5, ncol=2, frameon=False)

    # (b) the point: total observed P-sites vs RNA abundance, all 35
    tpms = np.array([r["tpm"] for r in rows])
    tots = np.array([r["total"] for r in rows], dtype=float)
    ax2.scatter(tpms, tots, s=26, c="#2C6FBB", zorder=3)
    for r in rows:
        if r["tpm"] > 10 or r["total"] == tots.max():
            ax2.annotate(r["tx"], (r["tpm"], r["total"]), fontsize=6.2,
                         xytext=(4, 4), textcoords="offset points")
    ax2.set_xscale("log")
    ax2.set_xlabel("RNA-seq abundance, TPM (log scale)")
    ax2.set_ylabel("total observed P-sites")
    ax2.set_ylim(0, tots.max() * 1.12)
    rng = 100 * (tots.max() - tots.min()) / tots.max()
    ax2.set_title(f"(b) all {len(rows)} isoforms. RNA spans {tpms.max()/tpms.min():.0f}x; "
                  f"observed Ribo-seq spans only {rng:.0f}% -- posture A gives them all the "
                  f"same footprints", fontsize=10, loc="left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"X_gtf2i_isoforms.{ext}", dpi=190, bbox_inches="tight")

    vals = {
        "model": None,
        "source": "data/packed_union_Hepatocytes/{target_counts,offsets,lengths,tx_order}"
                  " + data/union_universe.tsv (RNA TPM)",
        "provenance_note": "MODEL-FREE: observed Ribo-seq only, no prediction. Hepatocytes is the "
                           "LOTO held-out tissue; these are its measured P-sites under the "
                           "posture-A multimapper convention (methods.md 2.3).",
        "gene": GENE, "n_isoforms": len(rows),
        "summed_psites_over_isoforms": int(tots.sum()),
        "per_isoform_psites_min": int(tots.min()), "per_isoform_psites_max": int(tots.max()),
        "psite_spread_pct_of_max": round(rng, 2),
        "rna_tpm_min": float(tpms.min()), "rna_tpm_max": float(tpms.max()),
        "rna_tpm_fold_range": round(float(tpms.max() / tpms.min()), 1),
        "inflation_vs_single_isoform": round(float(tots.sum() / tots.max()), 1),
        "highest_ribo_isoform": rows[int(np.argmax(tots))]["tx"],
        "highest_ribo_isoform_tpm": float(rows[int(np.argmax(tots))]["tpm"]),
        "highest_rna_isoform": rows[0]["tx"], "highest_rna_isoform_tpm": float(rows[0]["tpm"]),
        "isoforms": [{"tx": r["tx"], "len": r["len"], "tpm": r["tpm"], "total_psites": r["total"]}
                     for r in rows],
    }
    (OUT / "X_gtf2i_values.json").write_text(json.dumps(vals, indent=2))
    print(f"wrote {OUT}/X_gtf2i_isoforms.png  (+pdf, values json)")
    print(f"  {len(rows)} isoforms, summed {int(tots.sum()):,} P-sites, "
          f"{tots.sum()/tots.max():.1f}x inflation over a single isoform")


if __name__ == "__main__":
    raise SystemExit(main())
