#!/usr/bin/env python3
"""S_metagene -- positional P-site metagene around the start and stop codons, training tissues.

Built for the poster's training-data panel, which had per-library QC SCALARS
(`S_riboseq_qc`) but no positional profile: RiboCode's `metaplots` computes these
quantities to pick P-site offsets and does not persist them.

POSITION IS THE P-SITE, NOT THE READ 5' END -- state this in any caption. The packs'
`target_counts.npy` is RiboCode `process_bam` output with the per-read-length offset
already applied, so offset 0 is the first nucleotide of the start codon (the A of AUG),
and for the stop anchor the first nucleotide of the stop codon. A 5'-end metagene would
sit ~12 nt to the left.

The data is produced by `scripts/build_metagene.py`; this only draws it. Both are portable
TSV, so the panel can be restyled off-cluster without re-reading any pack.

cas12a env. Regenerate:
    python3 scripts/build_metagene.py --out results/metagene
    cd figures/S_metagene && python3 make_metagene.py
"""
from __future__ import annotations

import csv
import json
import pathlib
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
                   "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/S_metagene"
SRC = NEW / "results/metagene"
FRAME_C = {0: "#2b6cb0", 1: "#a0aec0", 2: "#cbd5e0"}


def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        r["offset_nt"] = int(r["offset_nt"])
        r["psites"] = int(r["psites"])
        r["n_tx"] = int(r["n_tx"])
        r["psites_norm"] = float(r["psites_norm"])
    return rows


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    pooled = load(SRC / "metagene_pooled.tsv")
    bytis = load(SRC / "metagene_by_tissue.tsv")
    meta = json.loads((SRC / "metagene_meta.json").read_text())

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 6.4), constrained_layout=True)

    for col, anchor in enumerate(("start", "stop")):
        sel = sorted((r for r in pooled if r["anchor"] == anchor), key=lambda r: r["offset_nt"])
        x = np.array([r["offset_nt"] for r in sel])
        # PER-TRANSCRIPT, not the raw sum: n_tx varies with offset because a transcript contributes
        # to an offset only if that offset exists in it. Plotting the raw sum makes the window edges
        # decay purely from missing UTR.
        y = np.array([r["psites"] / max(1, r["n_tx"]) for r in sel])
        ax = axes[0][col]
        ax.bar(x, y, width=1.0, color=[FRAME_C[int(v) % 3] for v in (x if anchor == "start"
                                                                     else x)], zorder=3)
        ax.axvline(0, color="#c05621", lw=1.1, ls="--", zorder=4)
        lab = "start codon (A of AUG)" if anchor == "start" else "stop codon"
        ax.set_title(f"{'AB'[col]}  Pooled, 8 training tissues -- {anchor} anchor",
                     fontsize=9.5, loc="left")
        ax.set_xlabel(f"nt from the first base of the {lab}   (P-site)", fontsize=8.5)
        ax.set_ylabel("P-sites per transcript", fontsize=8.5)
        ax.set_yscale("log")
        ax.grid(axis="y", alpha=0.25, zorder=0)

        ax2 = axes[1][col]
        for t in sorted({r["tissue"] for r in bytis}):
            s = sorted((r for r in bytis if r["tissue"] == t and r["anchor"] == anchor),
                       key=lambda r: r["offset_nt"])
            xs = np.array([r["offset_nt"] for r in s])
            # Each tissue normalised within its own window, so 8 libraries of very different depth
            # (28.8M to 330.8M windowed P-sites) are comparable in shape.
            ys = np.array([r["psites"] / max(1, r["n_tx"]) for r in s], dtype=float)
            ys = ys / ys.sum()
            ax2.plot(xs, ys, lw=0.9, alpha=0.85, label=t)
        ax2.axvline(0, color="#c05621", lw=1.1, ls="--")
        ax2.set_title(f"{'CD'[col]}  Per tissue, shape-normalised", fontsize=9.5, loc="left")
        ax2.set_xlabel(f"nt from the {lab}   (P-site)", fontsize=8.5)
        ax2.set_ylabel("fraction of windowed P-sites", fontsize=8.5)
        ax2.grid(alpha=0.25)
        if col == 0:
            ax2.legend(fontsize=6.4, frameon=False, ncol=2)

    fr = meta.get("cds_frame_fraction_first150nt", {})
    fig.suptitle("Ribo-seq metagene, 8 final-recipe training tissues (71,116 transcripts with an "
                 "annotated start codon)\n"
                 f"P-SITE positions, offset correction already applied. CDS frame-0 = "
                 f"{fr.get('frame0', float('nan')):.1%} over the first 150 nt "
                 f"(frame1 {fr.get('frame1', 0):.1%}, frame2 {fr.get('frame2', 0):.1%})",
                 fontsize=10)
    fig.savefig(HERE / "S_metagene.pdf")
    fig.savefig(HERE / "S_metagene.png", dpi=200)
    plt.close(fig)

    # Ship the portable data beside the figure so the package is self-contained.
    for f in ("metagene_pooled.tsv", "metagene_by_tissue.tsv", "metagene_meta.json"):
        shutil.copy2(SRC / f, HERE / f)

    st = [r for r in pooled if r["anchor"] == "start"]
    utr = [r["psites"] / max(1, r["n_tx"]) for r in st if -50 <= r["offset_nt"] <= -4]
    cds = [r["psites"] / max(1, r["n_tx"]) for r in st if 0 <= r["offset_nt"] <= 149]
    peak = next(r["psites"] / max(1, r["n_tx"]) for r in st if r["offset_nt"] == 0)
    vals = {"position_semantics": meta["position_semantics"], "anchor_zero": meta["anchor_zero"],
            "n_tx": meta["packs"][0]["n_tx"], "n_tissues": len(meta["packs"]),
            "tissues": [m["tissue"] for m in meta["packs"]],
            "cds_frame_fraction_first150nt": fr,
            "start_peak_psites_per_tx": peak,
            "mean_5utr_psites_per_tx": float(np.mean(utr)),
            "mean_cds_psites_per_tx": float(np.mean(cds)),
            "start_peak_over_utr": peak / float(np.mean(utr)),
            "cds_over_utr": float(np.mean(cds)) / float(np.mean(utr)),
            "start_window": meta["start_window"], "stop_window": meta["stop_window"],
            "source": str(SRC)}
    (HERE / "S_metagene_values.json").write_text(json.dumps(vals, indent=2))
    print(f"  wrote {HERE}/S_metagene.{{pdf,png}}")
    print(f"    start-codon peak {peak:.1f} P-sites/tx = {vals['start_peak_over_utr']:.0f}x the "
          f"5'UTR baseline; CDS/UTR {vals['cds_over_utr']:.1f}x; frame0 {fr.get('frame0', 0):.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
