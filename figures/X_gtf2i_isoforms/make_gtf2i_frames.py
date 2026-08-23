#!/usr/bin/env python3
"""GTF2I, held-out hepatocytes: one panel per isoform, observed P-sites coloured by reading frame.

FRAME IS RELATIVE TO EACH ISOFORM'S OWN ANNOTATED CDS START, not to the transcript 5' end. For
position i (0-based) with CDS starting at 1-based `cds_start`, frame = (i - (cds_start-1)) % 3.
Frame 0 is in-frame with that isoform's CDS, and a genuinely translated ORF sits overwhelmingly in
frame 0. Because each isoform has a DIFFERENT cds_start, the same shared footprints get assigned
DIFFERENT frames from isoform to isoform -- which is the whole point of the figure.

Positions are P-sites: target_counts.npy is RiboCode process_bam output with the per-read-length
offset already applied, so index 0 is the first nucleotide of the transcript.

  make_gtf2i_frames.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

NEW = Path(__file__).resolve().parents[2]
PACK = NEW / "data/packed_union_Hepatocytes"
UNIV = NEW / "data/union_universe.tsv"
CDS = NEW / "data/human_ribocode_annot_primary/transcripts_cds.txt"
OUT = Path(__file__).resolve().parent
GENE = "GTF2I"
FCOL = {0: "#2C6FBB", 1: "#C8A45B", 2: "#B4574A"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=0,
                    help="keep only the N highest-RNA isoforms (0 = all)")
    ap.add_argument("--ncol", type=int, default=5)
    ap.add_argument("--tag", default="", help="output filename suffix")
    a = ap.parse_args()
    order = (PACK / "tx_order.txt").read_text().split()
    off, lens = np.load(PACK / "offsets.npy"), np.load(PACK / "lengths.npy")
    tc = np.load(PACK / "target_counts.npy")
    idx = {t: i for i, t in enumerate(order)}

    tpm = {}
    with open(UNIV) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] == GENE:
                tpm[f[0]] = float(f[6])
    cds = {}
    with open(CDS) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in tpm:
                cds[f[0]] = (int(f[1]), int(f[2]))   # 1-based inclusive

    rows = []
    for t, v in sorted(tpm.items(), key=lambda kv: -kv[1]):
        i = idx[t]
        prof = tc[off[i]:off[i + 1]].astype(float)
        cs, ce = cds.get(t, (1, len(prof)))
        pos = np.nonzero(prof)[0]
        frame = (pos - (cs - 1)) % 3
        inc = (pos >= cs - 1) & (pos <= ce - 1)          # inside this isoform's CDS
        f0 = float(prof[pos[inc & (frame == 0)]].sum())
        tot_in = float(prof[pos[inc]].sum())
        rows.append({"tx": t, "tpm": v, "len": int(lens[i]), "cds": (cs, ce), "prof": prof,
                     "pos": pos, "frame": frame,
                     "total": float(prof.sum()), "in_cds": tot_in,
                     "f0_frac": (f0 / tot_in) if tot_in else float("nan")})

    if a.top:
        rows = rows[:a.top]                     # rows are already sorted by TPM descending
    n = len(rows)
    ncol = min(a.ncol, n)
    nrow = int(np.ceil(n / ncol))
    # One transcript per row gives the full page width to a ~5 kb profile, which is what makes
    # individual codon peaks readable; the 5-wide grid is for the all-35 overview.
    pw, ph = (15.5 / ncol) * ncol, (1.55 if ncol > 2 else 1.15)
    fig, axes = plt.subplots(nrow, ncol, figsize=(pw, ph * nrow + 0.6), sharex=True, squeeze=False)
    axes = np.asarray(axes).ravel()
    ymax = max(r["prof"].max() for r in rows) * 1.08
    for ax, r in zip(axes, rows):
        cs, ce = r["cds"]
        ax.axvspan(cs - 1, ce - 1, color="#000000", alpha=0.055, lw=0)
        for fr in (1, 2, 0):                              # draw frame 0 last, on top
            m = r["frame"] == fr
            if m.any():
                ax.vlines(r["pos"][m], 0, r["prof"][r["pos"][m]], color=FCOL[fr],
                          lw=(0.55 if ncol > 2 else 0.75))
        ax.set_ylim(0, ymax)
        ax.set_xlim(0, 5300)
        ax.set_title(f'{r["tx"]}\nTPM {r["tpm"]:.1f} | {int(r["total"]):,} P-sites | '
                     f'frame-0 {100*r["f0_frac"]:.0f}%', fontsize=(6.0 if ncol > 2 else 8.0), pad=2.2)
        ax.tick_params(labelsize=(5.5 if ncol > 2 else 7.0))
    for ax in axes[n:]:
        ax.axis("off")
    h = [plt.Line2D([], [], color=FCOL[f], lw=2.4, label=f"frame {f}") for f in (0, 1, 2)]
    h.append(plt.Line2D([], [], color="#000000", alpha=0.13, lw=8, label="annotated CDS"))
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.012))
    fig.suptitle(f"{GENE}, held-out hepatocytes: observed P-sites per isoform"
                 f"{f' -- top {a.top} by RNA TPM' if a.top else ''}, coloured by frame "
                 f"relative to THAT isoform's own CDS start", fontsize=11, y=1.003)
    fig.supxlabel("position along transcript (nt)", fontsize=9, y=0.012)
    fig.tight_layout(rect=[0, 0.022, 1, 0.995])
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"X_gtf2i_frames{a.tag}.{ext}", dpi=175, bbox_inches="tight")

    f0 = np.array([r["f0_frac"] for r in rows])
    vals = {
        "model": None,
        "source": "data/packed_union_Hepatocytes/target_counts.npy (observed P-sites); CDS from "
                  "data/human_ribocode_annot_primary/transcripts_cds.txt; TPM from "
                  "data/union_universe.tsv",
        "provenance_note": "MODEL-FREE: observed Ribo-seq only. Frame is relative to each "
                           "isoform's OWN annotated CDS start, so identical shared footprints "
                           "receive different frame labels across isoforms.",
        "gene": GENE, "n_isoforms": n, "top_n_by_tpm": (a.top or None),
        "frame0_frac_min": round(float(np.nanmin(f0)), 4),
        "frame0_frac_max": round(float(np.nanmax(f0)), 4),
        "frame0_frac_median": round(float(np.nanmedian(f0)), 4),
        "isoforms": [{"tx": r["tx"], "tpm": r["tpm"], "cds_start": r["cds"][0],
                      "cds_end": r["cds"][1], "total_psites": int(r["total"]),
                      "psites_in_cds": int(r["in_cds"]),
                      "frame0_frac_in_cds": round(r["f0_frac"], 4)} for r in rows],
    }
    (OUT / f"X_gtf2i_frames{a.tag}_values.json").write_text(json.dumps(vals, indent=2))
    print(f"wrote {OUT}/X_gtf2i_frames.png")
    print(f"  frame-0 fraction in CDS: min {np.nanmin(f0):.3f}  median {np.nanmedian(f0):.3f}  "
          f"max {np.nanmax(f0):.3f}")
    for r in sorted(rows, key=lambda r: -r["f0_frac"])[:3]:
        print(f"    best  {r['tx']}  TPM {r['tpm']:>6.2f}  frame-0 {100*r['f0_frac']:.1f}%")
    for r in sorted(rows, key=lambda r: r["f0_frac"])[:3]:
        print(f"    worst {r['tx']}  TPM {r['tpm']:>6.2f}  frame-0 {100*r['f0_frac']:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
