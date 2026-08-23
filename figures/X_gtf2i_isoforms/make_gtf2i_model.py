#!/usr/bin/env python3
"""GTF2I: observed P-sites coloured by reading frame, with each isoform's exon structure drawn
underneath in the SAME transcript coordinates.

WHY THE MODEL LOOKS CONTINUOUS. The x-axis is transcript position, not genome position, so exons
tile end-to-end with no introns between them. What the track carries is therefore WHERE THE
JUNCTIONS FALL (alternating shading + ticks), which is exactly what differs between isoforms and
what makes the same genomic footprint land at different transcript offsets in different isoforms.
Box height encodes CDS vs UTR, the usual gene-model convention.

Junction positions are cumulative exon lengths from the GENCODE v49 primary GTF, ordered by
`exon_number`. Verified: for all 35 GTF2I isoforms the exon lengths sum EXACTLY to the transcript
length recorded in the pack, so the transcript-coordinate mapping is correct.

  make_gtf2i_model.py --top 10
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.gridspec import GridSpec             # noqa: E402
import numpy as np                                    # noqa: E402

NEW = Path(__file__).resolve().parents[2]
PACK = NEW / "data/packed_union_Hepatocytes"
UNIV = NEW / "data/union_universe.tsv"
CDSF = NEW / "data/human_ribocode_annot_primary/transcripts_cds.txt"
GTF = Path("/private/groups/carpenterlab/emalekos/genomes/gencode.v49.primary_assembly.annotation.gtf")
EXCACHE = Path("/tmp/gtf2i_exons.gtf")
OUT = Path(__file__).resolve().parent
GENE = "GTF2I"
FCOL = {0: "#2C6FBB", 1: "#C8A45B", 2: "#B4574A"}


def exon_lengths():
    """tx -> [exon lengths in transcript order]. Uses the pre-filtered cache if present."""
    src = EXCACHE if EXCACHE.exists() else GTF
    ex: dict[str, list] = {}
    with open(src) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f'gene_name "{GENE}"' not in f[8]:
                continue
            tx = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            en = int(re.search(r"exon_number (\d+)", f[8]).group(1))
            ex.setdefault(tx, []).append((en, int(f[3]), int(f[4])))
    return {t: [e - s + 1 for _, s, e in sorted(v)] for t, v in ex.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--tag", default="")
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
    with open(CDSF) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in tpm:
                cds[f[0]] = (int(f[1]), int(f[2]))
    exl = exon_lengths()

    rows = []
    for t, v in sorted(tpm.items(), key=lambda kv: -kv[1]):
        i = idx[t]
        prof = tc[off[i]:off[i + 1]].astype(float)
        cs, ce = cds.get(t, (1, len(prof)))
        pos = np.nonzero(prof)[0]
        frame = (pos - (cs - 1)) % 3
        inc = (pos >= cs - 1) & (pos <= ce - 1)
        tin = float(prof[pos[inc]].sum())
        el = exl.get(t, [])
        assert not el or sum(el) == int(lens[i]), f"exon sum != tx length for {t}"
        rows.append({"tx": t, "tpm": v, "len": int(lens[i]), "cds": (cs, ce), "prof": prof,
                     "pos": pos, "frame": frame, "exons": el,
                     "junctions": list(np.cumsum(el)[:-1]) if el else [],
                     "total": float(prof.sum()),
                     "f0": float(prof[pos[inc & (frame == 0)]].sum()) / tin if tin else np.nan})
    if a.top:
        rows = rows[:a.top]
    n = len(rows)

    XMAX = 5300
    fig = plt.figure(figsize=(14.5, 1.62 * n + 0.7))
    gs = GridSpec(2 * n, 1, figure=fig, height_ratios=[3.4, 1.0] * n, hspace=0.0)
    ymax = max(r["prof"].max() for r in rows) * 1.08

    for k, r in enumerate(rows):
        cs, ce = r["cds"]
        axp = fig.add_subplot(gs[2 * k, 0])
        for fr in (1, 2, 0):
            m = r["frame"] == fr
            if m.any():
                axp.vlines(r["pos"][m], 0, r["prof"][r["pos"][m]], color=FCOL[fr], lw=0.75)
        axp.set_xlim(0, XMAX); axp.set_ylim(0, ymax)
        axp.set_xticklabels([]); axp.tick_params(labelsize=7)
        axp.set_ylabel("P-sites", fontsize=7)
        axp.set_title(f'{r["tx"]}   RNA TPM {r["tpm"]:.2f}   |   {int(r["total"]):,} P-sites   |   '
                      f'{len(r["exons"])} exons   |   frame-0 {100*r["f0"]:.1f}%',
                      fontsize=8.5, loc="left", pad=3)

        # ---- gene model, same transcript coordinates -------------------------------------
        axm = fig.add_subplot(gs[2 * k + 1, 0])
        start = 0
        for j, L in enumerate(r["exons"]):
            end = start + L
            # UTR thin / CDS thick; an exon may span the CDS boundary, so draw both parts
            axm.add_patch(plt.Rectangle((start, 0.40), L, 0.20, lw=0,
                                        color="#9AA5B1" if j % 2 else "#6B7684"))
            lo, hi = max(start, cs - 1), min(end, ce - 1)
            if hi > lo:
                axm.add_patch(plt.Rectangle((lo, 0.24), hi - lo, 0.52, lw=0,
                                            color="#2C6FBB" if j % 2 else "#1F4F87"))
            start = end
        for j in r["junctions"]:
            axm.plot([j, j], [0.12, 0.88], color="#111111", lw=0.5, alpha=0.55, zorder=4)
        axm.set_xlim(0, XMAX); axm.set_ylim(0, 1); axm.set_yticks([])
        for s in ("top", "right", "left"):
            axm.spines[s].set_visible(False)
        axm.tick_params(labelsize=7)
        if k < n - 1:
            axm.set_xticklabels([])
        else:
            axm.set_xlabel("position along transcript (nt)  --  exons tile end to end; "
                           "ticks are splice junctions", fontsize=9)

    h = [plt.Line2D([], [], color=FCOL[f], lw=2.6, label=f"P-sites, frame {f}") for f in (0, 1, 2)]
    h += [plt.Rectangle((0, 0), 1, 1, color="#1F4F87", label="CDS exon"),
          plt.Rectangle((0, 0), 1, 1, color="#6B7684", label="UTR exon"),
          plt.Line2D([], [], color="#111111", lw=1.0, label="splice junction")]
    fig.legend(handles=h, loc="lower center", ncol=6, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.006))
    fig.suptitle(f"{GENE}, held-out hepatocytes: top {n} isoforms by RNA TPM, observed P-sites over "
                 f"each isoform's own exon structure", fontsize=11.5, y=0.998)
    fig.tight_layout(rect=[0, 0.016, 1, 0.99])
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"X_gtf2i_model{a.tag}.{ext}", dpi=175, bbox_inches="tight")

    vals = {
        "model": None,
        "source": "data/packed_union_Hepatocytes/target_counts.npy (observed P-sites); CDS from "
                  "data/human_ribocode_annot_primary/transcripts_cds.txt; exon structure from "
                  "genomes/gencode.v49.primary_assembly.annotation.gtf; TPM from union_universe.tsv",
        "provenance_note": "MODEL-FREE: observed Ribo-seq only. Exon lengths verified to sum "
                           "EXACTLY to the packed transcript length for all 35 isoforms.",
        "gene": GENE, "n_isoforms": n, "top_n_by_tpm": a.top or None,
        "isoforms": [{"tx": r["tx"], "tpm": r["tpm"], "len": r["len"], "n_exons": len(r["exons"]),
                      "cds_start": r["cds"][0], "cds_end": r["cds"][1],
                      "total_psites": int(r["total"]),
                      "frame0_frac_in_cds": round(float(r["f0"]), 4),
                      "junctions_nt": [int(x) for x in r["junctions"]]} for r in rows],
    }
    (OUT / f"X_gtf2i_model{a.tag}_values.json").write_text(json.dumps(vals, indent=2))
    print(f"wrote {OUT}/X_gtf2i_model{a.tag}.png")
    print(f"  exons: {min(len(r['exons']) for r in rows)} to {max(len(r['exons']) for r in rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
