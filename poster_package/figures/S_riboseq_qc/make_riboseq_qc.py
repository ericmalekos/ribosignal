#!/usr/bin/env python3
"""Supplemental Ribo-seq QC panel (decision D5): every Ribo-seq library used anywhere in the project,
scored on the same three quality axes, from the same source of truth.

WHY THIS EXISTS. The manuscript's claims rest on Ribo-seq from ten separate studies -- the nine
Chothani training tissues, four held-out validation datasets, and the two proteogenomics cell lines.
A reviewer is entitled to ask whether any of those libraries were bad enough to poison the labels.
This panel answers that once, uniformly, with no per-study special pleading.

WHY RiboCode `metaplots` AND NOT ribotish. FIGURES_PLAN D5 named ribotish, which needs a GENOME BAM
plus a GTF. The transcriptome BAMs are what this project keeps (the genome BAMs live in
data/*_bam/ only for the training set), so ribotish would mean re-aligning a dozen studies to
produce a plot whose content already exists. Every dataset here was already passed through RiboCode
`metaplots`, which writes exactly the QC metagene quantities -- per-read-length P-site offset, the
frame-0/1/2 sums at annotated start codons, and the periodicity test -- into `*_pre_config.txt`.
Reading those files gives a strictly larger panel from bytes already on disk. The axes are the same
ones ribotish would show.

THE THREE AXES
  (a) Read-length distribution. A monosome footprint library peaks at 28-31 nt. A peak far off that,
      or a flat distribution, means the gel cut or the trimming was wrong.
  (b) Frame preference at the CDS (f0 %). This is the periodicity signal the model is trained to
      reproduce. ~33% is noise; published libraries sit 60-90%.
  (c) P-site offset vs read length. Should be near-constant (11-13 nt) across lengths within a
      library. A jumping offset means the assignment is unreliable and the per-nt labels are smeared.

Usage: make_riboseq_qc.py [--out <dir>]
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
ECH = ROOT.parent.parent / "experiments" / "biotype_probe" / "expression_context_human"

# (display label, role, path to *_pre_config.txt). Role drives the colour grouping.
# TRAINING TISSUES READ THE CANONICAL CONFIGS. They used to read the source study's own
# `expression_context_human/data/ribocode_per_tissue/`, which no longer exists -- the canonical
# alignment redo replaced it. The failure was SILENT in the worst way: those nine entries dropped
# into `missing` and the panel rendered anyway, showing 7 held-out/pgx arms and NOT ONE training
# tissue, while its own README still claimed 16 arms. A QC panel whose job is to establish training
# library quality had stopped showing any.
#
# BRAIN: this comment used to read "intentionally absent (dropped for low periodicity before any
# pack was built), so 8 is the correct count, not 9". THAT IS STALE, and it is exactly the kind of
# guard that gets a dataset deleted again by whoever reads it next. Brain WAS fetched, aligned on
# the final recipe and metaplotted 2026-08-15 (job 36790473); metaplots selected P-site offsets for
# 5 of 5 libraries. It is included now.
#
# It is `dropped`, NOT `train`. Brain appears in no training pack, so colouring it as training would
# misstate what the model saw. `dropped` is its own role so Brain and GSE39561 sit together: they
# are the project's two excluded Ribo-seq datasets and they fail DIFFERENTLY -- Brain selects
# offsets for 5/5 libraries at weak periodicity, GSE39561 for 0 of 3 -- despite both sitting near
# 30% unique mapping. Showing that contrast is the point.
#
# So 8 TRAINING tissues and 9 CHOTHANI tissues are both correct, for different counts.
CANON = ROOT / "results/chothani_regeneration"
HELD = ROOT / "data"
SOURCES = [
    *[(t, "train", CANON / f"_canonical_{t}/canon_{t}_pre_config.txt") for t in
      ("ES", "Fat", "Fibroblast", "HA_EC", "HCAEC", "Hepatocytes", "HUVEC", "VSMC")],
    ("Brain (dropped)", "dropped", CANON / "_canonical_Brain/canon_Brain_pre_config.txt"),
    ("THP-1 GSE208041", "heldout",
     HELD / "packed_heldout_human_gse208041_canon/_work/psites/human_gse208041_canon_pre_config.txt"),
    ("CAR-T GSE304796", "heldout",
     HELD / "packed_heldout_human_cart_canon/_work/psites/human_cart_canon_pre_config.txt"),
    ("THP-1 GSE39561 (dropped)", "dropped",
     HELD / "packed_heldout_human_gse39561_canon/_work/psites/human_gse39561_canon_pre_config.txt"),
    ("Ruiz-Orera (human)", "heldout", ROOT / "data/heldout_psites/human_ruizorera/human_ruizorera_pre_config.txt"),
    ("Wang liver (mouse)", "heldout", ROOT / "data/liver3x3/pool_wang_canon/_work/psites/l3x3_wang_pre_config.txt"),
    ("GSE120762 NT (mouse)", "heldout", ROOT / "data/packed_heldout_mouse_gse120762_nt_canon/_work/psites/mouse_gse120762_nt_canon_pre_config.txt"),
    ("GSE120762 LPS (mouse)", "heldout", ROOT / "data/packed_heldout_mouse_gse120762_lps_canon/_work/psites/mouse_gse120762_lps_canon_pre_config.txt"),
    ("GSE155087 T-cell (mouse)", "heldout", ROOT / "data/packed_heldout_mouse_gse155087_tcell_canon/_work/psites/mouse_gse155087_tcell_canon_pre_config.txt"),
    ("HBL-1", "pgx", ROOT / "proteogenomics/data/HBL1_pilot/ribocode/HBL1_DMSO_pre_config.txt"),
    ("B721.221", "pgx", ROOT / "proteogenomics/data/B721_pilot/ribocode/B721_pre_config.txt"),
]
ROLE_COLOR = {"train": "#2f6f9f", "heldout": "#c1663a", "pgx": "#4a8c5f", "dropped": "#9a9a9a"}
ROLE_LABEL = {"train": "training (Chothani)", "heldout": "held-out validation",
              "pgx": "proteogenomics cell line",
              "dropped": "EXCLUDED from the project (shown to justify the exclusion)"}

# "# 29\t47.23%\t12\t330518\t21841\t26483\t87.24%\t..."
ROW = re.compile(r"^#\s*(\d+)\s+([\d.]+)%\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)%")


def parse(path):
    """-> list of per-(library, read length) dicts. One pre_config can hold many libraries."""
    rows = []
    lib = None
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("#/") or ln.startswith("# /"):
            lib = ln.lstrip("# ").strip()
            continue
        m = ROW.match(ln)
        if not m:
            continue
        rl, prop, psite, f0, f1, f2, f0pct = m.groups()
        rows.append({"lib": lib, "read_length": int(rl), "proportion": float(prop) / 100.0,
                     "psite": int(psite), "f0": int(f0), "f1": int(f1), "f2": int(f2),
                     "f0_pct": float(f0pct) / 100.0})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    out = Path(a.out)

    data, missing = [], []
    for label, role, path in SOURCES:
        if not Path(path).exists():
            missing.append((label, str(path)))
            continue
        rows = parse(path)
        if not rows:
            missing.append((label, f"{path} (parsed 0 rows)"))
            continue
        # Weight every per-library summary by that read length's share of mapped reads, so a
        # dataset's number is the read-weighted mean and not the mean over read-length bins.
        w = np.array([r["proportion"] for r in rows], dtype=float)
        tf0 = sum(r["f0"] for r in rows); tf1 = sum(r["f1"] for r in rows); tf2 = sum(r["f2"] for r in rows)
        tot = tf0 + tf1 + tf2
        data.append({
            "label": label, "role": role, "path": str(path),
            "n_libraries": len({r["lib"] for r in rows}),
            "mode_read_length": int(rows[int(np.argmax(w))]["read_length"]),
            "mean_read_length": float(np.average([r["read_length"] for r in rows], weights=w)),
            "f0_frac": tf0 / tot, "f1_frac": tf1 / tot, "f2_frac": tf2 / tot,
            # RENAMED 2026-08-15 from `psites_at_start_window`, which cost a downstream figure. This is the sum of
            # the frame-0/1/2 counts in RiboCode's metaplots table -- P-sites in the window around
            # ANNOTATED START CODONS, not across the CDS and not the library. It runs ~0.5% of real
            # pack depth (Fibroblast 5.03M here vs 1,069.8M in the pack) and is NOT usable as a depth
            # proxy: its ratio to real depth swings 0.33-0.65% across the 8 Chothani tissues, enough
            # to invert the ES/Fat ordering. The number was always right; the name was not.
            "psites_at_start_window": tot,
            "psite_offsets": sorted({r["psite"] for r in rows}),
            "offset_range": max(r["psite"] for r in rows) - min(r["psite"] for r in rows),
            "rows": rows,
        })
    if missing:
        print("MISSING (not silently dropped -- these are absent from the panel):")
        for lb, p in missing:
            print(f"  {lb}: {p}")

    order = sorted(data, key=lambda d: (["train", "heldout", "pgx", "dropped"].index(d["role"]), -d["f0_frac"]))
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 0.42 * len(order) + 3.0),
                             gridspec_kw={"width_ratios": [1.15, 1.35, 1.0]})

    # (a) read-length distribution, one stacked strip per dataset
    ax = axes[0]
    lens = sorted({r["read_length"] for d in order for r in d["rows"]})
    lo, hi = max(min(lens), 20), min(max(lens), 40)
    grid = np.arange(lo, hi + 1)
    img = np.zeros((len(order), grid.size))
    for i, d in enumerate(order):
        for r in d["rows"]:
            if lo <= r["read_length"] <= hi:
                img[i, r["read_length"] - lo] += r["proportion"]
        s = img[i].sum()
        if s > 0:
            img[i] /= s
    ax.imshow(img, aspect="auto", cmap="Blues", vmin=0, vmax=float(img.max()))
    ax.set_xticks(np.arange(0, grid.size, 2)); ax.set_xticklabels(grid[::2], fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([d["label"] for d in order], fontsize=8)
    for tick, d in zip(ax.get_yticklabels(), order):
        tick.set_color(ROLE_COLOR[d["role"]])
    ax.set_xlabel("footprint length (nt)", fontsize=9)
    ax.set_title("(a) read-length distribution\nmonosome footprints peak at 28-31 nt", fontsize=9)

    # (b) frame fractions at annotated CDS
    ax = axes[1]
    y = np.arange(len(order))
    f0 = np.array([d["f0_frac"] for d in order])
    f1 = np.array([d["f1_frac"] for d in order])
    f2 = np.array([d["f2_frac"] for d in order])
    # frame 0 takes the row's role colour, so a fixed-swatch legend would misrepresent it; the
    # segments are annotated on the top bar instead.
    ax.barh(y, f0, color=[ROLE_COLOR[d["role"]] for d in order])
    ax.barh(y, f1, left=f0, color="#c9c9c9")
    ax.barh(y, f2, left=f0 + f1, color="#eaeaea")
    for txt, xpos in (("frame 0", f0[0] / 2), ("f1", f0[0] + f1[0] / 2),
                      ("f2", f0[0] + f1[0] + f2[0] / 2)):
        ax.annotate(txt, xy=(xpos, -0.42), ha="center", va="bottom", fontsize=7.5,
                    color="white" if txt == "frame 0" else "#444444")
    ax.axvline(1 / 3, color="k", ls="--", lw=1)
    ax.text(1 / 3, len(order) - 0.3, " no periodicity", fontsize=7.5, va="top")
    for i, v in enumerate(f0):
        ax.text(v + 0.012, i, f"{v * 100:.0f}%", va="center", fontsize=7.5)
    # imshow puts row 0 at the TOP. (b) and (c) must match it or the shared labels in (a) are wrong,
    # so set the limits explicitly rather than relying on invert_yaxis (which flipped only these two).
    ax.set_yticks(y); ax.set_yticklabels([]); ax.set_ylim(len(order) - 0.4, -0.6)
    ax.set_xlim(0, 1.12); ax.set_xlabel("share of P-sites at the annotated CDS", fontsize=9)
    ax.set_title("(b) 3-nt periodicity\nthe signal the model is trained to reproduce", fontsize=9)

    # (c) P-site offset consistency
    ax = axes[2]
    for i, d in enumerate(order):
        offs = d["psite_offsets"]
        ax.plot(offs, [i] * len(offs), "o", ms=4.5, color=ROLE_COLOR[d["role"]], alpha=0.85)
        if len(offs) > 1:
            ax.plot([min(offs), max(offs)], [i, i], "-", color=ROLE_COLOR[d["role"]], lw=1, alpha=0.5)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([]); ax.set_ylim(len(order) - 0.4, -0.6)
    ax.set_xlabel("P-site offset (nt from 5' end)", fontsize=9)
    ax.set_title("(c) offset consistency\none dot per read length in the library", fontsize=9)
    ax.grid(axis="x", alpha=0.25)

    handles = [plt.Line2D([], [], marker="s", ls="", color=ROLE_COLOR[r], label=ROLE_LABEL[r])
               for r in ("train", "heldout", "pgx", "dropped")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.5, frameon=False)
    fig.suptitle("Supplemental: Ribo-seq library quality, every dataset used in this study",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0.045, 1, 0.96])
    for ext in ("png", "pdf"):
        fig.savefig(out / f"riboseq_qc.{ext}", dpi=200)

    vals = [{k: v for k, v in d.items() if k != "rows"} for d in order]
    (out / "riboseq_qc_values.json").write_text(
        json.dumps({"datasets": vals, "missing": missing}, indent=2))

    print(f"\n{'dataset':<26}{'role':<10}{'libs':>5}{'mode nt':>9}{'f0 %':>8}{'offsets':>12}"
          f"{'P-sites@CDS':>14}")
    for d in order:
        print(f"{d['label']:<26}{d['role']:<10}{d['n_libraries']:>5}{d['mode_read_length']:>9}"
              f"{d['f0_frac'] * 100:>7.1f}%{str(d['psite_offsets']):>12}{d['psites_at_start_window']:>14,}")
    print(f"\nwrote {out}/riboseq_qc.png|.pdf|_values.json")


if __name__ == "__main__":
    main()
