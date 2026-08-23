#!/usr/bin/env python3
"""Assemble main Figures 1 and 2 from the per-panel PDFs (task 64).

VECTOR-PRESERVING. Each panel is placed as its original PDF page, scaled and translated, not
re-rastered -- so text stays selectable and lines stay crisp at any zoom. Panel letters are drawn into
a matplotlib overlay PDF of the same page size and merged on top. Requires `pypdf` (installed into the
cas12a env 2026-08-09).

WHY A LAYOUT SCRIPT RATHER THAN A HAND-ASSEMBLED FILE. Every panel here is regenerated from result
JSONs whenever the model changes, and five of six Fig 1 panels were recently found to have been built
from a superseded checkpoint. A figure assembled once by hand goes stale the same way and cannot be
re-derived; this can be re-run after any retrain, and it records the exact source PDF and its mtime for
every panel in `<fig>_manifest.json`.

The layout is a plain data structure (FIGURES below): a list of rows, each row a list of
(panel_path, relative_width). Re-specifying a figure means editing that list, nothing else.

Usage:  make_main_figures.py [--width IN] [--only fig1|fig2]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pypdf import PageObject, PdfReader, PdfWriter, Transformation

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
FIGDIR = NEW / "figures"
HERE = FIGDIR / "main"
PT = 72.0
GAP = 0.14          # inches between panels, both axes
MARGIN = 0.10       # inches around the whole page

# Each figure: rows of (relative path under figures/, relative width within the row).
# Relative widths inside a row are normalised, so [1, 1] is a half/half split.
FIGURES = {
    "Figure1_validity": {
        "caption": "the model predicts translation and generalizes",
        "rows": [
            [("prediction_examples/prediction_examples.pdf", 1.0)],
            [("A2_ribocode_dropin/A2_ribocode_dropin.pdf", 1.0),
             ("A1_localization_ceiling/A1_localization_ceiling.pdf", 0.82)],
            # Equal widths here on purpose. depth_crossover is nearly square (8.4 x 8.0 in), so
            # giving it the narrower half of an unequal split shrank it to 2.6 in and its axis text
            # became unreadable. A3 is wide-and-short and tolerates the same width fine.
            [("depth_crossover/depth_crossover.pdf", 1.0),
             ("A3_loto_spread/A3_loto_spread.pdf", 1.0)],
        ],
    },
    "Figure2_utility": {
        "caption": "the immunopeptidome application",
        "rows": [
            [("F2b_discovery_forest/F2b_discovery_forest.pdf", 1.0)],
            [("D12_cpat_cpc2/D12_cpat_cpc2.pdf", 1.0)],
            [("D13_discovery_errorbars/D13_discovery_errorbars.pdf", 1.0)],
        ],
    },
}
LETTERS = "abcdefghijklmnop"


def page_size_in(pdf: Path):
    b = PdfReader(str(pdf)).pages[0].mediabox
    return float(b.width) / PT, float(b.height) / PT


def plan(rows, page_w):
    """Resolve rows -> placements in inches, origin TOP-left (converted to PDF coords later)."""
    inner_w = page_w - 2 * MARGIN
    placed, y = [], MARGIN
    for row in rows:
        n = len(row)
        avail = inner_w - GAP * (n - 1)
        tot = sum(rw for _, rw in row)
        row_h, cells = 0.0, []
        x = MARGIN
        for rel, rw in row:
            src = FIGDIR / rel
            if not src.exists():
                raise SystemExit(f"missing panel: {src}")
            w_in, h_in = page_size_in(src)
            target_w = avail * (rw / tot)
            s = target_w / w_in
            cells.append({"src": src, "x": x, "w": target_w, "h": h_in * s, "scale": s})
            x += target_w + GAP
            row_h = max(row_h, h_in * s)
        for c in cells:
            c["y"] = y            # top of the row; panels are top-aligned within a row
        placed.extend(cells)
        y += row_h + GAP
    return placed, y - GAP + MARGIN


def letter_overlay(placements, page_w, page_h, out):
    """A same-size PDF containing only the panel letters, merged over the panels."""
    fig = plt.figure(figsize=(page_w, page_h))
    for i, c in enumerate(placements):
        # figure coords are fraction of page, origin bottom-left
        fx = (c["x"] - 0.02) / page_w
        fy = 1.0 - (c["y"] - 0.02) / page_h
        fig.text(max(fx, 0.002), min(fy, 0.998), LETTERS[i], fontsize=13, fontweight="bold",
                 va="top", ha="left", family="sans-serif")
    fig.savefig(out, format="pdf", transparent=True)
    plt.close(fig)


def build(name, spec, page_w):
    placements, page_h = plan(spec["rows"], page_w)
    page = PageObject.create_blank_page(width=page_w * PT, height=page_h * PT)
    for c in placements:
        src = PdfReader(str(c["src"])).pages[0]
        # PDF origin is BOTTOM-left; `y` above is measured from the TOP.
        ty = (page_h - c["y"] - c["h"]) * PT
        op = Transformation().scale(c["scale"]).translate(c["x"] * PT, ty)
        page.merge_transformed_page(src, op)

    ov = HERE / f".{name}_letters.pdf"
    letter_overlay(placements, page_w, page_h, ov)
    page.merge_page(PdfReader(str(ov)).pages[0])
    ov.unlink()

    w = PdfWriter(); w.add_page(page)
    out = HERE / f"{name}.pdf"
    with open(out, "wb") as fh:
        w.write(fh)

    man = {
        "figure": name, "caption": spec["caption"],
        "page_inches": [round(page_w, 3), round(page_h, 3)],
        "panels": [{"letter": LETTERS[i], "source": str(c["src"].relative_to(NEW)),
                    "source_mtime": _dt.datetime.fromtimestamp(
                        c["src"].stat().st_mtime).isoformat(timespec="seconds"),
                    "placed_inches": {"x": round(c["x"], 3), "y_from_top": round(c["y"], 3),
                                      "w": round(c["w"], 3), "h": round(c["h"], 3)},
                    "scale": round(c["scale"], 4)}
                   for i, c in enumerate(placements)],
    }
    (HERE / f"{name}_manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    print(f"{name}: {page_w:.2f} x {page_h:.2f} in, {len(placements)} panels -> {out.name}")
    for i, c in enumerate(placements):
        print(f"   ({LETTERS[i]}) {c['src'].relative_to(FIGDIR)}  "
              f"{c['w']:.2f}x{c['h']:.2f} in  scale {c['scale']:.3f}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=float, default=7.2,
                    help="page width in inches (default 7.2 = 183 mm, full journal width)")
    ap.add_argument("--only", choices=sorted(FIGURES))
    a = ap.parse_args()
    HERE.mkdir(parents=True, exist_ok=True)
    for name, spec in FIGURES.items():
        if a.only and name != a.only:
            continue
        build(name, spec, a.width)
    return 0


if __name__ == "__main__":
    sys.exit(main())
