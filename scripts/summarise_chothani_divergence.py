#!/usr/bin/env python3
"""Collect the per-tissue canonical divergence JSONs into one table.

Turns the reproducibility claim from one tissue into eight. Emits a markdown table for results.md
and a TSV for figures, and reports the depth-vs-divergence ordering that separates the regeneration
noise floor from the recipe term (see results.md, "Regeneration divergence separates into a depth
noise floor and a recipe term").

Divergence is `sum|new - old| / sum(old)` over per-nt P-sites on each pack's own universe -- the same
metric as the original HUVEC test, which this script reproduces exactly (3.04%, r=0.9567) when run
over the HUVEC output. Net total change is reported alongside it because the two answer different
questions: net says how much signal was gained or lost, |divergence| says how much MOVED.
"""
import argparse
import csv
import json
import pathlib
import statistics
import sys

TISSUES = ["ES", "Fat", "Fibroblast", "HA_EC", "HCAEC", "Hepatocytes", "HUVEC", "VSMC"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True,
                    help="results/chothani_regeneration (holds _canonical_<T>/divergence_canonical.json)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-md", required=True)
    a = ap.parse_args()

    root = pathlib.Path(a.dir)
    rows, missing = [], []
    for t in TISSUES:
        f = root / f"_canonical_{t}" / "divergence_canonical.json"
        if not f.exists():
            missing.append(t)
            continue
        rows.append(json.loads(f.read_text()))

    if not rows:
        sys.exit(f"no divergence JSONs under {root}")
    # Report what is absent rather than quietly tabulating a subset: a partial table that looks
    # complete is how a coverage gap turns into a claim (no_silent_caps).
    if missing:
        print(f"  INCOMPLETE: {len(missing)} tissue(s) missing -> {', '.join(missing)}", file=sys.stderr)

    rows.sort(key=lambda r: -r["old_psites"])
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    lines = [
        "| tissue | pack | n bam | training P-sites | canonical | net | \\|divergence\\| | per-nt r |",
        "|---|---|--:|--:|--:|--:|--:|--:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['tissue']} | `{r['pack']}` | {r['n_bams']} | {r['old_psites']:,} | "
            f"{r['new_psites']:,} | {r['pct_total_change']:+.2f}% | **{r['pct_of_signal']:.2f}%** | "
            f"{r['pearson']:.4f} |")
    if missing:
        lines.append("")
        lines.append(f"*Missing: {', '.join(missing)} -- table covers {len(rows)} of {len(TISSUES)} tissues.*")
    pathlib.Path(a.out_md).write_text("\n".join(lines) + "\n")

    div = [r["pct_of_signal"] for r in rows]
    # statistics.median, NOT sorted(div)[len(div)//2]: with an even count the latter returns the
    # UPPER of the two middle values (4.01 instead of 3.86 on the 8-tissue set) and that wrong
    # number was briefly propagated into results.md and the figure caption.
    print(f"  {len(rows)} tissue(s): divergence median {statistics.median(div):.2f}% "
          f"range {min(div):.2f}-{max(div):.2f}%")
    print(f"  per-nt r range {min(r['pearson'] for r in rows):.4f}-{max(r['pearson'] for r in rows):.4f}")
    print("  depth vs divergence (deepest first):")
    for r in rows:
        print(f"    {r['tissue']:<12} {r['old_psites']:>13,}  {r['pct_of_signal']:>6.2f}%  r={r['pearson']:.4f}")
    print(f"  wrote {a.out_tsv} and {a.out_md}")


if __name__ == "__main__":
    main()
