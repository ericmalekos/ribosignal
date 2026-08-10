#!/usr/bin/env python3
"""Write a plain TSV of the plotted values next to every figure.

Purpose: the figures are going to be restyled in other tools, and a PDF is not a data source. Every
generator already writes a `*_values.json` holding exactly the numbers it plotted; this converts each
of those into a TSV in the same folder, so the numbers can be read without parsing nested JSON or
re-running any model.

ONE CONVERTER RATHER THAN EIGHTEEN EDITS. Adding a bespoke `to_tsv` to each generator would mean
eighteen chances to emit something subtly different from what was drawn. This reads the same JSON the
generator wrote after drawing, so the TSV cannot drift from the figure. New figures are picked up
automatically -- nothing to register.

Two output shapes, chosen per file:
  WIDE  when the JSON's substance is a record table (a list of dicts, or a dict of same-shaped dicts):
        one row per record, one column per field, with any surrounding scalars repeated as context
        columns. This is the friendly form and covers most figures.
  LONG  otherwise: one row per leaf value, as `path<TAB>value`, with the dotted key path preserved.
        Universal and lossless; pivot it if you want a matrix.

The shape used is printed and recorded in the header comment of each TSV.

Usage:  export_figure_tsvs.py [--figdir DIR] [--quiet]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

FIGDIR = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
              "RNAZoo/experiments/riboseq_signal_model/figures")
SKIP_SUBSTR = ("manifest",)          # main-figure manifests are provenance, not plotted data


def is_record_list(v):
    return (isinstance(v, list) and len(v) >= 2 and all(isinstance(e, dict) for e in v))


def is_record_dict(v):
    """dict of same-shaped dicts, e.g. {"uORF": {"recall": .., "n_real": ..}, ...}"""
    if not isinstance(v, dict) or len(v) < 2:
        return False
    vals = list(v.values())
    if not all(isinstance(e, dict) for e in vals):
        return False
    keysets = [frozenset(e) for e in vals]
    # allow ragged records but require real overlap, else it is just a nested namespace
    common = set.intersection(*[set(k) for k in keysets])
    return len(common) >= 1


def scalar(v):
    return v is None or isinstance(v, (str, int, float, bool))


def find_table(obj):
    """Return (path, rows) for the largest record table anywhere in the JSON, else None."""
    best = None

    def walk(o, path):
        nonlocal best
        if is_record_list(o):
            rows = [dict(e) for e in o]
            if best is None or len(rows) > len(best[1]):
                best = (path, rows)
        elif is_record_dict(o):
            rows = []
            for k, e in o.items():
                r = {"key": k}
                r.update(e)
                rows.append(r)
            if best is None or len(rows) > len(best[1]):
                best = (path, rows)
        if isinstance(o, dict):
            for k, v in o.items():
                walk(v, f"{path}.{k}" if path else k)
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]")

    walk(obj, "")
    return best


def flatten(o, path="", out=None):
    out = [] if out is None else out
    if isinstance(o, dict):
        for k, v in o.items():
            flatten(v, f"{path}.{k}" if path else k, out)
    elif isinstance(o, list):
        if all(scalar(e) for e in o):
            out.append((path, ";".join("" if e is None else str(e) for e in o)))
        else:
            for i, v in enumerate(o):
                flatten(v, f"{path}[{i}]", out)
    else:
        out.append((path, o))
    return out


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:.8g}"
    if isinstance(v, (dict, list)):
        return json.dumps(v, separators=(",", ":"))
    return str(v)


def context_scalars(obj, table_path):
    """Top-level scalars that describe the whole file (model, n, recipe...). Repeated on every row."""
    if not isinstance(obj, dict):
        return {}
    head = table_path.split(".")[0] if table_path else None
    return {k: v for k, v in obj.items() if scalar(v) and k != head}


def convert(src: Path, quiet=False):
    try:
        obj = json.loads(src.read_text())
    except Exception as e:
        print(f"  SKIP {src.name}: unreadable ({e})", file=sys.stderr)
        return None
    dst = src.with_suffix(".tsv")
    tab = find_table(obj)

    if tab and len(tab[1]) >= 2:
        path, rows = tab
        ctx = context_scalars(obj, path)
        cols, seen = [], set()
        for r in rows:
            for k in r:
                if k not in seen:
                    seen.add(k); cols.append(k)
        ctx_cols = [k for k in ctx if k not in seen]
        with dst.open("w") as fh:
            fh.write(f"# {src.name} -> WIDE (record table at '{path or '<root>'}', "
                     f"{len(rows)} rows)\n")
            fh.write("\t".join(ctx_cols + cols) + "\n")
            for r in rows:
                fh.write("\t".join([fmt(ctx[c]) for c in ctx_cols]
                                   + [fmt(r.get(c)) for c in cols]) + "\n")
        shape = f"WIDE {len(rows)}x{len(ctx_cols) + len(cols)}"
    else:
        leaves = flatten(obj)
        with dst.open("w") as fh:
            fh.write(f"# {src.name} -> LONG (no record table; one row per value)\n")
            fh.write("path\tvalue\n")
            for k, v in leaves:
                fh.write(f"{k}\t{fmt(v)}\n")
        shape = f"LONG {len(leaves)} rows"
    if not quiet:
        print(f"  {src.parent.name:<28}{src.name:<34}-> {dst.name:<34}{shape}")
    return dst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figdir", default=str(FIGDIR))
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    root = Path(a.figdir)
    srcs = sorted(p for p in root.glob("*/*.json")
                  if not any(s in p.name.lower() for s in SKIP_SUBSTR))
    n = sum(1 for s in srcs if convert(s, a.quiet))
    print(f"\nwrote {n} TSV(s) from {len(srcs)} values JSON(s) under {root}")

    # Figures whose plotted data is NOT fully captured by a values JSON deserve a name-check rather
    # than silence: a missing TSV that nobody notices is the same failure mode as a stale number.
    missing = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        if d.name == "main":
            continue
        if not any(d.glob("*.pdf")) and not any(d.glob("*.png")):
            continue
        if not any(d.glob("*.tsv")):
            missing.append(d.name)
    if missing:
        print("\nFigures with output but still NO tsv (their generator writes no values JSON):")
        for m in missing:
            print(f"  {m}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
