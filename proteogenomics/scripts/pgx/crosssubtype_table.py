#!/usr/bin/env python3
"""Pivot the per-population pgx reports into one cross-subtype table.

Each population is searched against the SAME BMDM-derived databases (real Ribo-seq NT-only, real
Ribo-seq NT+LPS) plus its own model-predicted database, so the shared arms are directly comparable
across rows and any difference between subtypes is biological rather than bookkeeping.

The null arms are the exception and are labelled as such: each population's nulls derive from its
OWN expressed universe, so their database sizes differ per row by construction. Compare nulls down
a column with care, and never treat a null count as evidence of tissue specificity.

  python -m pgx.crosssubtype_table --report-dir <dir> --out table.md
"""
from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

ARMS = [
    ("model_predicted", "model predicted"),
    ("model_ribocode_bmdm_nt", "real Ribo-seq NT"),
    ("model_ribocode_bmdm", "real Ribo-seq NT+LPS"),
    ("null_atg", "null ATG"),
    ("null_nc", "null near-cognate"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-dir", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    data = {}
    for f in sorted(glob.glob(str(Path(a.report_dir) / "report_*.json"))):
        d = json.loads(Path(f).read_text())
        data[d["label"]] = d

    present = [(k, lab) for k, lab in ARMS if any(k in d["arms"] for d in data.values())]
    L = []

    for metric, title, fmt in [
        ("novel_peptides", "Novel peptides (unique, 1% class-specific FDR)", "{:,}"),
        ("novel_psms", "Novel PSMs (total)", "{:,}"),
        ("db_novel_seqs", "Database novel sequences", "{:,}"),
        ("d_gencode_psms", "Canonical cost, dPSM vs GENCODE-only baseline", "{:+,}"),
    ]:
        L += ["", f"### {title}", "",
              "| population | " + " | ".join(lab for _k, lab in present) + " |",
              "|---" + "|--:" * len(present) + "|"]
        for pop in sorted(data):
            cells = []
            for k, _lab in present:
                v = data[pop]["arms"].get(k, {}).get(metric)
                cells.append("-" if v is None else fmt.format(v))
            L.append(f"| {pop} | " + " | ".join(cells) + " |")
        med = []
        for k, _lab in present:
            vals = [d["arms"][k][metric] for d in data.values() if k in d["arms"]]
            med.append(fmt.format(int(st.median(vals))) if vals else "-")
        L.append("| **median** | " + " | ".join(f"**{m}**" for m in med) + " |")

    # discovery density: peptides per 1,000 database sequences, the size-normalised comparison
    L += ["", "### Discovery density (novel peptides per 1,000 database sequences)", "",
          "| population | " + " | ".join(lab for _k, lab in present) + " |",
          "|---" + "|--:" * len(present) + "|"]
    dens = {k: [] for k, _ in present}
    for pop in sorted(data):
        cells = []
        for k, _lab in present:
            arm = data[pop]["arms"].get(k, {})
            n, s = arm.get("novel_peptides"), arm.get("db_novel_seqs")
            if not n or not s:
                cells.append("-")
            else:
                v = n / s * 1000
                dens[k].append(v)
                cells.append(f"{v:.1f}")
        L.append(f"| {pop} | " + " | ".join(cells) + " |")
    L.append("| **median** | " + " | ".join(
        f"**{st.median(dens[k]):.1f}**" if dens[k] else "-" for k, _ in present) + " |")

    L += ["",
          "All arms are 1% FDR filtered: novel columns class-specific (vs `REV_nuORF|` decoys "
          "only), canonical columns global. PSM columns are total rank-1 spectra; peptide columns "
          "are distinct sequences.",
          "The two real-Ribo-seq arms and the null arms differ in one respect that matters: the "
          "Ribo-seq databases are BMDM-derived and shared by every row, so they are comparable "
          "across subtypes; each population's NULLS come from its own expressed universe, so null "
          "database sizes differ per row and the null columns are not cross-comparable.",
          "The real-Ribo-seq arms contain no N-terminal extensions, because RiboCode reports them "
          "without testing whether the upstream start is used. They are ATG-only by construction, "
          "so their ncStart is structurally zero and that is a property of the caller, not "
          "evidence about non-AUG initiation.",
          "BMDM is the MATCHED case for the Ribo-seq arms (the database came from BMDM Ribo-seq); "
          "the other 11 rows are the transfer test."]
    md = "\n".join(L) + "\n"
    print(md)
    if a.out:
        Path(a.out).write_text(md)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
