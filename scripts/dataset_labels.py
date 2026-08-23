#!/usr/bin/env python3
"""Resolve a dataset's internal key to its display label. Import this; never hardcode a label.

    from dataset_labels import display, info
    display("janich")                 -> "mouse_liver_GSE67305"
    display("Hepatocytes")            -> "human_hepatocytes_GSE182371"
    display("gse243134", short=True)  -> "liver GSE243134"

WHY. Labels grew as a mix of surnames, GEO accessions, cell types and product names, and the
surname-labelled ones were surname-labelled because no accession had been recorded anywhere. Two
datasets are also both THP-1 (GSE208041, GSE39561), so a cell-type label is ambiguous in our own
data and a surname breaks on GSE243134, whose paper has multiple first authors.

INTERNAL KEYS ARE NOT RENAMED. They are load-bearing in directory names, pack paths, script
arguments and every values JSON. This maps key -> label at RENDER time, which gets consistency
without a rename whose blast radius spans the whole project.

CHOTHANI IS PER-TISSUE. One study (GSE182371) supplied 9 cell types, so its registry row carries a
`<tissue>` placeholder and this module substitutes the tissue key: `Hepatocytes` ->
`human_hepatocytes_GSE182371`. Passing a bare Chothani tissue name works directly.

Run as a script to print the full cross-reference table.
"""
from __future__ import annotations

import csv
import functools
import pathlib

REGISTRY = pathlib.Path(__file__).resolve().parents[1] / "data" / "dataset_registry.tsv"

# The 9 Chothani cell types. Brain is included because it EXISTS in the source study; it was dropped
# from training for low periodicity, which is a modelling decision, not an identity question.
CHOTHANI_TISSUES = ("Brain", "ES", "Fat", "Fibroblast", "HA_EC", "HCAEC",
                    "Hepatocytes", "HUVEC", "VSMC")


@functools.lru_cache(maxsize=1)
def _rows():
    with open(REGISTRY) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


@functools.lru_cache(maxsize=1)
def _by_key():
    return {r["key"]: r for r in _rows()}


def info(key):
    """-> the registry row for a dataset key, or for a Chothani tissue. KeyError if unknown."""
    reg = _by_key()
    if key in reg:
        return reg[key]
    # Chothani tissues are not their own rows -- they share one accession.
    for t in CHOTHANI_TISSUES:
        if key.lower() == t.lower():
            r = dict(reg["chothani"])
            r["tissue"] = t
            r["display"] = f"human_{t.lower()}_{r['accession_ribo']}"
            return r
    raise KeyError(f"{key!r} is not in {REGISTRY.name} and is not a Chothani tissue. "
                   f"Add a row rather than hardcoding a label.")


def display(key, short=False):
    """Canonical display label. short=True drops the species prefix, for crowded axes."""
    r = info(key)
    d = r["display"]
    if short:
        parts = d.split("_", 1)
        return parts[1].replace("_", " ") if len(parts) > 1 else d
    return d


def accession(key):
    """Primary accession: the Ribo-seq series where one exists, else the SRA/PX project."""
    r = info(key)
    for f in ("accession_ribo", "sra"):
        if r.get(f, "-").strip() not in ("-", ""):
            return r[f]
    return "-"


def main():
    rows = _rows()
    w = max(len(r["key"]) for r in rows) + 2
    print(f"{'INTERNAL KEY':<{w}}{'DISPLAY LABEL':<34}{'ACCESSION':<14}{'VERIFIED'}")
    print("-" * (w + 60))
    for r in sorted(rows, key=lambda x: (x["species"], x["tissue"])):
        print(f"{r['key']:<{w}}{r['display']:<34}{accession(r['key']):<14}{r['verified']}")
    print(f"\nChothani (GSE182371) expands per tissue:")
    for t in CHOTHANI_TISSUES:
        note = "  [dropped from training: low periodicity]" if t == "Brain" else ""
        print(f"  {t:<{w}}{display(t)}{note}")


if __name__ == "__main__":
    main()
