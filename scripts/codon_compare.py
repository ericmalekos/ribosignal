#!/usr/bin/env python3
"""Phase 4.2-4.4: is the model's codon occupancy real, and does it generalise?

Three questions, in the only order that makes them answerable:

  4.2  THE EMPIRICAL CEILING FIRST. Correlate the OBSERVED 61-codon vectors between datasets. Codon
       dwell is protocol-sensitive (cycloheximide artefacts, digestion, depth), so a model-vs-observed
       correlation of, say, 0.22 means something completely different depending on whether two real
       experiments agree at 0.95 or at 0.30. Without this number the model comparison is
       uninterpretable -- the same logic as the ORF-call between-experiment ceiling.

  4.3  MODEL vs OBSERVED, per dataset. The model never receives codon identity, amino acid or tRNA
       abundance as input; only one-hot sequence, an ORF track and RNA coverage. Any recovered codon
       structure was learned from sequence.

  4.4  DOES IT GENERALISE. Recovering codon dwell on a training tissue is partly just fitting the
       training data. Recovering it on MOUSE -- different species, no transcript ever seen -- is
       evidence of a learned rule. Reported separately and never averaged.

Reads the per-dump tables from codon_occupancy.py.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import pathlib
from collections import defaultdict

import numpy as np


def load(dirpath, site):
    """{label: {source: 61-vector}} for one site convention."""
    out = defaultdict(dict)
    for f in sorted(pathlib.Path(dirpath).glob("*.tsv")):
        rows = list(csv.DictReader(open(f), delimiter="\t"))
        for src in ("observed", "predicted"):
            sel = {r["codon"]: float(r["occupancy"]) for r in rows
                   if r["source"] == src and r["site"] == site
                   and r["occupancy"] not in ("nan", "")}
            if len(sel) >= 55:
                out[f.stem][src] = sel
    return out


def vec(d, codons):
    return np.array([d.get(c, np.nan) for c in codons])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--site", default="asite", choices=["asite", "psite"])
    a = ap.parse_args()

    data = load(a.dir, a.site)
    labels = sorted(data)
    codons = sorted({c for lab in labels for src in data[lab] for c in data[lab][src]})
    print(f"  {len(labels)} dumps, {len(codons)} codons, site={a.site}\n")

    def r(x, y):
        ok = np.isfinite(x) & np.isfinite(y)
        return float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 10 else float("nan")

    species = lambda l: "mouse" if l.startswith(("mouse", "released")) else "human"  # noqa: E731

    def biodataset(lab):
        """The underlying Ribo-seq dataset, independent of model and universe.

        Two dumps of the SAME dataset share one `obs_flat`, so their observed 61-vectors are
        identical by construction and correlate at exactly 1.000. Including such pairs in the
        empirical ceiling would inflate it with a tautology -- the ceiling must be measured between
        INDEPENDENT experiments.
        """
        for k in ("janich", "gse243134", "wang", "cart", "gse208041"):
            if k in lab:
                return k
        return lab

    # ---- 4.2 empirical ceiling: observed vs observed, between datasets ----
    print("  === 4.2 EMPIRICAL CEILING (observed vs observed) ===")
    ceil = []
    n_self = 0
    for x, y in itertools.combinations(labels, 2):
        if "observed" not in data[x] or "observed" not in data[y]:
            continue
        if biodataset(x) == biodataset(y):
            n_self += 1          # same experiment seen twice; r == 1.000 by construction
            continue
        v = r(vec(data[x]["observed"], codons), vec(data[y]["observed"], codons))
        ceil.append(dict(a=x, b=y, r=v, same_species=species(x) == species(y)))
    print(f"    (excluded {n_self} same-dataset pairs -- identical obs_flat, r=1.000 by construction)")
    within = [c["r"] for c in ceil if c["same_species"]]
    across = [c["r"] for c in ceil if not c["same_species"]]
    if within:
        print(f"    within-species  n={len(within):>3}  r = {np.median(within):.3f} "
              f"(range {min(within):.3f} to {max(within):.3f})")
    if across:
        print(f"    across-species  n={len(across):>3}  r = {np.median(across):.3f} "
              f"(range {min(across):.3f} to {max(across):.3f})")

    # ---- 4.3 / 4.4 model vs observed, split by whether transcripts were ever seen ----
    print("\n  === 4.3/4.4 MODEL vs OBSERVED (per dump) ===")
    mv = []
    for lab in labels:
        if "observed" not in data[lab] or "predicted" not in data[lab]:
            continue
        v = r(vec(data[lab]["observed"], codons), vec(data[lab]["predicted"], codons))
        sp = species(lab)
        mv.append(dict(label=lab, species=sp, r=v))
        print(f"    {lab:<32} {sp:<6} r = {v:+.3f}")
    m_mouse = [x["r"] for x in mv if x["species"] == "mouse"]
    m_human = [x["r"] for x in mv if x["species"] == "human"]
    print()
    if m_mouse:
        print(f"    MOUSE  (transcripts never seen)  median r = {np.median(m_mouse):+.3f}  n={len(m_mouse)}")
    if m_human:
        print(f"    HUMAN  (99.7% seen in training)  median r = {np.median(m_human):+.3f}  n={len(m_human)}")

    res = dict(site=a.site, n_dumps=len(labels), ceiling=ceil, model_vs_observed=mv,
               ceiling_within_species_median=float(np.median(within)) if within else None,
               ceiling_across_species_median=float(np.median(across)) if across else None,
               model_mouse_median=float(np.median(m_mouse)) if m_mouse else None,
               model_human_median=float(np.median(m_human)) if m_human else None)
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
