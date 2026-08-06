#!/usr/bin/env python3
"""Fraction-level replication check for novel peptides (`pgx` credibility diagnostic).

A target-decoy FDR is an aggregate statement about a list; it does not say whether an INDIVIDUAL
identification is trustworthy. That matters most exactly where the FDR model is least reliable: a
search space of millions of heavily nested near-cognate ORFs, whose target and decoy sides have
very different redundancy structure than the procedure assumes.

Replication across fractions is an orthogonal check that does not rely on the decoy model at all.
A genuinely present peptide is usually sampled in more than one fraction or replicate; a chance
match to one spectrum in a very large search space is usually seen once and never again. So the
share of FDR-passing novel peptides supported by a single fraction is a direct, assumption-light
credibility signal, and it is comparable across arms because the spectra are identical -- only the
database differs.

Read it as a comparison, never as an absolute threshold: low-abundance real peptides do legitimately
appear once, so a single-fraction fraction is only interpretable RELATIVE to another arm searched
over the same spectra.

  python -m pgx.replication --search-root <dir> --arms model_poisson,null_atg,null_nc
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import sys
from pathlib import Path

from .report import FDR, klass, load_rank1, novel_at_class_fdr, split_accs


def fractions_per_peptide(search_dir, keep):
    """{peptide: set(fraction file stems)} over novel rank-1 PSMs for peptides in `keep`."""
    out = collections.defaultdict(set)
    for t in sorted(glob.glob(str(Path(search_dir) / "*.tsv"))):
        frac = Path(t).stem
        with open(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("hit_rank") != "1":
                    continue
                p = r["peptide"]
                if p in keep and klass(split_accs(r.get("proteins", ""))) == "novel_t":
                    out[p].add(frac)
    return out


def summarize(search_dir, fdr=FDR):
    rows = load_rank1(search_dir)
    if not rows:
        return None
    n_psm, peps, accs, _ = novel_at_class_fdr(rows, fdr)
    if not peps:
        return {"novel_peptides": 0}
    pf = fractions_per_peptide(search_dir, peps)
    counts = collections.Counter(len(v) for v in pf.values())
    tot = sum(counts.values())
    return {
        "novel_peptides": tot,
        "novel_psms": n_psm,
        "novel_seqs": len(accs),
        "one_fraction": counts.get(1, 0),
        "ge2_fractions": sum(v for k, v in counts.items() if k >= 2),
        "ge3_fractions": sum(v for k, v in counts.items() if k >= 3),
        "max_fractions": max(counts) if counts else 0,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search-root", required=True)
    ap.add_argument("--arms", default=None, help="comma list; default = every subdirectory")
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    root = Path(a.search_root)
    arms = ([x.strip() for x in a.arms.split(",") if x.strip()] if a.arms
            else sorted(p.name for p in root.iterdir() if p.is_dir()))
    lines = ["| arm | novel pept | 1 fraction | >=2 fractions | >=3 fractions | max |",
             "|---|--:|--:|--:|--:|--:|"]
    for arm in arms:
        s = summarize(root / arm, a.fdr)
        if not s or not s["novel_peptides"]:
            print(f"[skip] {arm}: no novel peptides", file=sys.stderr)
            continue
        t = s["novel_peptides"]
        lines.append(f"| {arm} | {t} | {s['one_fraction']} ({s['one_fraction'] / t * 100:.1f}%) | "
                     f"{s['ge2_fractions']} ({s['ge2_fractions'] / t * 100:.1f}%) | "
                     f"{s['ge3_fractions']} ({s['ge3_fractions'] / t * 100:.1f}%) | "
                     f"{s['max_fractions']} |")
        print(f"  {arm} done", file=sys.stderr)
    lines += ["",
              "Share of FDR-passing novel peptides seen in only ONE fraction. Spectra are identical",
              "across arms, so this is comparable between them and does not depend on the decoy",
              "model. A high single-fraction share signals identifications that the aggregate FDR",
              "admitted but that nothing else corroborates. Interpret RELATIVE to another arm, never",
              "against an absolute cutoff: genuine low-abundance peptides do appear once."]
    md = "\n".join(lines) + "\n"
    print("\n" + md)
    if a.out:
        Path(a.out).write_text(md)
        print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
