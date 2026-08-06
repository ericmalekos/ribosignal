#!/usr/bin/env python3
"""Compare raw-depth vs depth-normalized (global_mean) RNA-seq input, pooled over 5 folds.

Confirms the depth normalization (train.py --cov_norm global_mean) preserves the within-dataset
results: it removes absolute sequencing depth from the coverage input for transferability, and
should NOT change the pooled pc / lncRNA / uORF / dORF medians (coverage barely drives shape, and
the count head only sees a rescaled abundance). Raw = results/improve/<cfg>_f*, normalized =
results/covnorm/<cfg>_f*, for baseline and orf_v2_attn. Prints raw vs normalized per region.
"""
from __future__ import annotations

import csv
import glob
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
CONFIGS = ["baseline", "orf_v2_attn"]
# region -> (pertx column, biotype filter)
REGIONS = [("pc", "profile_pearson", "protein_coding"),
           ("lncRNA", "profile_pearson", "lncRNA"),
           ("uORF", "utr5_pearson", "protein_coding"),
           ("dORF", "utr3_pearson", "protein_coding")]


def pool(run_glob, col, bt):
    vals = []
    for fp in sorted(glob.glob(run_glob)):
        with open(fp) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["biotype"] != bt:
                    continue
                try:
                    v = float(r[col])
                except (ValueError, KeyError):
                    continue
                if v == v:
                    vals.append(v)
    return (float(np.median(vals)), len(vals)) if vals else (None, 0)


def main():
    hdr = ("config", "region", "raw", "norm", "diff", "n_raw", "n_norm")
    print(f"{hdr[0]:<13} {hdr[1]:<7} {hdr[2]:>7} {hdr[3]:>7} {hdr[4]:>7} {hdr[5]:>7} {hdr[6]:>7}")
    for cfg in CONFIGS:
        for name, col, bt in REGIONS:
            rm, rn = pool(str(NEW / "results" / "improve" / f"{cfg}_f*" / "pertx.tsv"), col, bt)
            nm, nn = pool(str(NEW / "results" / "covnorm" / f"{cfg}_f*" / "pertx.tsv"), col, bt)
            diff = f"{nm - rm:+.3f}" if isinstance(rm, float) and isinstance(nm, float) else "  .  "
            rs = f"{rm:.3f}" if isinstance(rm, float) else "  .  "
            ns = f"{nm:.3f}" if isinstance(nm, float) else "  .  "
            print(f"{cfg:<13} {name:<7} {rs:>7} {ns:>7} {diff:>7} {rn:>7} {nn:>7}")
    print("\n(diff ~ 0 => depth normalization preserves within-dataset results while making the "
          "input depth-invariant; see the depth x10 check in results.md Task 10)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
