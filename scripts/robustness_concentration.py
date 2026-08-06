#!/usr/bin/env python3
"""Robustness of the config rankings to profile concentration (spikiness).

Profile Pearson is inflated by spiky transcripts: if one nucleotide holds most of the P-sites
(e.g. a strong initiation pause), matching that single spike gives r ~ 0.99 for free. This asks
whether the improvement-matrix conclusions (non-AUG > ATG on uORF, attention helps, orf_v2_attn
best) survive when restricted to transcripts whose signal is genuinely distributed.

Concentration is measured on the SAME window as each metric: whole transcript for pc/lncRNA, the
5'UTR window for uORF, the 3'UTR window for dORF. The metric is max single-nt fraction
max(c)/sum(c) over the window (config-independent, computed once per transcript from the packed
target counts). Each config's per-transcript Pearson (from pertx.tsv, pooled over 5 folds) is then
bucketed by that concentration and the pooled median recomputed per stratum.

Strata: distributed (<0.10), moderate (0.10-0.30), spiky (>=0.30). Writes
results/improve/robustness_concentration.{tsv,json} and prints per-region tables.
"""
from __future__ import annotations

import csv
import glob
import json
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
PACK = NEW / "data" / "packed"
CONFIGS = ["baseline", "orf", "orf_v2", "orf_ctg", "attn", "orf_attn",
           "orf_v2_attn", "orf_ctg_attn"]
# region -> (pertx pearson column, window key for concentration)
REGIONS = {"pc": ("profile_pearson", "whole"), "lncRNA": ("profile_pearson", "whole"),
           "uorf": ("utr5_pearson", "u5"), "dorf": ("utr3_pearson", "u3")}
DIST, SPIKY = 0.10, 0.30


def load_counts_index():
    offsets = np.load(PACK / "offsets.npy")
    counts = np.load(PACK / "target_counts.npy", mmap_mode="r")
    order = (PACK / "tx_order.txt").read_text().split()
    row = {tx: i for i, tx in enumerate(order)}
    return offsets, counts, row


def load_cds():
    out = {}
    with (NEW / "data" / "tx2cds.tsv").open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len", "utr3_len")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[ix["tx_id"]]] = (int(f[ix["utr5_len"]]), int(f[ix["utr3_len"]]))
    return out


def main():
    offsets, counts, row = load_counts_index()
    cds = load_cds()
    conc_cache: dict[str, dict] = {}

    def concentration(tx):
        """window -> max single-nt fraction (max/sum), for whole / u5 / u3 windows."""
        if tx in conc_cache:
            return conc_cache[tx]
        r = row.get(tx)
        if r is None:
            conc_cache[tx] = {}
            return {}
        c = np.asarray(counts[offsets[r]:offsets[r + 1]], dtype=np.float64)
        L = c.shape[0]
        out = {}
        tot = c.sum()
        if tot > 0:
            out["whole"] = float(c.max() / tot)
        cc = cds.get(tx)
        if cc:
            u5, u3 = cc
            if 0 < u5 <= L:
                w = c[:u5]
                if w.sum() > 0:
                    out["u5"] = float(w.max() / w.sum())
            if 0 < u3 <= L:
                w = c[L - u3:]
                if w.sum() > 0:
                    out["u3"] = float(w.max() / w.sum())
        conc_cache[tx] = out
        return out

    # collect per (config, region): list of (pearson, concentration)
    data = {c: {reg: [] for reg in REGIONS} for c in CONFIGS}
    for cfg in CONFIGS:
        files = sorted(glob.glob(str(NEW / "results" / "improve" / f"{cfg}_f*" / "pertx.tsv")))
        for fp in files:
            with open(fp) as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    bt = rec["biotype"]
                    conc = concentration(rec["tx_id"])
                    for reg, (col, wkey) in REGIONS.items():
                        if reg in ("pc", "lncRNA") and bt != ("protein_coding" if reg == "pc"
                                                              else "lncRNA"):
                            continue
                        if reg in ("uorf", "dorf") and bt != "protein_coding":
                            continue
                        try:
                            v = float(rec[col])
                        except (ValueError, KeyError):
                            continue
                        if v != v or wkey not in conc:
                            continue
                        data[cfg][reg].append((v, conc[wkey]))

    def stratum_median(pairs, lo, hi):
        vals = [p for p, m in pairs if lo <= m < hi]
        return (float(np.median(vals)), len(vals)) if vals else (None, 0)

    out_rows = []
    for cfg in CONFIGS:
        rec = {"config": cfg}
        for reg in REGIONS:
            pairs = data[cfg][reg]
            allm, alln = stratum_median(pairs, 0.0, 1.01)
            dm, dn = stratum_median(pairs, 0.0, DIST)
            mm, mn = stratum_median(pairs, DIST, SPIKY)
            sm, sn = stratum_median(pairs, SPIKY, 1.01)
            rec[reg] = {"all": allm, "all_n": alln, "dist": dm, "dist_n": dn,
                        "mod": mm, "mod_n": mn, "spiky": sm, "spiky_n": sn}
        out_rows.append(rec)

    base = next(r for r in out_rows if r["config"] == "baseline")

    def fmt(x):
        return f"{x:.3f}" if isinstance(x, float) else "  .  "

    def dvs(r, reg, key):
        a, b = r[reg][key], base[reg][key]
        return f"{a - b:+.3f}" if isinstance(a, float) and isinstance(b, float) else "  .  "

    # composition: what fraction of scored tx are spiky, per region (from baseline's pairs)
    print("=== profile concentration composition (max single-nt fraction of window) ===")
    for reg in REGIONS:
        p = data["baseline"][reg]
        n = len(p)
        if not n:
            continue
        d = sum(1 for _, m in p if m < DIST)
        s = sum(1 for _, m in p if m >= SPIKY)
        print(f"  {reg:<7} n={n:>6}  distributed(<{DIST})={d/n:5.1%}  "
              f"spiky(>={SPIKY})={s/n:5.1%}  median max-nt={np.median([m for _, m in p]):.1%}")

    for reg in REGIONS:
        print(f"\n=== {reg}: pooled median profile Pearson by concentration stratum "
              f"(delta vs baseline) ===")
        print(f"{'config':<13} {'ALL':>7} {'n':>6}  {'DISTRIB':>8} {'dvs':>7} {'n':>6}  "
              f"{'SPIKY':>7} {'n':>6}")
        for r in out_rows:
            b = r[reg]
            print(f"{r['config']:<13} {fmt(b['all']):>7} {b['all_n']:>6}  "
                  f"{fmt(b['dist']):>8} {dvs(r, reg, 'dist'):>7} {b['dist_n']:>6}  "
                  f"{fmt(b['spiky']):>7} {b['spiky_n']:>6}")

    (NEW / "results" / "improve" / "robustness_concentration.json").write_text(
        json.dumps(out_rows, indent=2))
    cols = ["config"] + [f"{reg}_{k}" for reg in REGIONS
                         for k in ("all", "all_n", "dist", "dist_n", "spiky", "spiky_n")]
    with (NEW / "results" / "improve" / "robustness_concentration.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(
                [r["config"]] + [str(r[reg][k]) for reg in REGIONS
                                 for k in ("all", "all_n", "dist", "dist_n", "spiky", "spiky_n")]
            ) + "\n")
    print(f"\nwrote {NEW/'results'/'improve'/'robustness_concentration.tsv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
