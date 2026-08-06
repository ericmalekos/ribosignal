#!/usr/bin/env python3
"""Is the cross-tissue (LOTO) profile Pearson genuine, or inflated by spiky transcripts?

Same idea as robustness_concentration.py but for the held-out-tissue runs: profile Pearson is
inflated when one nucleotide holds most of a transcript's P-sites (matching a single spike gives
r ~ 0.99 for free). This restricts each LOTO run's per-transcript Pearson to transcripts whose
signal is genuinely distributed, and checks whether the headline (orf_v2_attn ~= within-tissue on
pc, baseline worse) survives.

Concentration = max single-nt fraction over the metric's own window, measured on the HELD-OUT
tissue's observed P-site target (data/packed_<Tissue>/target_counts.npy) -- config-independent.
Strata: distributed (<0.10), moderate (0.10-0.30), spiky (>=0.30).

Usage: robustness_concentration_loto.py [Tissue=Hepatocytes]
Writes results/loto/robustness_concentration_<Tissue>.{tsv,json}, prints per-region tables.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
REGIONS = {"pc": ("profile_pearson", "whole", "protein_coding"),
           "lncRNA": ("profile_pearson", "whole", "lncRNA"),
           "uorf": ("utr5_pearson", "u5", "protein_coding"),
           "dorf": ("utr3_pearson", "u3", "protein_coding")}
DIST, SPIKY = 0.10, 0.30


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
    tissue = sys.argv[1] if len(sys.argv) > 1 else "Hepatocytes"
    pack = NEW / "data" / ("packed" if tissue == "Fibroblast" else f"packed_{tissue}")
    offsets = np.load(pack / "offsets.npy")
    counts = np.load(pack / "target_counts.npy", mmap_mode="r")
    row = {tx: i for i, tx in enumerate((pack / "tx_order.txt").read_text().split())}
    cds = load_cds()
    conc_cache: dict[str, dict] = {}

    def concentration(tx):
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
            if 0 < u5 <= L and c[:u5].sum() > 0:
                out["u5"] = float(c[:u5].max() / c[:u5].sum())
            if 0 < u3 <= L and c[L - u3:].sum() > 0:
                out["u3"] = float(c[L - u3:].max() / c[L - u3:].sum())
        conc_cache[tx] = out
        return out

    runs = sorted(d for d in (NEW / "results" / "loto").glob(f"*_holdout_{tissue}")
                  if (d / "pertx.tsv").exists())
    if not runs:
        print(f"no LOTO pertx for holdout={tissue}", file=sys.stderr)
        return 1

    def label(name):
        be = "orthrus" if "_orthrus_" in name else "rinalmo"
        cfg = name.split("_orthrus_")[0] if be == "orthrus" else name.split("_holdout_")[0]
        return f"{be}/{cfg}"

    data = {}
    for d in runs:
        pairs = {reg: [] for reg in REGIONS}
        with (d / "pertx.tsv").open() as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                bt = rec["biotype"]
                conc = concentration(rec["tx_id"])
                for reg, (col, wkey, want_bt) in REGIONS.items():
                    if bt != want_bt or wkey not in conc:
                        continue
                    try:
                        v = float(rec[col])
                    except (ValueError, KeyError):
                        continue
                    if v == v:
                        pairs[reg].append((v, conc[wkey]))
        data[label(d.name)] = pairs

    def smed(pairs, lo, hi):
        vals = [p for p, m in pairs if lo <= m < hi]
        return (float(np.median(vals)), len(vals)) if vals else (None, 0)

    def fmt(x):
        return f"{x:.3f}" if isinstance(x, float) else "  .  "

    # composition (from the first run's pairs; concentration is config-independent)
    first = next(iter(data.values()))
    print(f"=== holdout={tissue}: profile concentration composition ===")
    for reg in REGIONS:
        p = first[reg]
        if p:
            d_ = sum(1 for _, m in p if m < DIST)
            s_ = sum(1 for _, m in p if m >= SPIKY)
            print(f"  {reg:<7} n={len(p):>6}  distributed(<{DIST})={d_/len(p):5.1%}  "
                  f"spiky(>={SPIKY})={s_/len(p):5.1%}  "
                  f"median max-nt={np.median([m for _, m in p]):.1%}")

    out_rows = []
    for reg in REGIONS:
        print(f"\n=== {reg}: median profile Pearson by concentration stratum ===")
        print(f"{'run':<22} {'ALL':>7} {'n':>6}  {'DISTRIB':>8} {'n':>6}  {'SPIKY':>7} {'n':>6}")
        for lab, pairs in data.items():
            allm, alln = smed(pairs[reg], 0.0, 1.01)
            dm, dn = smed(pairs[reg], 0.0, DIST)
            sm, sn = smed(pairs[reg], SPIKY, 1.01)
            print(f"{lab:<22} {fmt(allm):>7} {alln:>6}  {fmt(dm):>8} {dn:>6}  {fmt(sm):>7} {sn:>6}")
            out_rows.append({"run": lab, "region": reg, "all": allm, "all_n": alln,
                             "dist": dm, "dist_n": dn, "spiky": sm, "spiky_n": sn})

    out = NEW / "results" / "loto" / f"robustness_concentration_{tissue}"
    out.with_suffix(".json").write_text(json.dumps(out_rows, indent=2))
    with out.with_suffix(".tsv").open("w") as fh:
        fh.write("run\tregion\tall\tall_n\tdist\tdist_n\tspiky\tspiky_n\n")
        for r in out_rows:
            fh.write("\t".join(str(r[k]) for k in
                     ("run", "region", "all", "all_n", "dist", "dist_n", "spiky", "spiky_n"))
                     + "\n")
    print(f"\nwrote {out.with_suffix('.tsv')}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
