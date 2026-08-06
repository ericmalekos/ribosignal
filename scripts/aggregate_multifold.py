#!/usr/bin/env python3
"""Aggregate the multi-fold improvement matrix into per-fold and POOLED views.

Pools all 5 gene-disjoint folds (each scorable transcript held out once), so deltas rest on the
full pooled set, not one noisy fold. Four scored regions, each a profile Pearson:
  - protein_coding : whole transcript (CDS-dominated)
  - lncRNA         : whole transcript (the lncRNA-ORF target); pooled ~1,239
  - uorf (5'UTR)   : pc 5'UTR window, translated-uORF transcripts; pooled ~22,000. uORFs are the
                     most common alternative ORF and the whole-tx pc metric barely reflects them.
  - dorf (3'UTR)   : pc 3'UTR window (downstream ORFs)

Two views per config x region: per-fold (mean +/- std of the per-fold medians) and pooled (median
over the concatenated per-transcript values from pertx.tsv, which eval_extra writes with columns
profile_pearson / utr5_pearson / utr3_pearson). Whole-tx per-fold medians come from
test_metrics.json; UTR per-fold medians from extra_metrics.json (uorf_5utr / dorf_3utr blocks).
Writes results/improve/summary_multifold.{tsv,json}.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
IMPROVE = NEW / "results" / "improve"
CONFIGS = ["baseline", "orf", "orf_v2", "orf_ctg", "attn", "orf_attn",
           "orf_v2_attn", "orf_ctg_attn"]
FOLDS = [0, 1, 2, 3, 4]
# regions pooled from pertx.tsv columns
POOL = {"protein_coding": "profile_pearson", "lncRNA": "profile_pearson",
        "uorf": "utr5_pearson", "dorf": "utr3_pearson"}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def read_pertx(path: Path):
    """region -> list of per-tx values (NaN dropped). uorf/dorf are pc-only UTR-window Pearsons."""
    out = {k: [] for k in POOL}
    if not path.exists():
        return out
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            b = row["biotype"]
            if b in ("protein_coding", "lncRNA"):
                v = _f(row.get("profile_pearson"))
                if v == v:
                    out[b].append(v)
            if b == "protein_coding":
                for key, col in (("uorf", "utr5_pearson"), ("dorf", "utr3_pearson")):
                    v = _f(row.get(col))
                    if v == v:
                        out[key].append(v)
    return out


def perfold_median(tm, em, region):
    """Per-fold median for a region from test_metrics (whole-tx) or extra_metrics (UTR)."""
    if region == "protein_coding":
        return (tm.get("protein_coding") or {}).get("pearson_median")
    if region == "lncRNA":
        return (tm.get("lncRNA") or {}).get("pearson_median")
    block = (em.get("uorf_5utr") if region == "uorf" else em.get("dorf_3utr")) or {}
    return block.get("profile_pearson_median")


def main():
    rows = []
    missing = []
    for cfg in CONFIGS:
        per_fold = {k: [] for k in POOL}
        pooled = {k: [] for k in POOL}
        n_folds = 0
        for k in FOLDS:
            rd = IMPROVE / f"{cfg}_f{k}"
            tm_p = rd / "test_metrics.json"
            if not tm_p.exists():
                missing.append(f"{cfg}_f{k}")
                continue
            n_folds += 1
            tm = json.loads(tm_p.read_text())
            em_p = rd / "extra_metrics.json"
            em = json.loads(em_p.read_text()) if em_p.exists() else {}
            for region in POOL:
                v = perfold_median(tm, em, region)
                if isinstance(v, int | float):
                    per_fold[region].append(v)
            px = read_pertx(rd / "pertx.tsv")
            for region in POOL:
                pooled[region].extend(px[region])

        rec = {"config": cfg, "n_folds": n_folds}
        for region in POOL:
            tag = {"protein_coding": "pc", "lncRNA": "lnc",
                   "uorf": "uorf", "dorf": "dorf"}[region]
            pf = np.array(per_fold[region], dtype=float)
            pl = np.array(pooled[region], dtype=float)
            rec[f"{tag}_perfold_mean"] = float(pf.mean()) if pf.size else None
            rec[f"{tag}_perfold_std"] = float(pf.std(ddof=0)) if pf.size else None
            rec[f"{tag}_pooled_median"] = float(np.median(pl)) if pl.size else None
            rec[f"{tag}_pooled_n"] = int(pl.size)
        rows.append(rec)

    base = next((r for r in rows if r["config"] == "baseline"), None)

    def d(r, key):
        if base is None or r is base:
            return ""
        a, b = r.get(key), base.get(key)
        if isinstance(a, int | float) and isinstance(b, int | float):
            return f"{a - b:+.3f}"
        return ""

    def f3(x):
        return f"{x:.3f}" if isinstance(x, int | float) else "  .  "

    print("=== multi-fold matrix: POOLED median profile Pearson (delta vs baseline) ===")
    print(f"{'config':<13} {'pc':>7} {'lnc':>7} {'dlnc':>8} {'uorf':>7} {'duorf':>8} "
          f"{'dorf':>7} {'lnc_n':>6} {'uorf_n':>7} {'folds':>5}")
    for r in rows:
        print(f"{r['config']:<13} {f3(r['pc_pooled_median']):>7} {f3(r['lnc_pooled_median']):>7} "
              f"{d(r, 'lnc_pooled_median'):>8} {f3(r['uorf_pooled_median']):>7} "
              f"{d(r, 'uorf_pooled_median'):>8} {f3(r['dorf_pooled_median']):>7} "
              f"{r['lnc_pooled_n']:>6} {r['uorf_pooled_n']:>7} {r['n_folds']:>5}")

    if missing:
        print(f"\n[pending] {len(missing)} cells: {', '.join(missing)}", file=sys.stderr)

    cols = ["config", "n_folds",
            "pc_perfold_mean", "pc_perfold_std", "pc_pooled_median", "pc_pooled_n",
            "lnc_perfold_mean", "lnc_perfold_std", "lnc_pooled_median", "lnc_pooled_n",
            "uorf_perfold_mean", "uorf_perfold_std", "uorf_pooled_median", "uorf_pooled_n",
            "dorf_perfold_mean", "dorf_perfold_std", "dorf_pooled_median", "dorf_pooled_n"]
    (IMPROVE / "summary_multifold.json").write_text(json.dumps(rows, indent=2))
    with (IMPROVE / "summary_multifold.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"\nwrote {IMPROVE/'summary_multifold.tsv'} and .json", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
