#!/usr/bin/env python3
"""Aggregate the RiNALMo-vs-Orthrus-vs-concat fold-0 comparison into one table.

Reads results/backend_cmp/<backend>_f0/{test_metrics.json, extra_metrics.json, args.json}
for each backend and prints a per-backend row: profile Pearson (all / protein_coding /
lncRNA), observed vs predicted periodicity, count-head magnitude Pearson, model params,
and d_emb. Writes results/backend_cmp/summary.tsv and summary.json.

The three runs share the identical gene-disjoint fold-0 split (usable() intersects the
per-backend indices), so differences are attributable to the embedding, not the tx set.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
CMP = NEW / "results" / "backend_cmp"
BACKENDS = ["rinalmo", "orthrus", "concat"]
D_EMB = {"rinalmo": 1280, "orthrus": 512, "concat": 1792}


def g(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def fmt(x):
    return f"{x:.3f}" if isinstance(x, int | float) else "  .  "


def main():
    rows = []
    for b in BACKENDS:
        run = CMP / f"{b}_f0"
        tm = run / "test_metrics.json"
        em = run / "extra_metrics.json"
        if not tm.exists():
            print(f"[skip] {b}: no test_metrics.json yet ({run})", file=sys.stderr)
            continue
        t = json.loads(tm.read_text())
        e = json.loads(em.read_text()) if em.exists() else {}
        row = {
            "backend": b,
            "d_emb": D_EMB[b],
            "n_test": g(t, "n", default=0),
            "all_pearson": g(t, "all", "pearson_median"),
            "all_spearman": g(t, "all", "spearman_median"),
            "pc_pearson": g(t, "protein_coding", "pearson_median"),
            "pc_spearman": g(t, "protein_coding", "spearman_median"),
            "pc_period_pred": g(t, "protein_coding", "period_pred_median"),
            "pc_period_obs": g(t, "protein_coding", "period_obs_median"),
            "pc_n": g(t, "protein_coding", "n"),
            "lnc_pearson": g(t, "lncRNA", "pearson_median"),
            "lnc_n": g(t, "lncRNA", "n"),
            "count_pearson": g(e, "protein_coding", "count_pearson"),
        }
        rows.append(row)

    if not rows:
        print("no completed backend runs found", file=sys.stderr)
        return 1

    cols = ["backend", "d_emb", "pc_pearson", "pc_spearman", "pc_period_pred",
            "pc_period_obs", "count_pearson", "lnc_pearson", "pc_n", "lnc_n"]
    hdr = ["backend", "d_emb", "pc_prof_P", "pc_prof_S", "pc_per_pred",
           "pc_per_obs", "count_P", "lnc_prof_P", "pc_n", "lnc_n"]
    w = [10, 6, 10, 10, 11, 10, 8, 10, 6, 6]

    def line(vals):
        return "  ".join(str(v).ljust(wi) for v, wi in zip(vals, w, strict=False))

    print(line(hdr))
    for r in rows:
        vals = [r["backend"], r["d_emb"], fmt(r["pc_pearson"]), fmt(r["pc_spearman"]),
                fmt(r["pc_period_pred"]), fmt(r["pc_period_obs"]), fmt(r["count_pearson"]),
                fmt(r["lnc_pearson"]), r["pc_n"], r["lnc_n"]]
        print(line(vals))

    CMP.mkdir(parents=True, exist_ok=True)
    (CMP / "summary.json").write_text(json.dumps(rows, indent=2))
    with (CMP / "summary.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {CMP/'summary.tsv'} and summary.json", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
