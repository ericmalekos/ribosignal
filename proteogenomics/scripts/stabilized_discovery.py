#!/usr/bin/env python3
"""Stabilize the mokapot discovery counts by ENSEMBLING the rescored score across the N replicate runs,
then applying the (deterministic) class-specific FDR once. This removes the run-to-run bimodality at its
source: instead of one noisy mokapot roll, every PSM gets a consensus score.

Cross-run scale caveat: mokapot's semi-supervised discriminant is not calibrated to a common scale across
runs, so we ensemble in two scale-robust ways and report both:
  - avg_pep    : mean posterior error probability (calibrated 0-1); score = -mean(pep)   [PRIMARY]
  - avg_rank   : mean within-run percentile rank of the raw score; score = mean(percentile)
plus the literal avg_score (mean raw score) for reference. A PSM is keyed by (spectrum_id, peptidoform);
rank-1 per spectrum is re-derived from the ensemble score; then compare_rescored.class_fdr_novel counts
novel peptides at 1% class FDR against REV_nuORF| decoys only.

Writes results/mokapot_stabilized.json. cas12a env. Usage: stabilized_discovery.py
"""
from __future__ import annotations

import ast
import csv
import glob
import json
import statistics as st
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "proteogenomics/scripts"))
from compare_rescored import bare, class_fdr_novel  # noqa: E402

DATASETS = {"HBL1": "HBL1_pilot", "DoHH2": "DoHH2_pilot", "SUDHL4": "SUDHL4_pilot"}
DBS = ["model", "null"]


def parse_accs(s):
    try:
        return tuple(ast.literal_eval(s))
    except Exception:
        return tuple(a.strip(" '\"[]") for a in s.split(",") if a.strip(" '\"[]"))


def load_runs(pilot, db):
    """Per replicate run: {(spectrum_id, peptidoform): (score, pep, accs)} and the run's score list
    (for percentile ranking). Returns list-of-runs."""
    runs = []
    for sd in sorted(glob.glob(str(NEW / f"proteogenomics/data/{pilot}/rescore_seed*"))):
        files = glob.glob(f"{sd}/{db}/*.psms.tsv")
        if not files:
            continue
        rows = {}
        for f in files:
            for r in csv.DictReader(open(f), delimiter="\t"):
                try:
                    sc = float(r["score"]); pep = float(r["pep"])
                except (ValueError, KeyError):
                    continue
                rows[(r["spectrum_id"], r["peptidoform"])] = (sc, pep, parse_accs(r["protein_list"]))
        if rows:
            runs.append(rows)
    return runs


def percentile_map(scores):
    """value -> percentile in [0,1] (average-rank of equal values)."""
    order = sorted(scores)
    n = len(order)
    import bisect
    return lambda v: (bisect.bisect_left(order, v) + bisect.bisect_right(order, v)) / (2.0 * n)


def ensemble_count(runs, score_of, fdr=0.01):
    """score_of(psm_key) -> ensemble score (higher=better). Re-derive rank-1 per spectrum, class-FDR."""
    # collect all keys present in >=1 run; ensemble score per key
    best = {}   # spectrum_id -> (ensemble_score, accs)
    keys = set().union(*[set(r) for r in runs])
    for key in keys:
        sc = score_of(key)
        if sc is None:
            continue
        accs = next(r[key][2] for r in runs if key in r)
        spec, pf = key
        if spec not in best or sc > best[spec][0]:
            best[spec] = (sc, accs, pf)
    psms = [(bare(pf), sc, accs) for (sc, accs, pf) in best.values()]
    keep, t, d = class_fdr_novel(psms, fdr)
    return len(keep)


def main():
    out = {"datasets": {}}
    for name, pilot in DATASETS.items():
        ds = {}
        for db in DBS:
            cm = NEW / f"proteogenomics/data/{pilot}/db/{db}_class_map.tsv"
            n_orfs = (sum(1 for _ in open(cm)) - 1) if cm.exists() else 0
            runs = load_runs(pilot, db)
            if not runs:
                continue
            # per-run percentile maps for the rank ensemble
            pmaps = [percentile_map([v[0] for v in r.values()]) for r in runs]

            def avg_pep(key):
                vals = [r[key][1] for r in runs if key in r]
                return -st.mean(vals) if vals else None

            def avg_rank(key):
                vals = [pm(r[key][0]) for r, pm in zip(runs, pmaps) if key in r]
                return st.mean(vals) if vals else None

            def avg_score(key):
                vals = [r[key][0] for r in runs if key in r]
                return st.mean(vals) if vals else None

            ds[db] = {
                "n_runs": len(runs), "n_orfs": n_orfs,
                "count_avg_pep": ensemble_count(runs, avg_pep),
                "count_avg_rank": ensemble_count(runs, avg_rank),
                "count_avg_score": ensemble_count(runs, avg_score)}
        out["datasets"][name] = ds
        for db in DBS:
            if db in ds:
                e = ds[db]
                r = 1000.0 / (e["n_orfs"] or 1)
                print(f"  {name} {db}: pep={e['count_avg_pep']} rank={e['count_avg_rank']} "
                      f"score={e['count_avg_score']}  (rate/1k pep {e['count_avg_pep']*r:.2f})")
    dst = NEW / "results/mokapot_stabilized.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
