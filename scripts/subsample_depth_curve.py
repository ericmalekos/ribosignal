#!/usr/bin/env python3
"""Depth-crossover analysis: does predicting the Ribo-seq profile from RNA-seq (no experiment) beat
MEASURING it at a given sequencing depth, for ORF calling?

Builds two things, both scored against the SAME deep ground truth (the official Hepatocytes RiboCode
calls, restricted to the model's held-out test-tx genes, genomic locus key):

  * MEASUREMENT CURVE  F1_measure(f): RiboCode calls on the REAL P-site profile binomially thinned to
    a fraction f of depth (subsample_depth.sbatch outputs real_depth<f>_seed<s>_collapsed.txt), vs the
    deep official calls. Decreasing as f -> 0. f=1 is the un-thinned real_collapsed.txt.
  * PREDICTION LINES (flat, depth-independent): the model's pred_preddepth (predicted shape x count-head
    predicted depth -- fully standalone, uses RNA-seq + sequence, ZERO Ribo-seq) and pred_obsdepth
    (predicted shape at full real depth), each vs the deep official calls.

The crossover f* where F1_measure(f) drops below the pred_preddepth line is the depth below which
predicting from RNA-seq recovers more of the deep truth than measuring at that depth. Reported as a
fraction and as an absolute test-set P-site count (f * full-depth).

Reuses compare_dropin_calls.{load_tx2gene,load_calls,compare}. cas12a env (numpy only).
Usage: subsample_depth_curve.py --dropin_dir <dir> --official <Hepatocytes_collapsed.txt> \
    --profiles <pred_profiles.npz> [--subsample_dir <dir/subsample>] [--out <dir>]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_dropin_calls import TX2GENE, compare, load_calls, load_tx2gene  # noqa: E402

DEPTH_RE = re.compile(r"real_depth([0-9.eE+-]+)_seed(\d+)_collapsed\.txt$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dropin_dir", required=True, help="the run's dropin/ dir (real + pred collapsed)")
    ap.add_argument("--official", required=True)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--subsample_dir", default=None, help="default <dropin_dir>/subsample")
    ap.add_argument("--tx2gene", default=TX2GENE)
    ap.add_argument("--out", default=None, help="default <dropin_dir>")
    ap.add_argument("--min_len", type=int, default=90)
    ap.add_argument("--max_pval", type=float, default=0.05)
    ap.add_argument("--min_enrichment", type=float, default=0.5)
    args = ap.parse_args()

    d = Path(args.dropin_dir)
    sub = Path(args.subsample_dir) if args.subsample_dir else d / "subsample"
    out = Path(args.out) if args.out else d
    out.mkdir(parents=True, exist_ok=True)

    npz = np.load(args.profiles, allow_pickle=False)
    keep_tx = set(npz["tx_ids"].tolist())
    lengths = npz["lengths"]
    pflat = npz["pred_flat"]
    obs_flat = npz["obs_flat"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    idx = {str(t): i for i, t in enumerate(npz["tx_ids"])}
    full_depth = float(obs_flat.sum())

    def enrich(tx, a, e):
        i = idx.get(str(tx))
        if i is None:
            return float("nan")
        p = pflat[off[i]:off[i + 1]]
        if a < 0 or e > len(p) or e <= a:
            return float("nan")
        return float(p[a:e].mean() * len(p))

    t2g = load_tx2gene(args.tx2gene)
    keep_genes = {t2g[t] for t in keep_tx if t in t2g}
    if keep_tx and not keep_genes:
        raise SystemExit(f"ERROR: 0/{len(keep_tx)} test tx matched {args.tx2gene}; wrong assembly?")
    print(f"test tx {len(keep_tx):,}; test genes {len(keep_genes):,}; full-depth P-sites {full_depth:.3e}",
          file=sys.stderr)

    def lc(path, is_pred=False):
        return load_calls(path, "genomic", keep_tx, keep_genes, args.max_pval, args.min_len,
                          enrich if is_pred else None, args.min_enrichment)

    official = lc(d / "official_collapsed.txt") if (d / "official_collapsed.txt").exists() \
        else lc(Path(args.official))

    def scored(path, is_pred=False):
        c = compare(lc(path, is_pred), official, path.name)
        nc = c["recall_grouped"]["noncanonical"]["recall"]
        return {"n_pred": c["n_pred"], "n_real": c["n_real"], "n_match": c["n_match"],
                "precision": c["precision"], "recall": c["recall"], "f1": c["f1"], "noncanon_recall": nc}

    rows = []
    # prediction lines (flat)
    for kind, fname, isp in (("pred_preddepth", "pred_preddepth_collapsed.txt", True),
                             ("pred_obsdepth", "pred_obsdepth_collapsed.txt", True)):
        p = d / fname
        if p.exists():
            rows.append({"kind": kind, "depth_frac": None, "abs_psites": None, "seed": None,
                         **scored(p, isp)})
    # measurement curve: f=1 (un-thinned) + thinned depths
    if (d / "real_collapsed.txt").exists():
        rows.append({"kind": "measure", "depth_frac": 1.0, "abs_psites": full_depth, "seed": None,
                     **scored(d / "real_collapsed.txt")})
    for p in sorted(sub.glob("real_depth*_seed*_collapsed.txt")):
        m = DEPTH_RE.search(p.name)
        if not m:
            continue
        f = float(m.group(1)); s = int(m.group(2))
        rows.append({"kind": "measure", "depth_frac": f, "abs_psites": f * full_depth, "seed": s,
                     **scored(p)})

    # write per-run TSV
    cols = ["kind", "depth_frac", "abs_psites", "seed", "n_pred", "n_real", "n_match",
            "precision", "recall", "f1", "noncanon_recall"]
    tsv = out / "depth_curve.tsv"
    with tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else (f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]))
                               for c in cols) + "\n")

    # crossover: mean measurement F1 per depth vs the pred_preddepth line, on log10(abs_psites)
    pred_line = next((r["f1"] for r in rows if r["kind"] == "pred_preddepth"), None)
    meas = {}
    for r in rows:
        if r["kind"] == "measure":
            meas.setdefault(r["depth_frac"], []).append(r["f1"])
    curve = sorted(((f, float(np.mean(v))) for f, v in meas.items()), key=lambda x: x[0])
    crossover = None
    if pred_line is not None and len(curve) >= 2:
        xs = [np.log10(f * full_depth) for f, _ in curve]
        ys = [y for _, y in curve]
        for i in range(len(xs) - 1):
            lo, hi = ys[i], ys[i + 1]
            if (lo - pred_line) * (hi - pred_line) <= 0 and hi != lo:
                frac = (pred_line - lo) / (hi - lo)
                xstar = xs[i] + frac * (xs[i + 1] - xs[i])
                depth_star = 10 ** xstar
                crossover = {"abs_psites": depth_star, "depth_frac": depth_star / full_depth,
                             "pred_preddepth_f1": pred_line}
                break

    summary = {"full_depth_psites": full_depth, "pred_preddepth_f1": pred_line,
               "pred_obsdepth_f1": next((r["f1"] for r in rows if r["kind"] == "pred_obsdepth"), None),
               "crossover": crossover,
               "measurement_curve": [{"depth_frac": f, "abs_psites": f * full_depth, "mean_f1": y}
                                     for f, y in curve]}
    (out / "depth_curve_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\nfull-depth test-set P-sites: {full_depth:.3e}", file=sys.stderr)
    print(f"pred_preddepth F1 (standalone, no Ribo-seq) = {pred_line}", file=sys.stderr)
    print(f"pred_obsdepth  F1 (pred shape, full depth)  = {summary['pred_obsdepth_f1']}", file=sys.stderr)
    print("\nmeasurement curve (measure at depth f, vs deep official):", file=sys.stderr)
    print(f"{'depth_frac':>11} {'abs_psites':>12} {'mean_F1':>8}", file=sys.stderr)
    for f, y in sorted(curve, reverse=True):
        print(f"{f:>11.4g} {f*full_depth:>12.3e} {y:>8.3f}", file=sys.stderr)
    if crossover:
        print(f"\nCROSSOVER: prediction beats measurement below ~{crossover['abs_psites']:.2e} "
              f"P-sites (f~{crossover['depth_frac']:.4g}), where measurement F1 falls to the "
              f"pred_preddepth line {pred_line:.3f}", file=sys.stderr)
    else:
        print("\nNO CROSSOVER in the sampled range (measurement stays above the prediction line, "
              "or the grid does not bracket it -- extend DEPTHS downward).", file=sys.stderr)
    print(f"\nwrote {tsv} and {out/'depth_curve_summary.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
