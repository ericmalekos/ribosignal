#!/usr/bin/env python3
"""Score predicted vs observed per-nt Ribo-seq profiles for the cross-species packs.

Uses the project's OWN pearson/spearman from scripts/train.py rather than a fresh
implementation, so the numbers are comparable to everything already in results.md.

WHAT IS MEASURED, per transcript, then summarised as a median across transcripts:

  profile_r        pearson on log1p(counts). log1p because the raw P-site profile is
                   count data spanning orders of magnitude and a handful of very tall
                   peaks otherwise dominate the correlation; log1p is what the project
                   uses throughout (39 call sites).
  profile_rho      spearman on the raw profile, rank-based so it is insensitive to that
                   same scaling question. Reported alongside r deliberately: agreement
                   between the two is evidence the number is not an artifact of the
                   transform.
  total_r          pearson of log1p(pred_total) against log1p(observed total) across
                   transcripts, where pred_total is the model's COUNT head (pred_flat
                   cannot supply this: it is normalised per transcript and sums to 1). This is a
                   much easier task than the per-nt profile and is reported separately
                   so the two are never conflated -- getting a transcript's overall
                   translation level right says nothing about getting the shape right.

WHY A MEDIAN AND NOT A MEAN: the per-transcript distribution is long-tailed, and short
or low-count transcripts produce unstable correlations. The median is what results.md
reports as `pearson_median`.

A NULL IS INCLUDED, and it matters. `shuffled_r` recomputes profile_r after shuffling
each observed profile's positions. A model that has learned nothing about WHERE
ribosomes sit, but has learned a transcript's overall level, can still score a
respectable-looking profile_r; the shuffled control shows how much of the score is
positional information rather than magnitude.

Usage: python3 scripts/xspecies/eval_xspecies_profiles.py [--min-nt 100] [--min-obs 50]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "scripts"))
from train import pearson, spearman  # noqa: E402

ARCHS = ["attn", "mamba4"]
RUN = "orf_v2_{a}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"


def score(npz: Path, min_nt: int, min_obs: int, rng: np.random.Generator) -> dict:
    d = np.load(npz, allow_pickle=True)
    lens = d["lengths"]
    pred, obs = d["pred_flat"], d["obs_flat"]
    off = np.concatenate([[0], np.cumsum(lens)])
    # pred_flat is a normalised distribution over positions: every transcript's slice sums
    # to exactly 1.0. The per-transcript MAGNITUDE lives in the separate count head,
    # pred_total. Correlating log1p(pred) against log1p(obs) directly therefore compares a
    # probability against a count and is dominated by the scale difference, and summing
    # pred_flat per transcript yields the constant 1.0 -- which is why the first version of
    # this script reported total_r ~= 0 for every species. Scale the profile to the observed
    # depth before the log, and take the total from pred_total.
    ptot_head = d["pred_total"].astype(np.float64)
    rs, rhos, shuf = [], [], []
    ptot, ototo = [], []
    n_used = 0
    for i, L in enumerate(lens):
        p = pred[off[i]:off[i + 1]].astype(np.float64)
        o = obs[off[i]:off[i + 1]].astype(np.float64)
        osum = o.sum()
        ptot.append(ptot_head[i]); ototo.append(osum)
        if L < min_nt or osum < min_obs:
            continue
        ps = p * osum                      # same expected depth as the observation
        r = pearson(np.log1p(ps), np.log1p(o))
        rho = spearman(p, o)
        if np.isfinite(r):
            rs.append(r); n_used += 1
        if np.isfinite(rho):
            rhos.append(rho)
        os_ = o.copy(); rng.shuffle(os_)
        sr = pearson(np.log1p(ps), np.log1p(os_))
        if np.isfinite(sr):
            shuf.append(sr)
    ptot = np.asarray(ptot); ototo = np.asarray(ototo)
    tot_r = pearson(np.log1p(ptot), np.log1p(ototo))
    return {"n_tx": int(len(lens)), "n_scored": n_used,
            "profile_r": float(np.median(rs)) if rs else float("nan"),
            "profile_rho": float(np.median(rhos)) if rhos else float("nan"),
            "shuffled_r": float(np.median(shuf)) if shuf else float("nan"),
            "total_r": tot_r}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-nt", type=int, default=100)
    ap.add_argument("--min-obs", type=int, default=50)
    ap.add_argument("--tsv", type=Path, default=NEW / "results" / "xspecies_profile_eval.tsv")
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    rows = []
    for a in ARCHS:
        base = NEW / "results" / "loto" / RUN.format(a=a)
        for dd in sorted(base.glob("dropin_xsp_*")):
            npz = dd / "pred_profiles.npz"
            if not npz.exists():
                continue
            ds = dd.name.replace("dropin_", "")
            print(f"  scoring {a:7s} {ds} ...", file=sys.stderr, flush=True)
            rec = {"arch": a, "pack": ds}
            rec.update(score(npz, args.min_nt, args.min_obs, rng))
            rows.append(rec)

    if not rows:
        print("no prediction npz found", file=sys.stderr); return 1
    rows.sort(key=lambda r: (r["pack"], r["arch"]))
    hdr = f"  {'pack':<16}{'arch':<9}{'scored':>8}{'profile_r':>11}{'profile_rho':>13}{'shuffled_r':>12}{'total_r':>9}"
    print("\n" + hdr); print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['pack']:<16}{r['arch']:<9}{r['n_scored']:>8,}"
              f"{r['profile_r']:>11.3f}{r['profile_rho']:>13.3f}"
              f"{r['shuffled_r']:>12.3f}{r['total_r']:>9.3f}")
    args.tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\n  -> {args.tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
