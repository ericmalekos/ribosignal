#!/usr/bin/env python3
"""Score predicted vs observed per-nt Ribo-seq profiles for the cross-species packs.

Uses the project's OWN pearson/spearman definitions (copied verbatim from scripts/train.py).

CAUTION, corrected 2026-08-30: `profile_r` here is NOT the same statistic as the
`pearson_median` reported by train.py/eval_extra.py and tabulated throughout results.md.
Those compute pearson(p, c/N) on the RAW profile; this applies log1p first. The two are
not interchangeable: measured per transcript, they rank transcripts at Spearman 0.65-0.79
(human 0.734, worm 0.651) and share only 34-42% of their top decile. Do not compare
`profile_r` against an existing `pearson_median` figure.

WHAT IS MEASURED, per transcript, then summarised as a median across transcripts:

  profile_r        pearson on log1p(counts). log1p because the raw P-site profile is
                   count data spanning orders of magnitude and a handful of very tall
                   peaks otherwise dominate the correlation. NOTE: this transform is
                   LOCAL to this script. An earlier version of this docstring claimed
                   log1p was the project convention "throughout (39 call sites)"; that was
                   wrong. The project's other log1p uses are on transcript TOTALS
                   (abundance, count head) and on the RNA coverage INPUT, never on a
                   within-transcript profile before correlating.
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
import math
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "scripts"))
from scipy.stats import rankdata  # noqa: E402

# pearson and _rank/spearman are COPIED VERBATIM from scripts/train.py rather than
# imported, because importing train pulls in torch, which the riboseq env does not have.
# Copying is justified here only because both are three lines and are reproduced exactly;
# if train.py's definitions change, these must be updated with them.


def pearson(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b) / d if d > 0 else float("nan")


def _rank(x):
    r = np.empty_like(x)
    r[np.argsort(x, kind="stable")] = np.arange(len(x))
    return r


def spearman(a, b):
    return pearson(_rank(a).astype(np.float64), _rank(b).astype(np.float64))


def spearman_ties(a, b):
    """Spearman with AVERAGE ranks for ties, which is the standard definition.

    The project's train.spearman ranks via `_rank`, which uses
    `argsort(kind="stable")` and assigns ORDINAL ranks 0..n-1. That breaks ties by
    POSITION rather than averaging them. On a dense vector the two agree; on a sparse
    P-site profile they do not, and the difference is not small.

    A per-nt Ribo-seq profile is 74-94% zeros here. Under ordinal ranking every tied
    zero receives a distinct rank in positional order, so the observed rank vector
    becomes a proxy for POSITION along the transcript rather than for signal. The model
    predicts low density in the 3'UTR; in a species with long 3'UTRs those positions are
    both observed-zero (hence high ordinal rank, being last) and predicted-low, which
    manufactures an ANTI-correlation out of nothing. Measured, the effect tracks the zero
    fraction exactly:

        pack        zero%   ordinal   avg-ties
        yeast       78.0%     0.163      0.287
        celegans    73.8%     0.142      0.321
        gorilla     86.6%    -0.075      0.394
        human       86.9%    -0.190      0.394
        zebrafish   94.4%    -0.362      0.148

    So the negative `profile_rho` values this script reported before 2026-08-29 were a
    tie-handling artifact, NOT evidence of anti-correlated profiles. Both are emitted now:
    `profile_rho` (tie-corrected, the one to use) and `profile_rho_ordinal` (the old
    definition, kept so the change is auditable).

    NOTE, unfixed deliberately: train._rank is shared with eval_localization.py and
    inspect_distributions.py and feeds `spearman_median` columns in existing results
    tables. Correcting it there would silently move previously published numbers, so it
    is flagged rather than changed. Any Spearman computed on sparse data through
    train.spearman is affected.
    """
    return pearson(rankdata(a).astype(np.float64), rankdata(b).astype(np.float64))

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
    rs, rhos, rhos_ord, shuf = [], [], [], []
    ptot, ototo = [], []
    n_used = 0
    zfr = []
    for i, L in enumerate(lens):
        p = pred[off[i]:off[i + 1]].astype(np.float64)
        o = obs[off[i]:off[i + 1]].astype(np.float64)
        osum = o.sum()
        if L >= min_nt and osum >= min_obs:
            zfr.append(float((o == 0).mean()))
        ptot.append(ptot_head[i]); ototo.append(osum)
        if L < min_nt or osum < min_obs:
            continue
        ps = p * osum                      # same expected depth as the observation
        r = pearson(np.log1p(ps), np.log1p(o))
        rho = spearman_ties(p, o)
        rho_ord = spearman(p, o)
        if np.isfinite(r):
            rs.append(r); n_used += 1
        if np.isfinite(rho):
            rhos.append(rho)
        if np.isfinite(rho_ord):
            rhos_ord.append(rho_ord)
        os_ = o.copy(); rng.shuffle(os_)
        sr = pearson(np.log1p(ps), np.log1p(os_))
        if np.isfinite(sr):
            shuf.append(sr)
    ptot = np.asarray(ptot); ototo = np.asarray(ototo)
    tot_r = pearson(np.log1p(ptot), np.log1p(ototo))
    return {"n_tx": int(len(lens)), "n_scored": n_used,
            "profile_r": float(np.median(rs)) if rs else float("nan"),
            "profile_rho": float(np.median(rhos)) if rhos else float("nan"),
            "profile_rho_ordinal": float(np.median(rhos_ord)) if rhos_ord else float("nan"),
            "zero_frac": float(np.mean(zfr)) if zfr else float("nan"),
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
