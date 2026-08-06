#!/usr/bin/env python3
"""Consensus ORF calling over Poisson draws (`pgx`).

The calibrated arm samples `Poisson(rate)`, so a single run's call list is ONE SAMPLE from a
distribution, not a fixed object. Measured on BMDM at theta* = 0.05 across five seeds:

    total calls          10,823 +- 23          (CV 0.21%, counts are stable)
    pairwise Jaccard     0.822                 (membership is NOT)
    called in all 5      8,869  = 66.9% of the 13,250-call union
    called in exactly 1  1,673  = 12.6% of the union

So roughly 18% of any single draw's calls are draw-specific. The instability sits almost entirely in
the non-canonical classes -- `annotated` varies by 0.16%, while dORF varies ~17% -- which is
precisely where a proteogenomic search database gets its content. Requiring an ORF to recur across
independent draws removes that tail while keeping the reproducible core, exactly as fraction-level
replication does for peptides (pgx.replication).

Use this whenever the call list is an END PRODUCT (a search database, a published table) rather than
an intermediate. The deterministic arm (theta = 1, rounding) needs none of this: it is exactly
reproducible.

  python -m pgx.consensus --profiles p.npz --species mouse --theta 0.05 \
      --out <dir>/consensus --seeds 5 --min-seeds 3
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from . import rc_io
from .calibrate import collapsed_path, run_theta
from .call_orfs import COLS, write_calls
from .refs import species_refs


def orf_key(o):
    """Identity of a call, independent of which draw produced it."""
    return (o["tx"], o["start0"], o["end0"])


def gather(collapsed_paths, biotype, qvalue):
    """-> (counts per ORF key, best record per key, per-draw key sets)."""
    counts = Counter()
    best = {}
    per_draw = []
    for p in collapsed_paths:
        keys = set()
        for o in rc_io.read_collapsed(p):
            if o["qval"] > qvalue:
                continue
            k = orf_key(o)
            keys.add(k)
            prev = best.get(k)
            if prev is None or o["qval"] < prev["qval"]:
                best[k] = o
        per_draw.append(keys)
        counts.update(keys)
    return counts, best, per_draw


def jaccard(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--theta", type=float, required=True, help="calibrated theta (from theta.json)")
    ap.add_argument("--seeds", type=int, default=5, help="number of independent Poisson draws")
    ap.add_argument("--min-seeds", type=int, default=3, help="keep ORFs called in >= this many draws")
    ap.add_argument("--calib-dir", default=None,
                    help="calibration sweep dir; its seed-0 draw at --theta is reused if present")
    ap.add_argument("--orf-qvalue", type=float, default=0.05)
    ap.add_argument("--min-aa", type=int, default=5)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--annot", default=None)
    ap.add_argument("--tx2biotype", default=None)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    if not 1 <= a.min_seeds <= a.seeds:
        raise SystemExit(f"--min-seeds must be in [1, {a.seeds}]")
    refs = species_refs(a.species, ribocode_annot=a.annot, tx2biotype=a.tx2biotype)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    biotype = rc_io.load_biotype(refs["tx2biotype"])

    paths = []
    for s in range(a.seeds):
        cached = (collapsed_path(Path(a.calib_dir), a.theta)
                  if (a.calib_dir and s == 0) else None)
        if cached and cached.exists() and not a.force:
            print(f"[reuse] seed 0 <- {cached}", file=sys.stderr)
            paths.append(cached)
            continue
        d = out / f"seed_{s}"
        paths.append(run_theta(a.theta, a.profiles, refs["ribocode_annot"], d,
                               min_aa=a.min_aa, pval=a.pval, seed=s, poisson=True,
                               force=a.force))

    counts, best, per_draw = gather(paths, biotype, a.orf_qvalue)
    union = set(counts)
    core = {k for k, n in counts.items() if n == a.seeds}
    keep = {k for k, n in counts.items() if n >= a.min_seeds}
    once = sum(1 for n in counts.values() if n == 1)
    js = [jaccard(per_draw[i], per_draw[j])
          for i in range(len(per_draw)) for j in range(i + 1, len(per_draw))]

    rows = []
    for k in sorted(keep):
        o = best[k]
        if o["orf_type"] in rc_io.SKIP_TYPES:
            continue
        kl = rc_io.pgx_class(o, biotype.get(o["tx"]))
        if not kl or not o["aaseq"]:
            continue
        rows.append({**{c: o.get(c) for c in COLS if c in o}, "pgx_class": kl,
                     "n_seeds": counts[k]})
    write_calls(rows, out / "calls_consensus.tsv")

    stats = {
        "theta": a.theta, "seeds": a.seeds, "min_seeds": a.min_seeds,
        "orf_qvalue": a.orf_qvalue,
        "per_draw_calls": [len(s) for s in per_draw],
        "union": len(union), "core_all_seeds": len(core),
        "kept": len(keep), "called_once": once,
        "mean_pairwise_jaccard": round(sum(js) / len(js), 4) if js else None,
        "novel_rows": len(rows),
        "novel_by_class": dict(Counter(r["pgx_class"] for r in rows).most_common()),
        "collapsed": [str(p) for p in paths],
    }
    (out / "consensus_summary.json").write_text(json.dumps(stats, indent=2) + "\n")

    print(f"draws                : {[len(s) for s in per_draw]}")
    print(f"union                : {len(union):,}")
    print(f"core (all {a.seeds} draws)   : {len(core):,} ({len(core) / len(union) * 100:.1f}%)")
    print(f"called exactly once  : {once:,} ({once / len(union) * 100:.1f}%)")
    print(f"mean pairwise Jaccard: {stats['mean_pairwise_jaccard']}")
    print(f"kept (>= {a.min_seeds} draws)     : {len(keep):,} -> {len(rows):,} novel rows "
          f"{stats['novel_by_class']}")
    print(f"\nwrote {out}/calls_consensus.tsv + consensus_summary.json")


if __name__ == "__main__":
    main()
