#!/usr/bin/env python3
"""Aggregate the fold-0 improvement matrix into one table.

Compares, on the identical fold-0 split, the RiNALMo baseline against:
  - orf         sequence-derived ORF-candidate track, ATG-only start channel   (v1)
  - orf_v2      same track but graded start_ext over AUG + 9 near-cognate starts (non-AUG)
  - attn        + 2 self-attention layers (full-transcript context)
  - orf_attn    ORF (ATG) + attention
  - orf_v2_attn ORF (non-AUG) + attention

Answers: does the ORF track help, does attention help, and does non-AUG start information
(orf_v2 vs orf) matter. Reports protein_coding AND lncRNA (the goal is non-canonical ORFs, so
lncRNA is not secondary). Reads results/improve/<cfg>_f0 + the baseline
results/backend_cmp/rinalmo_f0; writes results/improve/summary.tsv + summary.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")

# (label, run dir, ORF track, attention layers)
ROWS = [
    ("baseline",    NEW / "results" / "backend_cmp" / "rinalmo_f0", "none",   0),
    ("orf",         NEW / "results" / "improve" / "orf_f0",         "ATG",    0),
    ("orf_v2",      NEW / "results" / "improve" / "orf_v2_f0",      "nonAUG", 0),
    ("attn",        NEW / "results" / "improve" / "attn_f0",        "none",   2),
    ("orf_attn",    NEW / "results" / "improve" / "orf_attn_f0",    "ATG",    2),
    ("orf_v2_attn", NEW / "results" / "improve" / "orf_v2_attn_f0", "nonAUG", 2),
]


def g(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def fmt(x):
    return f"{x:.3f}" if isinstance(x, int | float) else "  .  "


def main():
    out = []
    for label, run, orf, attn in ROWS:
        tm = run / "test_metrics.json"
        if not tm.exists():
            print(f"[pending] {label}: no test_metrics.json ({run})", file=sys.stderr)
            continue
        t = json.loads(tm.read_text())
        em = run / "extra_metrics.json"
        e = json.loads(em.read_text()) if em.exists() else {}
        out.append({
            "config": label, "orf_track": orf, "attn": attn,
            "pc_pearson": g(t, "protein_coding", "pearson_median"),
            "pc_period_pred": g(t, "protein_coding", "period_pred_median"),
            "pc_n": g(t, "protein_coding", "n"),
            "lnc_pearson": g(t, "lncRNA", "pearson_median"),
            "lnc_n": g(t, "lncRNA", "n"),
            "count_pearson": g(e, "protein_coding", "count_pearson"),
        })

    if not out:
        print("no completed improve runs yet", file=sys.stderr)
        return 1

    hdr = ["config", "orf", "attn", "pc_P", "pc_perP", "count_P", "lnc_P", "pc_n", "lnc_n"]
    w = [12, 7, 5, 7, 8, 8, 7, 6, 6]

    def line(vals):
        return "  ".join(str(v).ljust(wi) for v, wi in zip(vals, w, strict=False))

    print(line(hdr))
    base = next((r for r in out if r["config"] == "baseline"), None)
    for r in out:
        d = ""
        if base and r is not base and isinstance(r["pc_pearson"], int | float):
            d = f"  (pc {r['pc_pearson'] - base['pc_pearson']:+.3f} vs base)"
        print(line([r["config"], r["orf_track"], r["attn"], fmt(r["pc_pearson"]),
                    fmt(r["pc_period_pred"]), fmt(r["count_pearson"]), fmt(r["lnc_pearson"]),
                    r["pc_n"], r["lnc_n"]]) + d)

    cols = ["config", "orf_track", "attn", "pc_pearson", "pc_period_pred",
            "count_pearson", "lnc_pearson", "pc_n", "lnc_n"]
    outdir = NEW / "results" / "improve"
    (outdir / "summary.json").write_text(json.dumps(out, indent=2))
    with (outdir / "summary.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {outdir/'summary.tsv'} and summary.json", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
