#!/usr/bin/env python3
"""Consolidated cross-tissue (LOTO) table from each run's pertx.tsv.

For every results/loto/*_holdout_<Tissue>[_orthrus] run, computes the median per-transcript
profile Pearson for four regions -- pc (whole-tx, protein_coding), lncRNA (whole-tx),
uORF (5'UTR window), dORF (3'UTR window) -- directly from pertx.tsv, so all backends/configs
are scored identically and the Orthrus arm is picked up automatically when it lands.

Compared against the WITHIN-Fibroblast ceiling at the SAME coverage normalization
(cov_norm=global_mean), i.e. the Task 10 depth-normalized numbers -- the correct apples-to-
apples reference for the normalized LOTO models (not the raw Task 9 numbers).

Writes results/loto/loto_summary.tsv and prints the table.
"""
import csv
import statistics
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
LOTO = NEW / "results" / "loto"

# within-Fibroblast, cov_norm=global_mean (results.md Task 10 covnorm table) -- the matched ceiling
CEILING = {
    "baseline":    {"pc": 0.593, "lncRNA": 0.363, "uORF": 0.544, "dORF": 0.285},
    "orf_v2_attn": {"pc": 0.609, "lncRNA": 0.387, "uORF": 0.595, "dORF": 0.340},
}


def med(vals):
    vals = [v for v in vals if v == v]  # drop nan
    return statistics.median(vals) if vals else float("nan")


def _get(r, k):
    try:
        return float(r[k])
    except (ValueError, KeyError):
        return float("nan")


def region_medians(pertx):
    pc, lnc, u5, u3 = [], [], [], []
    with open(pertx) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["biotype"] == "protein_coding":
                pc.append(_get(r, "profile_pearson"))
                u5.append(_get(r, "utr5_pearson"))
                u3.append(_get(r, "utr3_pearson"))
            elif r["biotype"] == "lncRNA":
                lnc.append(_get(r, "profile_pearson"))
    return {"pc": med(pc), "lncRNA": med(lnc), "uORF": med(u5), "dORF": med(u3),
            "n_pc": len([v for v in pc if v == v]),
            "n_lnc": len([v for v in lnc if v == v])}


def parse_run(name):
    # Backend is encoded as a "_<backend>_" infix before "_holdout_"; rinalmo is the un-suffixed
    # default (e.g. orf_v2_attn_holdout_X = rinalmo, orf_v2_attn_onehot_holdout_X = onehot).
    for b in ("orthrus", "onehot", "concat"):
        if f"_{b}_" in name:
            return b, name.split(f"_{b}_")[0]
    return "rinalmo", name.split("_holdout_")[0]


def main():
    runs = sorted(d for d in LOTO.glob("*_holdout_*") if (d / "pertx.tsv").exists())
    rows = []
    for d in runs:
        backend, cfg = parse_run(d.name)
        m = region_medians(d / "pertx.tsv")
        rows.append((backend, cfg, m))

    hdr = ("backend", "config", "pc", "lncRNA", "uORF", "dORF", "n_pc", "n_lnc")
    print(f"{'backend':9} {'config':13} {'pc':>7} {'lncRNA':>7} {'uORF':>7} {'dORF':>7} "
          f"{'n_pc':>6} {'n_lnc':>6}")
    out_lines = ["\t".join(hdr)]
    for backend, cfg, m in rows:
        print(f"{backend:9} {cfg:13} {m['pc']:7.4f} {m['lncRNA']:7.4f} {m['uORF']:7.4f} "
              f"{m['dORF']:7.4f} {m['n_pc']:6d} {m['n_lnc']:6d}")
        out_lines.append(f"{backend}\t{cfg}\t{m['pc']:.4f}\t{m['lncRNA']:.4f}\t"
                         f"{m['uORF']:.4f}\t{m['dORF']:.4f}\t{m['n_pc']}\t{m['n_lnc']}")

    print("\n-- within-Fibroblast ceiling (cov_norm=global_mean, Task 10) --")
    for cfg, c in CEILING.items():
        print(f"{'within-Fib':9} {cfg:13} {c['pc']:7.4f} {c['lncRNA']:7.4f} {c['uORF']:7.4f} "
              f"{c['dORF']:7.4f}")
        out_lines.append(f"within-Fib\t{cfg}\t{c['pc']:.4f}\t{c['lncRNA']:.4f}\t"
                         f"{c['uORF']:.4f}\t{c['dORF']:.4f}\t.\t.")

    (LOTO / "loto_summary.tsv").write_text("\n".join(out_lines) + "\n")
    print(f"\nwrote {LOTO / 'loto_summary.tsv'}  ({len(rows)} run(s) scored)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
