#!/usr/bin/env python3
"""Validate the pooled Fibroblast RNAseq coverage against the P-site target and salmon.

Fast checks, mostly from the summary TSVs (no hd5 row reads except two spot controls):
  1. Axis: pooled coverage tx set == target tx set (order may differ; join by id).
  2. Controls: ACTB/GAPDH pooled coverage is large and its per-nt length equals the
     target per-nt length (join by id) -- proves the per-nt axis is the transcript axis.
  3. Coverage vs abundance over the expressed universe:
       - mean per-nt coverage (total/length) vs salmon mean TPM (length-normalized, the
         clean abundance check), log1p Pearson -- expect strongly positive.
       - total coverage vs total P-sites, log1p Pearson -- expect ~0.7-0.8 (both are
         per-nt totals, so the length confound largely cancels).
Writes logs/validate_pooled_coverage.txt.
"""
import math
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
COV_HD5 = NEW / "data" / "rnaseq_coverage" / "Fibroblast_rnaseq_coverage_pooled.hd5"
COV_SUM = NEW / "data" / "rnaseq_coverage" / "Fibroblast_rnaseq_coverage_summary.tsv"
TGT_HD5 = NEW / "data" / "target" / "Fibroblast_psites_pooled.hd5"
TGT_SUM = NEW / "data" / "target" / "Fibroblast_psites_summary.tsv"
UNIV = NEW / "data" / "fibroblast_universe.tsv"
SALMON = NEW / "data" / "fibroblast_salmon_mean_tpm.tsv"
OUT = NEW / "logs" / "validate_pooled_coverage.txt"


def col_map(path):
    with path.open() as fh:
        return {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}


def load_col(path, key_col, val_col, cast=float):
    c = col_map(path)
    d = {}
    with path.open() as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[c[key_col]]] = cast(f[c[val_col]])
    return d


def logpearson(xs, ys):
    x = np.log1p(np.asarray(xs, dtype=np.float64))
    y = np.log1p(np.asarray(ys, dtype=np.float64))
    x -= x.mean()
    y -= y.mean()
    d = math.sqrt(float(x @ x) * float(y @ y))
    return float(x @ y) / d if d > 0 else float("nan")


def main():
    R = []

    def p(s=""):
        R.append(s)
        print(s, file=sys.stderr)

    universe = set(load_col(UNIV, "tx_id", "length", int))
    p(f"universe transcripts: {len(universe):,}")

    # [2] controls from hd5 (axis + length join)
    p("\n[controls] pooled coverage vs target per-nt length (join by id)")
    with h5py.File(COV_HD5, "r") as hc, h5py.File(TGT_HD5, "r") as ht:
        cid = hc["transcript_ids"].asstr()[:]
        tid = ht["transcript_ids"].asstr()[:]
        p(f"  coverage tx: {len(cid):,}  target tx: {len(tid):,}  "
          f"same set: {set(cid) == set(tid)}")
        cix = {t: i for i, t in enumerate(cid)}
        tix = {t: i for i, t in enumerate(tid)}
        for name, tx in [("ACTB", "ENST00000674681.1"), ("GAPDH", "ENST00000229239.10")]:
            if tx in cix and tx in tix:
                cc = np.asarray(hc["coverage"][cix[tx]], dtype=np.int64)
                tt = np.asarray(ht["p_sites"][tix[tx]], dtype=np.int64)
                p(f"  {name} {tx}: cov_len={cc.shape[0]} tgt_len={tt.shape[0]} "
                  f"len_match={cc.shape[0] == tt.shape[0]}  cov_total={int(cc.sum()):,}  "
                  f"psite_total={int(tt.sum()):,}")

    # [3] correlations over the universe, from summary TSVs
    cov_tot = load_col(COV_SUM, "tx_id", "total_coverage", int)
    cov_len = load_col(COV_SUM, "tx_id", "length", int)
    tgt_tot = load_col(TGT_SUM, "tx_id", "total_psites", int)
    # salmon file is 2 columns: tx_id, mean_tpm -> take column 1
    salmon = {}
    with SALMON.open() as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            salmon[f[0]] = float(f[1])

    u = [t for t in universe if t in cov_tot and cov_len.get(t, 0) > 0]
    mean_cov = [cov_tot[t] / cov_len[t] for t in u]
    tpm = [salmon.get(t, 0.0) for t in u]
    r1 = logpearson(mean_cov, tpm)
    p(f"\n[abundance] mean per-nt coverage vs salmon meanTPM (log1p Pearson, n={len(u)}): "
      f"{r1:.3f}  (length-normalized; expect strongly positive)")

    u2 = [t for t in u if t in tgt_tot]
    r2 = logpearson([cov_tot[t] for t in u2], [tgt_tot[t] for t in u2])
    p(f"[coupling] total coverage vs total P-sites (log1p Pearson, n={len(u2)}): {r2:.3f}"
      f"  (expect ~0.7-0.8)")

    OUT.write_text("\n".join(R) + "\n")
    p(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
