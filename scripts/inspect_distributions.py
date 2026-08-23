#!/usr/bin/env python3
"""Empirically characterize the P-site target and RNAseq coverage to ground model design.

Answers the design questions that decide the loss and the input normalization:
  1. Target sparsity: what fraction of transcript-nt carry zero P-sites?
  2. Target overdispersion: variance/mean of per-nt counts on expressed tx (Poisson vs NB)?
  3. 3-nt periodicity: is there a lag-3/6/9 peak in the P-site autocorrelation (the
     biological signal the model must reproduce)?
  4. Coverage dynamic range: justify log1p on the RNAseq input.
  5. Coupling: do per-nt coverage and per-nt P-sites correlate WITHIN a transcript
     (they should not strongly -- RNAseq is ~uniform along the tx, P-sites concentrate on
     the CDS -- so the embedding, not the coverage, must localize translation), while
     per-transcript totals SHOULD correlate (abundance prior)?

Reads: the target hd5, one completed per-sample coverage hd5, the target/coverage summary
TSVs, the universe list. Writes a compact text report to logs/inspect_distributions.txt.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
TARGET = NEW / "data" / "target" / "Fibroblast_psites_pooled.hd5"
TARGET_SUM = NEW / "data" / "target" / "Fibroblast_psites_summary.tsv"
COV1 = NEW / "data" / "rnaseq_coverage" / "per_sample" / "SRR15513236_coverage.hd5"
UNIV = NEW / "data" / "fibroblast_universe.tsv"
SALMON = NEW / "data" / "tpm" / "fibroblast_salmon_mean_tpm.tsv"

SEED = 42
N_SAMPLE = 1500       # transcripts sampled for per-nt stats
MIN_SIGNAL = 50       # min total P-sites to count a tx as "expressed & translated"


def load_universe():
    """tx_id -> (transcript_type, length) for the 36,668 expressed pc+lncRNA universe."""
    d = {}
    with UNIV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ttype = f[col["transcript_type"]] if "transcript_type" in col else "NA"
            d[f[col["tx_id"]]] = ttype
    return d


def id_to_index(h):
    ids = h["transcript_ids"].asstr()[:]
    return ids, {t: i for i, t in enumerate(ids)}


def autocorr_lags(x, max_lag=12):
    """Normalized autocorrelation of a 1-D per-nt track at lags 1..max_lag."""
    x = x.astype(np.float64)
    x = x - x.mean()
    denom = np.dot(x, x)
    if denom <= 0:
        return None
    out = []
    for lag in range(1, max_lag + 1):
        out.append(float(np.dot(x[:-lag], x[lag:]) / denom))
    return np.array(out)


def summarize_tsv_numeric(path, value_col, type_col=None):
    """Return per-tx totals array (+ optional type list) from a summary TSV."""
    vals, types = [], []
    with path.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(hdr)}
        vi = col[value_col]
        ti = col.get(type_col) if type_col else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals.append(float(f[vi]))
            if ti is not None:
                types.append(f[ti])
    return np.array(vals), types


def main():
    rng = np.random.default_rng(SEED)
    univ = load_universe()
    univ_ids = set(univ)
    R = []  # report lines

    def p(s=""):
        R.append(s)
        print(s, file=sys.stderr)

    p("=" * 72)
    p("P-SITE TARGET + RNASEQ COVERAGE DISTRIBUTION REPORT (Fibroblast)")
    p(f"universe: {len(univ_ids):,} expressed pc+lncRNA tx  |  sample N={N_SAMPLE} "
      f"(min total P-sites >= {MIN_SIGNAL})")
    p("=" * 72)

    # ---- per-transcript level from summary TSVs (cheap, whole-universe) ----
    p("\n[1] PER-TRANSCRIPT TOTALS (whole universe, from summary TSVs)")
    tp, ttypes = summarize_tsv_numeric(TARGET_SUM, "total_psites", "transcript_type")
    tsum_ids, _ = summarize_tsv_numeric(TARGET_SUM, "length")  # placeholder for count
    # restrict target-summary rows to universe
    tgt_ids = []
    with TARGET_SUM.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = hdr.index("tx_id")
        for line in fh:
            tgt_ids.append(line.split("\t")[ci])
    tgt_ids = np.array(tgt_ids)
    in_univ = np.array([t in univ_ids for t in tgt_ids])
    tp_u = tp[in_univ]
    p(f"  target total P-sites/tx (universe): median={np.median(tp_u):,.0f}  "
      f"mean={tp_u.mean():,.0f}  p90={np.percentile(tp_u,90):,.0f}  "
      f"max={tp_u.max():,.0f}")
    p(f"  tx with >= {MIN_SIGNAL} P-sites: {(tp_u>=MIN_SIGNAL).sum():,} / {len(tp_u):,} "
      f"({100*(tp_u>=MIN_SIGNAL).mean():.1f}%)")
    p(f"  tx with zero P-sites: {(tp_u==0).sum():,} ({100*(tp_u==0).mean():.1f}%)")

    # ---- sample expressed+translated tx for per-nt stats ----
    elig = tgt_ids[in_univ & (tp >= MIN_SIGNAL)]
    n = min(N_SAMPLE, len(elig))
    sample_ids = rng.choice(elig, size=n, replace=False)
    p(f"\n[2] PER-NUCLEOTIDE STATS (sample of {n} expressed+translated tx)")

    tgt = h5py.File(TARGET, "r")
    tids, tidx = id_to_index(tgt)
    cov = h5py.File(COV1, "r")
    cids, cidx = id_to_index(cov)

    # read sampled rows (sort by target index for chunk locality)
    order = sorted(sample_ids, key=lambda t: tidx[t])
    tgt_ps = tgt["p_sites"]
    cov_ds = cov["coverage"]

    # aggregate per-nt value histogram + moments for target and coverage
    tgt_nt_total = 0
    tgt_zero = 0
    tgt_bins = {"1": 0, "2-5": 0, "6-20": 0, "21+": 0}
    tgt_sum = 0.0
    tgt_sqsum = 0.0
    tgt_max = 0
    cov_zero = 0
    cov_nt_total = 0
    cov_sum = 0.0
    cov_max = 0
    lag_accum = np.zeros(12)
    lag_n = 0
    per_nt_spearman = []   # within-tx coverage vs psite rank corr
    per_tx_total_t = []
    per_tx_total_c = []

    def spearman(a, b):
        if a.size < 8:
            return None
        ra = np.argsort(np.argsort(a))
        rb = np.argsort(np.argsort(b))
        ra = ra - ra.mean()
        rb = rb - rb.mean()
        d = np.sqrt(np.dot(ra, ra) * np.dot(rb, rb))
        return float(np.dot(ra, rb) / d) if d > 0 else None

    for t in order:
        a = np.asarray(tgt_ps[tidx[t]], dtype=np.int64)
        tgt_nt_total += a.size
        tgt_zero += int((a == 0).sum())
        tgt_bins["1"] += int((a == 1).sum())
        tgt_bins["2-5"] += int(((a >= 2) & (a <= 5)).sum())
        tgt_bins["6-20"] += int(((a >= 6) & (a <= 20)).sum())
        tgt_bins["21+"] += int((a >= 21).sum())
        tgt_sum += float(a.sum())
        tgt_sqsum += float(np.dot(a, a))
        tgt_max = max(tgt_max, int(a.max()))
        per_tx_total_t.append(float(a.sum()))
        # periodicity on the P-site track
        ac = autocorr_lags(a, 12)
        if ac is not None:
            lag_accum += ac
            lag_n += 1
        # coverage (may be absent from universe? it isn't -- same set)
        if t in cidx:
            c = np.asarray(cov_ds[cidx[t]], dtype=np.int64)
            cov_nt_total += c.size
            cov_zero += int((c == 0).sum())
            cov_sum += float(c.sum())
            cov_max = max(cov_max, int(c.max()))
            per_tx_total_c.append(float(c.sum()))
            if c.size == a.size:
                sp = spearman(a.astype(np.float64), c.astype(np.float64))
                if sp is not None:
                    per_nt_spearman.append(sp)
        else:
            per_tx_total_c.append(np.nan)

    # target per-nt moments
    mean_t = tgt_sum / tgt_nt_total
    var_t = tgt_sqsum / tgt_nt_total - mean_t ** 2
    p(f"  TARGET per-nt: zero={100*tgt_zero/tgt_nt_total:.1f}%  "
      f"=1:{100*tgt_bins['1']/tgt_nt_total:.1f}%  "
      f"2-5:{100*tgt_bins['2-5']/tgt_nt_total:.1f}%  "
      f"6-20:{100*tgt_bins['6-20']/tgt_nt_total:.1f}%  "
      f"21+:{100*tgt_bins['21+']/tgt_nt_total:.2f}%  max={tgt_max:,}")
    p(f"  TARGET per-nt mean={mean_t:.3f}  var={var_t:.3f}  "
      f"var/mean={var_t/mean_t:.1f}  (>>1 => overdispersed => NB over Poisson)")
    mean_c = cov_sum / cov_nt_total
    p(f"  COVERAGE per-nt: zero={100*cov_zero/cov_nt_total:.1f}%  "
      f"mean={mean_c:.2f}  max={cov_max:,}  (one sample SRR15513236)")

    # periodicity
    if lag_n:
        lag_mean = lag_accum / lag_n
        p(f"\n[3] 3-nt PERIODICITY (mean P-site autocorr over {lag_n} tx, lags 1-12)")
        p("  lag:  " + "  ".join(f"{i:>5d}" for i in range(1, 13)))
        p("  ac:   " + "  ".join(f"{v:5.2f}" for v in lag_mean))
        frame = lag_mean[[2, 5, 8, 11]].mean()      # lags 3,6,9,12
        off = lag_mean[[0, 1, 3, 4, 6, 7]].mean()   # lags 1,2,4,5,7,8
        p(f"  mean ac at lags 3/6/9/12 = {frame:.3f}  vs off-frame (1/2/4/5/7/8) = "
          f"{off:.3f}  -> periodicity {'PRESENT' if frame>off+0.02 else 'weak'}")

    # coupling
    sp_arr = np.array(per_nt_spearman)
    tt = np.array(per_tx_total_t)
    tc = np.array(per_tx_total_c)
    ok = ~np.isnan(tc)
    # per-tx total corr in log space
    lt = np.log1p(tt[ok])
    lc = np.log1p(tc[ok])
    lt -= lt.mean()
    lc -= lc.mean()
    dd = np.sqrt(np.dot(lt, lt) * np.dot(lc, lc))
    tx_corr = float(np.dot(lt, lc) / dd) if dd > 0 else float("nan")
    p("\n[4] COUPLING coverage <-> P-sites")
    p(f"  WITHIN-tx per-nt Spearman: median={np.median(sp_arr):.3f} "
      f"IQR=[{np.percentile(sp_arr,25):.3f},{np.percentile(sp_arr,75):.3f}]  "
      f"(low => coverage does NOT localize translation; embedding must)")
    p(f"  ACROSS-tx total(P-sites) vs total(coverage), log1p Pearson r={tx_corr:.3f}  "
      f"(high => coverage is a good abundance prior)")

    # salmon cross-check of per-tx coverage total
    sal = {}
    with SALMON.open() as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sal[f[0]] = float(f[1])
    xs, ys = [], []
    for i, t in enumerate(order):
        if ok[i] and t in sal and sal[t] > 0:
            xs.append(np.log1p(sal[t]))
            ys.append(np.log1p(tc[i]))
    if len(xs) > 10:
        xs = np.array(xs)
        ys = np.array(ys)
        xs -= xs.mean()
        ys -= ys.mean()
        dd = np.sqrt(np.dot(xs, xs) * np.dot(ys, ys))
        r = float(np.dot(xs, ys) / dd) if dd > 0 else float("nan")
        p(f"  SANITY total(coverage) vs salmon meanTPM, log1p Pearson r={r:.3f} "
          f"(n={len(xs)}; should be strongly positive)")

    p("\n" + "=" * 72)
    out = NEW / "logs" / "inspect_distributions.txt"
    out.write_text("\n".join(R) + "\n")
    print(f"\nwrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
