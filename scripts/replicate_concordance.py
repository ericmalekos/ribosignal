#!/usr/bin/env python3
"""Replicate-concordance ceiling (Task 22): how reproducible is the Fibroblast Ribo-seq per-nt profile
between two INDEPENDENT measurements? This is the upper bound the model competes against -- no predictor
can match the observed pool better than the observed pool matches an independent re-measurement of itself.

Method: split the 32 Fibroblast Ribo-seq samples into two random half-pools (16 + 16), pool each half's
per-nt P-sites per transcript, and compute the SAME per-transcript profile Pearson the model reports
(pearson over the transcript; plus the 5'UTR window = uORF and 3'UTR window = dORF), stratified by biotype.
Repeated over K random splits (mean +/- sd). Each half is ~half the full 32-sample depth, so the raw
split-half r is Spearman-Brown corrected to full depth: R_full = 2r/(1+r) = the reproducibility of the
full pool the model actually predicts. Report raw (half-depth) + SB (full-depth) next to the model's
profile Pearson (Task 6/9/15: pc ~0.61-0.65, lncRNA ~0.39, uORF ~0.62, dORF ~0.33).

Computed on the model's 36,668-tx packed universe for a 1:1 comparison. cas12a env (numpy + h5py).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
           "biotype_probe/expression_context_human")
NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
SAMPLE_DIR = NEW / "data" / "ribocode_per_tissue" / "Fibroblast"


def pearson(a, b):
    a = a.astype(np.float64); b = b.astype(np.float64)
    if a.size < 3:
        return np.nan
    sa, sb = a.std(), b.std()
    if sa == 0 or sb == 0:
        return np.nan
    return float(((a - a.mean()) * (b - b.mean())).mean() / (sa * sb))


def load_meta():
    tx = [t for t in (NEW / "data/packed/tx_order.txt").read_text().split() if t]
    bt, cds = {}, {}
    with open(NEW / "data/tx2biotype.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t"); ti, gi = h.index("tx_id"), h.index("gene_type")
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); bt[f[ti]] = f[gi]
    with open(NEW / "data/tx2cds.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len", "utr3_len", "has_start_codon")}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            cds[f[ix["tx_id"]]] = (int(f[ix["utr5_len"]]), int(f[ix["cds_len"]]),
                                   int(f[ix["utr3_len"]]), int(f[ix["has_start_codon"]]))
    return tx, bt, cds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5, help="number of random half-splits")
    ap.add_argument("--min_pooled", type=int, default=50, help="min full-32 P-sites for a scorable tx")
    ap.add_argument("--min_half", type=int, default=10, help="min P-sites in EACH half for a valid pearson")
    ap.add_argument("--out", default=str(NEW / "results/replicate_concordance.json"))
    args = ap.parse_args()

    tx_order, bt, cds = load_meta()
    keep_set = set(tx_order)
    samples = sorted(SAMPLE_DIR.glob("*_psites.hd5"))
    assert len(samples) == 32, f"expected 32 Fibroblast samples, found {len(samples)}"

    # axis + index of the packed tx within the 509,650 sample axis (same order across samples)
    with h5py.File(samples[0]) as h0:
        axis = h0["transcript_ids"].asstr()[:]
    pos = {t: i for i, t in enumerate(axis)}
    idx = np.array([pos[t] for t in tx_order if t in pos])
    kept_tx = [t for t in tx_order if t in pos]
    print(f"packed universe {len(tx_order):,}; on sample axis {len(kept_tx):,}", flush=True)

    # load each sample's per-nt for the kept tx ONCE (read is the expensive step); store ragged
    n_tx = len(kept_tx)
    per_sample = []          # list of length 32, each a list of n_tx int32 arrays
    for si, sp in enumerate(samples):
        with h5py.File(sp) as h:
            ps = h["p_sites"]
            assert np.array_equal(h["transcript_ids"].asstr()[:], axis), f"axis mismatch {sp.name}"
            rows = [np.asarray(ps[i], dtype=np.int32) for i in idx]
        per_sample.append(rows)
        print(f"  loaded {si+1}/32 {sp.name}", flush=True)

    lengths = np.array([per_sample[0][i].shape[0] for i in range(n_tx)])
    rng_order = np.arange(32)

    def stratum(t, L):
        u5, cl, u3, hs = cds.get(t, (0, 0, 0, 0))
        return bt.get(t, "other"), (u5 if (hs and u5 >= 8 and u5 <= L) else 0), \
               (u3 if (u3 >= 8 and u3 <= L) else 0)

    results = {"config": vars(args), "n_tx_axis": n_tx, "per_seed": []}
    # accumulate per-seed medians
    agg = {k: [] for k in ("pc_whole", "lncRNA_whole", "pc_uorf5", "pc_dorf3", "pc_cds_codon")}
    for seed in range(args.seeds):
        rng = np.random.RandomState(1000 + seed)
        perm = rng.permutation(rng_order)
        A, B = perm[:16], perm[16:]
        rec = {"pc_whole": [], "lncRNA_whole": [], "pc_uorf5": [], "pc_dorf3": [],
               "pc_cds_codon": []}   # per-codon CDS Shape r = elongation ONLY (periodicity removed)
        for i in range(n_tx):
            t = kept_tx[i]; L = lengths[i]
            ha = np.zeros(L, np.int64); hb = np.zeros(L, np.int64)
            for s in A:
                ha += per_sample[s][i]
            for s in B:
                hb += per_sample[s][i]
            tot = ha.sum() + hb.sum()
            if tot < args.min_pooled or ha.sum() < args.min_half or hb.sum() < args.min_half:
                continue
            biotype, u5, u3 = stratum(t, L)
            r = pearson(ha, hb)
            if r == r:
                if biotype == "protein_coding":
                    rec["pc_whole"].append(r)
                elif biotype == "lncRNA":
                    rec["lncRNA_whole"].append(r)
            # per-codon CDS Shape r (elongation ONLY -- periodicity removed by summing in-frame
            # triplets, matching seq2ribo's target construction; comparable to its Shape r 0.05-0.19)
            if biotype == "protein_coding":
                u5c, clc, u3c, hsc = cds.get(t, (0, 0, 0, 0))
                if hsc and clc >= 90 and u5c + clc <= L:   # valid CDS, >= 30 codons
                    ncod = clc // 3
                    ca = ha[u5c:u5c + ncod * 3].reshape(ncod, 3).sum(1)
                    cb = hb[u5c:u5c + ncod * 3].reshape(ncod, 3).sum(1)
                    if ca.sum() >= args.min_half and cb.sum() >= args.min_half:
                        rc = pearson(ca, cb)
                        if rc == rc:
                            rec["pc_cds_codon"].append(rc)
            # UTR windows (pc only, matching eval_extra uorf_5utr / dorf_3utr)
            if biotype == "protein_coding":
                if u5 and ha[:u5].sum() >= args.min_half and hb[:u5].sum() >= args.min_half:
                    ru = pearson(ha[:u5], hb[:u5])
                    if ru == ru:
                        rec["pc_uorf5"].append(ru)
                if u3 and ha[L - u3:].sum() >= args.min_half and hb[L - u3:].sum() >= args.min_half:
                    rd = pearson(ha[L - u3:], hb[L - u3:])
                    if rd == rd:
                        rec["pc_dorf3"].append(rd)
        med = {k: (float(np.median(v)) if v else np.nan) for k, v in rec.items()}
        n = {k: len(v) for k, v in rec.items()}
        results["per_seed"].append({"seed": seed, "median": med, "n": n})
        for k in agg:
            agg[k].append(med[k])
        print(f"seed {seed}: " + "  ".join(f"{k} {med[k]:.3f} (n={n[k]})" for k in med), flush=True)

    def sb(r):
        return (2 * r) / (1 + r) if r == r else np.nan
    summary = {}
    for k, vals in agg.items():
        v = np.array(vals); m = float(np.nanmean(v)); s = float(np.nanstd(v))
        summary[k] = {"raw_halfdepth_mean": m, "raw_sd": s, "sb_fulldepth": sb(m)}
    results["summary"] = summary
    Path(args.out).write_text(__import__("json").dumps(results, indent=2))

    print("\n=== Replicate-concordance CEILING (Fibroblast, 32-sample split-half) ===")
    # pc_cds_codon = elongation-only (periodicity removed); ref = seq2ribo Shape r 0.05-0.19; our per-codon TBD
    modelnum = {"pc_whole": 0.638, "lncRNA_whole": 0.364, "pc_uorf5": 0.622, "pc_dorf3": 0.329,
                "pc_cds_codon": None}
    print(f"{'class':14} {'raw r (half)':>13} {'SB (full)':>10} {'model':>7} {'model/SB':>9}")
    for k in ("pc_whole", "lncRNA_whole", "pc_uorf5", "pc_dorf3", "pc_cds_codon"):
        m = summary[k]["raw_halfdepth_mean"]; sbf = summary[k]["sb_fulldepth"]; mo = modelnum[k]
        mo_s = f"{mo:>7.3f}" if mo is not None else f"{'TBD':>7}"
        rat = f"{mo/sbf:>8.0%}" if mo is not None else f"{'--':>8}"
        print(f"{k:14} {m:>13.3f} {sbf:>10.3f} {mo_s} {rat}")
    print("  (pc_cds_codon = elongation only, periodicity removed; seq2ribo's comparable Shape r = 0.05-0.19)")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
