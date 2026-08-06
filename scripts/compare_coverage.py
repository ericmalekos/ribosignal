#!/usr/bin/env python3
"""Compare two per-sample coverage hd5 on a probe set (validation, e.g. pysam vs deeptools).

Reads only a probe set of transcript rows from each file (housekeeping isoforms + a
random sample), so it is fast even on 3.7 GB files. Reports per-transcript total
Pearson r, max relative difference, ACTB/GAPDH totals, and exact per-nt agreement on a
few transcripts. Datasets named `coverage` or `p_sites` are both handled.

Usage: compare_coverage.py <A.hd5> <B.hd5> [tx2biotype.tsv]
"""
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def decode(arr):
    return np.array([t.decode() if isinstance(t, bytes | np.bytes_) else str(t)
                     for t in arr])


def dset_key(h):
    return "coverage" if "coverage" in h else "p_sites"


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    tx2b = sys.argv[3] if len(sys.argv) > 3 else str(NEW / "data" / "tx2biotype.tsv")

    ha = h5py.File(a_path, "r")
    hb = h5py.File(b_path, "r")
    ida = decode(ha["transcript_ids"][:])
    idb = decode(hb["transcript_ids"][:])
    ka, kb = dset_key(ha), dset_key(hb)
    posa = {t: i for i, t in enumerate(ida)}
    posb = {t: i for i, t in enumerate(idb)}
    common = [t for t in ida if t in posb]
    print(f"A={a_path} ({ka}, {len(ida):,} tx)")
    print(f"B={b_path} ({kb}, {len(idb):,} tx)")
    print(f"common tx: {len(common):,}")

    # housekeeping probes
    hk = {}
    for line in open(tx2b):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2 and f[2] in ("ACTB", "GAPDH", "FN1", "THBS1"):
            hk.setdefault(f[2], []).append(f[0])

    rng = np.random.default_rng(0)
    common_arr = np.array(common)
    sample = list(rng.choice(common_arr, size=min(3000, len(common)), replace=False))
    probe = sorted(set(sample) | {t for v in hk.values() for t in v if t in posb})

    ia = sorted(posa[t] for t in probe)
    ib_map = {t: posb[t] for t in probe}
    rows_a = {ida[i]: np.asarray(ha[ka][i], dtype=np.int64) for i in ia}
    ib = sorted(ib_map.values())
    idb_at = {i: idb[i] for i in ib}
    rows_b = {}
    for i in ib:
        rows_b[idb_at[i]] = np.asarray(hb[kb][i], dtype=np.int64)

    ta, tb, lendiff = [], [], 0
    exact = 0
    for t in probe:
        a = rows_a.get(t)
        b = rows_b.get(t)
        if a is None or b is None:
            continue
        if a.shape[0] != b.shape[0]:
            lendiff += 1
            continue
        ta.append(int(a.sum()))
        tb.append(int(b.sum()))
        if np.array_equal(a, b):
            exact += 1
    ta, tb = np.array(ta, float), np.array(tb, float)
    r = np.corrcoef(ta, tb)[0, 1] if len(ta) > 2 else float("nan")
    denom = np.where(ta + tb > 0, (ta + tb) / 2, 1)
    reldiff = np.abs(ta - tb) / denom
    print(f"\nprobe transcripts compared: {len(ta):,}  (length-mismatch: {lendiff})")
    print(f"per-tx total Pearson r : {r:.6f}")
    print(f"per-tx total: exact-array-equal {exact}/{len(ta)}  "
          f"median|rel diff| {np.median(reldiff):.4g}  max {np.max(reldiff):.4g}")
    print(f"grand total (probe): A={int(ta.sum()):,}  B={int(tb.sum()):,}  "
          f"ratio {tb.sum() / max(ta.sum(), 1):.4f}")

    print("\nhousekeeping gene totals (A vs B):")
    for g, txs in sorted(hk.items()):
        sa = sum(int(rows_a[t].sum()) for t in txs if t in rows_a)
        sb = sum(int(rows_b[t].sum()) for t in txs if t in rows_b)
        ratio = sb / max(sa, 1)
        print(f"  {g:<6} isoforms={len(txs):<3} A={sa:>12,}  B={sb:>12,}  ratio {ratio:.4f}")


if __name__ == "__main__":
    main()
