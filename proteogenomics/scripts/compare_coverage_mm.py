#!/usr/bin/env python3
"""Compare mm1 vs mm20 RNA-seq coverage on the DoHH2 universe: per-transcript Pearson + total-depth ratio,
to quantify how much the multimap posture changes the model's INPUT. If coverage is near-identical on most
transcripts (differing only on paralog/repeat-heavy ones), the model outputs will be near-identical and a
full mm1 retrain is unwarranted. cas12a env (numpy + h5py)."""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
P = NEW / "proteogenomics/data/mm_proxy"


def load(fp, wanted):
    """Load per-nt coverage for only the `wanted` transcripts (sorted-index read is far faster than
    reading all 509k vlen rows one-by-one)."""
    with h5py.File(fp, "r") as h:
        ids = h["transcript_ids"].asstr()[:]
        pos = {t: i for i, t in enumerate(ids)}
        idx = sorted(pos[t] for t in wanted if t in pos)
        dset = h["coverage"]
        cov = {}
        for i in idx:
            cov[ids[i]] = np.asarray(dset[i], dtype=np.float64)
    return cov


def main():
    f1, f20 = P / "DoHH2_mm1_coverage.hd5", P / "DoHH2_mm20_coverage.hd5"
    for f in (f1, f20):
        if not f.exists():
            sys.exit(f"missing {f}")
    uni = [t for t in (NEW / "proteogenomics/data/DoHH2_pilot/DoHH2_universe_tx.txt").read_text().split() if t]
    print(f"loading coverage for {len(uni):,} universe transcripts...", file=sys.stderr)
    wanted = set(uni)
    c1, c20 = load(f1, wanted), load(f20, wanted)

    tot1 = tot20 = 0.0
    pears, ratios, expressed = [], [], []
    top = []  # (pearson, tx, tot1, tot20)
    for t in uni:
        a, b = c1.get(t), c20.get(t)
        if a is None or b is None or a.shape != b.shape:
            continue
        s1, s20 = a.sum(), b.sum()
        tot1 += s1; tot20 += s20
        if s20 < 50:      # only assess transcripts with real coverage
            continue
        expressed.append(t)
        r = np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else (1.0 if np.allclose(a, b) else 0.0)
        pears.append(r); ratios.append(s1 / s20 if s20 else np.nan)
        top.append((r, t, s1, s20))

    pears = np.array(pears); ratios = np.array(ratios)
    print("\n===== mm1 vs mm20 RNA coverage (DoHH2 universe) =====")
    print(f"transcripts assessed (mm20 depth >=50): {len(pears):,} / {len(uni):,} universe")
    print(f"TOTAL coverage: mm1 {tot1:.3e}  mm20 {tot20:.3e}  ratio mm1/mm20 = {tot1/tot20:.4f}")
    print(f"per-tx Pearson(mm1,mm20): median {np.median(pears):.4f}  mean {pears.mean():.4f}  "
          f"p05 {np.percentile(pears,5):.4f}")
    print(f"  tx with Pearson >= 0.999: {(pears>=0.999).sum():,} ({100*(pears>=0.999).mean():.1f}%)")
    print(f"  tx with Pearson <  0.95 : {(pears<0.95).sum():,} ({100*(pears<0.95).mean():.1f}%)")
    print(f"per-tx depth ratio mm1/mm20: median {np.nanmedian(ratios):.4f}  "
          f"tx where mm20 > 1.5x mm1 (multimap-inflated): {(ratios<0.667).sum():,} "
          f"({100*(ratios<0.667).mean():.1f}%)")
    print("\n-- 10 most-divergent transcripts (lowest Pearson) --")
    for r, t, s1, s20 in sorted(top)[:10]:
        print(f"   {t}  pearson={r:.3f}  depth mm1={s1:.0f} mm20={s20:.0f} (mm20/mm1={s20/max(s1,1):.1f}x)")
    out = P / "coverage_mm_comparison.txt"
    out.write_text(f"tx_assessed={len(pears)} total_ratio_mm1_mm20={tot1/tot20:.4f} "
                   f"median_pearson={np.median(pears):.4f} frac_pearson_ge_0.999={100*(pears>=0.999).mean():.1f}%\n")
    # per-tx dump for the figure (B7)
    tsv = P / "coverage_mm_pertx.tsv"
    with open(tsv, "w") as fh:
        fh.write("tx_id\tpearson\tdepth_mm1\tdepth_mm20\tratio_mm1_mm20\n")
        for r, t, s1, s20 in top:
            fh.write(f"{t}\t{r:.6f}\t{s1:.0f}\t{s20:.0f}\t{s1/s20 if s20 else float('nan'):.6f}\n")
    print(f"\nwrote {out}\nwrote {tsv} ({len(top):,} tx)")


if __name__ == "__main__":
    main()
