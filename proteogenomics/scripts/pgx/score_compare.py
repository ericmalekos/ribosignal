#!/usr/bin/env python3
"""Post-FDR score comparison between database arms (`pgx` diagnostic).

Everything here is measured AFTER 1% class-specific FDR filtering. Pre-FDR score distributions are
dominated by the random-match floor and say little about the identifications actually reported.

THE CONFOUND, AND WHY THE NAIVE TEST IS MEANINGLESS. Each arm's FDR threshold is set by its own
target-decoy competition, and a larger novel space raises it: on BMDM the cuts are 20.20 (model),
28.10 (null_atg), 28.67 (null_nc). Post-FDR scores are therefore FLOORED AT DIFFERENT VALUES, so
"which arm has higher passing scores" is answered by the thresholds before any data is consulted --
the nulls must win, whatever the quality of their identifications.

Three tests are reported, in increasing order of how much they can be trusted:

  raw        passing PSMs per arm, compared directly. Shown only to quantify the confound; the
             difference it reports is largely the difference in thresholds.
  common     all arms re-thresholded at max(cut) so the floor is shared. Removes the confound but
             discards the model's advantage of a lower threshold, so it is conservative.
  matched    the same PEPTIDE where it passes in both arms, compared pairwise. Confound-free: the
             floor cannot differ for a fixed peptide present in both sets. This is the test that
             answers "is a given identification scored differently depending on the database".

Mann-Whitney U (two-sided, nonparametric) for the unpaired tests; Wilcoxon signed-rank for matched.
Effect size is reported as the rank-biserial correlation, because a p-value on tens of thousands of
spectra will be small regardless of whether the difference matters.

  python -m pgx.score_compare --search-root <dir> --arms model_poisson,null_atg,null_nc
"""
from __future__ import annotations

import argparse
import itertools
import statistics as st
from pathlib import Path

from .report import FDR, cut, load_rank1


def passing_psms(rows, fdr=FDR):
    """[(peptide, score)] for novel PSMs above the arm's own class-specific FDR threshold."""
    c = cut(rows, {"novel_t"}, {"novel_d"}, fdr)
    if c is None:
        return [], None
    return [(p, h) for (p, h, k, _a) in rows if k == "novel_t" and h >= c], c


def mwu(a, b):
    """Mann-Whitney U -> (p, rank-biserial r). |r| ~ 0.1 small, 0.3 medium, 0.5 large."""
    from scipy.stats import mannwhitneyu
    if not a or not b:
        return float("nan"), float("nan")
    res = mannwhitneyu(a, b, alternative="two-sided")
    r = 2.0 * res.statistic / (len(a) * len(b)) - 1.0
    return float(res.pvalue), float(r)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search-root", required=True)
    ap.add_argument("--arms", default="model_poisson,null_atg,null_nc")
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    root = Path(a.search_root)
    arms = [x.strip() for x in a.arms.split(",") if x.strip()]
    P, C = {}, {}
    for arm in arms:
        rows = load_rank1(root / arm)
        P[arm], C[arm] = passing_psms(rows, a.fdr)
        print(f"  {arm}: {len(P[arm]):,} passing novel PSMs, cut = {C[arm]:.2f}")

    L = ["", "### Post-FDR novel-PSM scores (1% class-specific FDR)", "",
         "| arm | passing PSMs | uniq peptides | FDR cut | median | p25 | p75 |",
         "|---|--:|--:|--:|--:|--:|--:|"]
    for arm in arms:
        s = sorted(h for _p, h in P[arm])
        if not s:
            continue
        peps = len({p for p, _h in P[arm]})
        L.append(f"| {arm} | {len(s):,} | {peps:,} | {C[arm]:.2f} | {st.median(s):.2f} | "
                 f"{s[len(s) // 4]:.2f} | {s[3 * len(s) // 4]:.2f} |")

    L += ["", "**raw** (confounded by the differing thresholds -- shown to size the confound):", "",
          "| pair | median A | median B | p (MWU) | rank-biserial |", "|---|--:|--:|--:|--:|"]
    for x, y in itertools.combinations(arms, 2):
        sa = [h for _p, h in P[x]]
        sb = [h for _p, h in P[y]]
        if not sa or not sb:
            continue
        p, r = mwu(sa, sb)
        L.append(f"| {x} vs {y} | {st.median(sa):.2f} | {st.median(sb):.2f} | {p:.3g} | {r:+.3f} |")

    floor = max(c for c in C.values() if c is not None)
    L += ["", f"**common floor** (all arms re-thresholded at {floor:.2f}, the strictest cut):", "",
          "| pair | n A | n B | median A | median B | p (MWU) | rank-biserial |",
          "|---|--:|--:|--:|--:|--:|--:|"]
    for x, y in itertools.combinations(arms, 2):
        sa = [h for _p, h in P[x] if h >= floor]
        sb = [h for _p, h in P[y] if h >= floor]
        if not sa or not sb:
            continue
        p, r = mwu(sa, sb)
        L.append(f"| {x} vs {y} | {len(sa):,} | {len(sb):,} | {st.median(sa):.2f} | "
                 f"{st.median(sb):.2f} | {p:.3g} | {r:+.3f} |")

    L += ["", "**matched peptides** (same peptide passing in both arms; confound-free):", "",
          "| pair | n peptides | median A | median B | median diff | p (Wilcoxon) |",
          "|---|--:|--:|--:|--:|--:|"]
    from scipy.stats import wilcoxon
    best = {arm: {} for arm in arms}
    for arm in arms:
        for p_, h in P[arm]:
            best[arm][p_] = max(best[arm].get(p_, 0.0), h)
    for x, y in itertools.combinations(arms, 2):
        shared = sorted(set(best[x]) & set(best[y]))
        if len(shared) < 5:
            L.append(f"| {x} vs {y} | {len(shared)} | -- | -- | -- | too few |")
            continue
        va = [best[x][p_] for p_ in shared]
        vb = [best[y][p_] for p_ in shared]
        d = [i - j for i, j in zip(va, vb)]
        pv = wilcoxon(va, vb).pvalue if any(abs(z) > 1e-12 for z in d) else 1.0
        L.append(f"| {x} vs {y} | {len(shared):,} | {st.median(va):.2f} | {st.median(vb):.2f} | "
                 f"{st.median(d):+.3f} | {pv:.3g} |")

    # Per-SPECTRUM matched test: the same scan, matched to the same peptide, in both arms. This is
    # the only fully paired comparison -- max-per-peptide (above) is biased by how many spectra each
    # arm passes, since a maximum over more draws runs higher.
    from .report import load_rank1_spec
    S = {}
    for arm in arms:
        rows = load_rank1_spec(root / arm)
        c = C[arm]
        S[arm] = {sid: (pep, h) for (sid, pep, h, k) in rows
                  if k == "novel_t" and c is not None and h >= c}
    L += ["", "**matched SPECTRA** (same scan AND same peptide passing in both arms; fully paired):",
          "", "| pair | shared scans | same peptide | diff peptide | median diff | p (Wilcoxon) |",
          "|---|--:|--:|--:|--:|--:|"]
    for x, y in itertools.combinations(arms, 2):
        both = set(S[x]) & set(S[y])
        same = [s for s in both if S[x][s][0] == S[y][s][0]]
        diff = len(both) - len(same)
        if len(same) < 5:
            L.append(f"| {x} vs {y} | {len(both):,} | {len(same):,} | {diff:,} | -- | too few |")
            continue
        va = [S[x][s][1] for s in same]
        vb = [S[y][s][1] for s in same]
        d = [i - j for i, j in zip(va, vb)]
        pv = wilcoxon(va, vb).pvalue if any(abs(z) > 1e-12 for z in d) else 1.0
        L.append(f"| {x} vs {y} | {len(both):,} | {len(same):,} | {diff:,} | {st.median(d):+.4f} | "
                 f"{pv:.3g} |")

    L += ["",
          "Read the matched-SPECTRA rows first: they are the only fully paired, confound-free "
          "comparison. A zero median difference there means hyperscore is a property of the "
          "spectrum-peptide match and not of the database, so every apparent score advantage in "
          "the raw rows is the threshold difference re-expressed.",
          "The `diff peptide` column counts scans that pass in both arms but were assigned "
          "DIFFERENT peptides -- rank-1 competition changing the winner as the search space grows.",
          "The max-per-peptide `matched peptides` rows above are retained for continuity but are "
          "biased when the arms pass different numbers of spectra per peptide; prefer the "
          "per-spectrum rows.",
          "Effect size (rank-biserial) is reported alongside p because at these n a p-value is "
          "small whether or not the difference is meaningful."]
    md = "\n".join(L) + "\n"
    print(md)
    if a.out:
        Path(a.out).write_text(md)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
