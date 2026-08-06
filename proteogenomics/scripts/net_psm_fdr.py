#!/usr/bin/env python3
"""FDR-filtered net PSM -- correcting the DB-tradeoff table's `net PSM` column.

THE BUG. `macro_churn_aggregate.analyze()` FDR-filters only the novel-PEPTIDE column. Its `canon`
and `total_t` (reported as "canon" and "net PSM") count EVERY rank-1 PSM with no FDR filter at all:

    canon   = sum(1 for (_,_,k) in sel if k == "canon_t")
    total_t = sum(1 for (_,_,k) in sel if k in ("canon_t","novel_t"))

So `net PSM` measures how many spectra a database ABSORBED, not how many identifications it made,
and a bigger database absorbs more spectra by chance. That inverts the column's meaning: on the
macrophage set the null DB (2.36M sequences) "wins" net PSM with 1,592,864 novel-target rank-1 PSMs
of which only 374 unique peptides survive 1% class FDR -- 4,259 raw matches per real peptide, against
399 for the CDS-anchored model DB. Worse, the null's novel gain (+1,592,864) is ~93% accounted for by
its canonical loss (-1,480,032): it is mostly RELABELLING canonical identifications as junk-novel.

THE FIX. Filter both components at 1% FDR before summing, each against its own decoy class:
    canonical : canon_t  vs  other_d      (REV_ accessions that are NOT REV_nuORF|)
    novel     : novel_t  vs  novel_d      (REV_nuORF| only -- the class-specific denominator)
    net       : the sum of the two filtered counts
Reported in PSMs (comparable to the old column) and in unique peptides (the honest unit).

Same rank-1 pooling, same klass() logic, same 1% target as the aggregator, so the only change is
that the canonical and net columns are now filtered like the novel column always was.

cas12a env (stdlib). Reads the existing search output -- no new MSFragger required.
"""
import argparse
import csv
import glob
import sys
from pathlib import Path

PD = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
          "riboseq_signal_model/proteogenomics/data/macrophage_tissue")
FDR = 0.01


def klass(accs):
    tgt = [a for a in accs if not a.startswith("REV_")]
    if tgt:
        return "novel_t" if all(a.startswith("nuORF|") for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith("REV_nuORF|") for a in accs) else "other_d"


def load_rank1(search_root, pop, db):
    out = []
    for t in sorted(glob.glob(str(search_root / pop / db / "*.tsv"))):
        with open(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("hit_rank") != "1":
                    continue
                accs = tuple(a for a in r["proteins"].split(",") if a)
                out.append((r["peptide"], float(r["hyperscore"]), accs))
    return out


def fdr_cut(psms, tkey, dkey):
    """Descending-hyperscore target-decoy cut; returns (n_psm, n_peptide) at <=FDR."""
    lab = [(p, h, klass(a)) for (p, h, a) in psms]
    ranked = sorted([(h, k) for (_, h, k) in lab if k in (tkey, dkey)], key=lambda x: -x[0])
    t = d = 0
    last_ok = None
    for h, k in ranked:
        if k == dkey:
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= FDR:
            last_ok = h
    if last_ok is None:
        return 0, 0
    n_psm = sum(1 for (_, h, k) in lab if k == tkey and h >= last_ok)
    peps = {p for (p, h, k) in lab if k == tkey and h >= last_ok}
    return n_psm, len(peps)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="", help="model tag, e.g. _mamba4_union_cal")
    ap.add_argument("--reuse-tag", default="", help="fallback tag for canonical/null")
    ap.add_argument("--dbs", default="canonical,model,null")
    a = ap.parse_args()

    root = PD / f"search{a.tag}"
    alt = PD / f"search{a.reuse_tag}"
    pops = sorted(p.name for p in root.iterdir() if p.is_dir())
    print(f"tag={a.tag or '<orig>'}  populations={len(pops)}  FDR={FDR}\n")
    print(f"  {'db':<10}{'canon PSM':>12}{'novel PSM':>11}{'NET PSM':>12}"
          f"{'canon pep':>11}{'novel pep':>11}{'NET pep':>10}")
    print("  " + "-" * 77)
    tot = {}
    for db in a.dbs.split(","):
        cP = cQ = nP = nQ = 0
        for pop in pops:
            src = root if any((root / pop / db).glob("*.tsv")) else alt
            psms = load_rank1(src, pop, db)
            if not psms:
                continue
            p1, q1 = fdr_cut(psms, "canon_t", "other_d")
            p2, q2 = fdr_cut(psms, "novel_t", "novel_d")
            cP += p1; cQ += q1; nP += p2; nQ += q2
        tot[db] = (cP, nP, cP + nP, cQ, nQ, cQ + nQ)
        print(f"  {db:<10}{cP:>12,}{nP:>11,}{cP + nP:>12,}{cQ:>11,}{nQ:>11,}{cQ + nQ:>10,}")
    if "canonical" in tot and "model" in tot:
        b, m = tot["canonical"], tot["model"]
        print(f"\n  vs canonical-only baseline:  net PSM {m[2] - b[2]:+,}   net peptides {m[5] - b[5]:+,}")
        if "null" in tot:
            n = tot["null"]
            print(f"  null vs baseline:            net PSM {n[2] - b[2]:+,}   net peptides {n[5] - b[5]:+,}")
            print(f"  MODEL - NULL:                net PSM {m[2] - n[2]:+,}   net peptides {m[5] - n[5]:+,}")


if __name__ == "__main__":
    main()
