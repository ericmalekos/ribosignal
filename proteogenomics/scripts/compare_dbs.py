#!/usr/bin/env python3
"""Generalized N-way search-DB comparison (A549 or HBL-1): novel peptides at class-specific FDR per DB,
plus pairwise discovery-rate ratios and the novel-peptide SET OVERLAP between DBs. Novel = a rank-1 PSM
whose peptide maps ONLY to nuORF| ORFs (not shared with any canonical protein); its FDR is estimated
against REV_nuORF| decoys ONLY (honest per-class denominator), not the global pool -- per the db1_tradeoff
"crazy FDR" lesson. The overlap block is the biology: for HBL-1, does the MODEL's ORF selection recover
immunopeptidome peptides that the authors' experimental PRICE Ribo-seq DB missed, and vice versa? stdlib."""
from __future__ import annotations

import argparse
import csv
import glob
from itertools import combinations
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def load_rank1(pilot, db, prefix="search_"):
    out = []
    for t in sorted(glob.glob(str(NEW / f"proteogenomics/data/{pilot}/{prefix}{db}/*.tsv"))):
        for r in csv.DictReader(open(t), delimiter="\t"):
            if r.get("hit_rank") != "1":
                continue
            accs = tuple(a for a in r["proteins"].split(",") if a)
            out.append((r["peptide"], float(r["hyperscore"]), accs))
    return out


def klass(accs):
    tgt = [a for a in accs if not a.startswith("REV_")]
    if tgt:
        return "novel_t" if all(a.startswith("nuORF|") for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith("REV_nuORF|") for a in accs) else "other_d"


def class_fdr_novel(psms, fdr=0.01):
    """Return (accepted novel peptide set, n_target, n_decoy) at class-specific FDR."""
    sel = [(p, h, klass(a)) for (p, h, a) in psms]
    sel = [(p, h, k) for (p, h, k) in sel if k in ("novel_t", "novel_d")]
    sel.sort(key=lambda x: -x[1])
    t = d = 0; last_ok_h = None
    for _, h, k in sel:
        if k == "novel_d":
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= fdr:
            last_ok_h = h
    if last_ok_h is None:
        return set(), 0, 0
    keep = set(); t = d = 0
    for p, h, k in sel:
        if h < last_ok_h:
            break
        if k == "novel_d":
            d += 1
        else:
            t += 1; keep.add(p)
    return keep, t, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default="HBL1_pilot")
    ap.add_argument("--dbs", nargs="+", default=["canonical", "model", "null", "price"])
    ap.add_argument("--search_prefix", default="search_", help="e.g. search_ or search_fresh_")
    ap.add_argument("--db_dir", default="db", help="dir holding <db>_class_map.tsv (db or db_fresh)")
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or str(NEW / f"proteogenomics/data/{args.pilot}/db_comparison.md")

    dbsize = {}
    for db in args.dbs:
        cm = NEW / f"proteogenomics/data/{args.pilot}/{args.db_dir}/{db}_class_map.tsv"
        dbsize[db] = (sum(1 for _ in open(cm)) - 1) if cm.exists() else 0

    L = [f"# {args.pilot}: N-way search-DB comparison (model vs null vs PRICE)\n",
         f"Novel peptides at {args.fdr:.0%} class-specific FDR (rank-1, novel-class target-decoy). "
         "Novel = peptide maps ONLY to non-canonical ORFs.\n",
         "| DB | novel ORFs in DB | novel peptides @1% class FDR | discovery rate (pep/1k ORFs) | canonical target PSMs |",
         "|----|---:|---:|---:|---:|"]
    keepset = {}; canon_psm = {}
    for db in args.dbs:
        psms = load_rank1(args.pilot, db, args.search_prefix)
        canon_t = sum(1 for (_, _, a) in psms if klass(a) == "canon_t")
        canon_psm[db] = canon_t
        if db == "canonical":
            L.append(f"| canonical (baseline) | 0 | -- | -- | {canon_t:,} |")
            keepset[db] = set(); continue
        keep, t, d = class_fdr_novel(psms, args.fdr)
        keepset[db] = keep
        n_orf = dbsize.get(db, 0)
        rate = 1000 * len(keep) / n_orf if n_orf else 0
        L.append(f"| {db} | {n_orf:,} | **{len(keep):,}** | {rate:.2f} | {canon_t:,} |")

    # ---- pairwise discovery-rate ratios ----
    L.append("\n## Discovery-rate ratios (novel pep / 1k ORFs)\n")
    novel_dbs = [d for d in args.dbs if d != "canonical"]
    rate_of = {d: (1000 * len(keepset[d]) / dbsize[d] if dbsize.get(d) else 0) for d in novel_dbs}
    for a, b in combinations(novel_dbs, 2):
        if rate_of[b]:
            L.append(f"- **{a} / {b}: {rate_of[a] / rate_of[b]:.2f}x**  "
                     f"({a} {len(keepset[a]):,} pep / {dbsize[a]:,} ORFs vs "
                     f"{b} {len(keepset[b]):,} pep / {dbsize[b]:,} ORFs)")

    # ---- novel-peptide SET OVERLAP (the biology) ----
    L.append("\n## Novel-peptide set overlap (unique peptide sequences @1% class FDR)\n")
    for a, b in combinations(novel_dbs, 2):
        A, B = keepset[a], keepset[b]
        L.append(f"- **{a} vs {b}**: {a}-only {len(A - B):,} | shared {len(A & B):,} | "
                 f"{b}-only {len(B - A):,}  (|{a}|={len(A):,}, |{b}|={len(B):,})")
    if "model" in keepset and "price" in keepset:
        m, p = keepset["model"], keepset["price"]
        L.append(f"\n### Headline (model vs authors' PRICE Ribo-seq DB)\n")
        L.append(f"- Model recovers **{len(m - p):,}** immunopeptidome peptides PRICE missed; "
                 f"PRICE finds {len(p - m):,} the model missed; {len(m & p):,} shared.")
        L.append(f"- Model DB is {dbsize.get('null',0)/max(dbsize['model'],1):.1f}x smaller than null and "
                 f"{dbsize.get('price',0)/max(dbsize['model'],1):.1f}x smaller than PRICE, yet finds "
                 f"{len(m):,} novel peptides (PRICE {len(p):,}).")

    md = "\n".join(L) + "\n"
    Path(out).write_text(md)
    print(md); print(f"wrote {out}")


if __name__ == "__main__":
    main()
