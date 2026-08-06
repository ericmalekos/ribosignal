#!/usr/bin/env python3
"""Class-specific FDR on MS2Rescore-rescored A549 PSMs (basic + ms2pip CID), mirroring compare_dbs.py but
using the rescored mokapot score instead of raw hyperscore. A NOVEL peptide is a rank-1 PSM whose peptide
maps ONLY to nuORF| ORFs; its FDR is estimated against REV_nuORF| decoys ONLY (honest per-class denominator).
Reads proteogenomics/data/A549_pilot/rescore/<db>.psms.tsv (pooled over the 12 fractions). Reports, per DB,
novel peptides at 1% class FDR + canonical PSMs, and the model-vs-null verdict -- the rescored analogue of
db_comparison_fresh.md. cas12a env (stdlib)."""
from __future__ import annotations

import argparse
import ast
import csv
import glob
import re
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
csv.field_size_limit(1 << 24)


def bare(peptidoform: str) -> str:
    """strip charge (/2) and [mod] annotations -> bare AA sequence (match the raw-hyperscore peptide count)."""
    s = peptidoform.split("/")[0]
    s = re.sub(r"\[[^\]]*\]", "", s)
    return re.sub(r"[^A-Z]", "", s.upper())


def load_rank1(db, pilot="A549_pilot", rescore_dir="rescore"):
    """rank-1 rescored PSMs pooled over the per-fraction files: (bare_peptide, rescored_score, accs-tuple).
    rescore_dir lets a caller point at a per-seed rescore dir (e.g. rescore_seed3) for the mokapot
    stochasticity study; defaults to the canonical `rescore`."""
    files = sorted(glob.glob(str(NEW / f"proteogenomics/data/{pilot}/{rescore_dir}/{db}/*.psms.tsv")))
    if not files:
        sys.exit(f"no rescored psms.tsv for db={db} under {pilot}/{rescore_dir}/{db}/")
    out = []
    for f in files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r.get("rank") not in ("1", "1.0"):
                continue
            try:
                accs = tuple(ast.literal_eval(r["protein_list"]))
            except Exception:
                accs = tuple(a.strip(" '\"[]") for a in r["protein_list"].split(",") if a.strip(" '\"[]"))
            out.append((bare(r["peptidoform"]), float(r["score"]), accs))
    return out


def klass(accs):
    tgt = [a for a in accs if not a.startswith("REV_")]
    if tgt:
        return "novel_t" if all(a.startswith("nuORF|") for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith("REV_nuORF|") for a in accs) else "other_d"


def class_fdr_novel(psms, fdr=0.01):
    sel = [(p, s, klass(a)) for (p, s, a) in psms]
    sel = [(p, s, k) for (p, s, k) in sel if k in ("novel_t", "novel_d")]
    sel.sort(key=lambda x: -x[1])
    t = d = 0; last_ok = None
    for _, s, k in sel:
        if k == "novel_d":
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= fdr:
            last_ok = s
    if last_ok is None:
        return set(), 0, 0
    keep = set(); t = d = 0
    for p, s, k in sel:
        if s < last_ok:
            break
        if k == "novel_d":
            d += 1
        else:
            t += 1; keep.add(p)
    return keep, t, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default="A549_pilot")
    ap.add_argument("--db_dir", default="db_fresh", help="dir with <db>_class_map.tsv (db_fresh or db)")
    ap.add_argument("--dbs", nargs="+", default=["canonical", "model", "null"])
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--out", default=None)
    ap.add_argument("--rescore_dir", default="rescore", help="per-seed rescore dir for the stochasticity study")
    args = ap.parse_args()
    out = args.out or str(NEW / f"proteogenomics/data/{args.pilot}/db_comparison_rescored.md")

    dbsize = {}
    for db in args.dbs:
        cm = NEW / f"proteogenomics/data/{args.pilot}/{args.db_dir}/{db}_class_map.tsv"
        dbsize[db] = (sum(1 for _ in open(cm)) - 1) if cm.exists() else 0

    L = [f"# {args.pilot}: MS2Rescore-rescored novel-peptide comparison\n",
         f"Novel peptides at {args.fdr:.0%} class-specific FDR on the **rescored mokapot score** (rank-1, "
         "novel-class target-decoy).\n",
         "| DB | novel ORFs | novel peptides @1% class FDR (rescored) | rate (pep/1k ORFs) | canonical PSMs |",
         "|----|---:|---:|---:|---:|"]
    keep = {}; canon = {}
    for db in args.dbs:
        psms = load_rank1(db, args.pilot, args.rescore_dir)
        canon[db] = sum(1 for (_, _, a) in psms if klass(a) == "canon_t")
        if db == "canonical":
            L.append(f"| canonical (baseline) | 0 | -- | -- | {canon[db]:,} |")
            keep[db] = set(); continue
        k, t, d = class_fdr_novel(psms, args.fdr)
        keep[db] = k
        n = dbsize.get(db, 0)
        L.append(f"| {db} | {n:,} | **{len(k):,}** | {1000*len(k)/n if n else 0:.2f} | {canon[db]:,} |")

    nov_dbs = [d for d in args.dbs if d != "canonical"]
    if len(nov_dbs) >= 2:
        L.append("\n## Discovery-rate ratios + overlap (rescored)\n")
        from itertools import combinations
        for a, b in combinations(nov_dbs, 2):
            ra = 1000*len(keep[a])/dbsize[a] if dbsize.get(a) else 0
            rb = 1000*len(keep[b])/dbsize[b] if dbsize.get(b) else 0
            ratio = f"{ra/rb:.2f}x" if rb else "n/a"
            L.append(f"- **{a}/{b}: {ratio}** ({a} {len(keep[a])} pep / {dbsize[a]:,} ORFs vs "
                     f"{b} {len(keep[b])} / {dbsize[b]:,}); overlap {a}-only {len(keep[a]-keep[b])} | "
                     f"shared {len(keep[a]&keep[b])} | {b}-only {len(keep[b]-keep[a])}")
        cb = canon.get("canonical", 0)
        L.append("\n## Canonical PSM churn (rescored)\n- canonical " + f"{cb:,}" + " -> " +
                 " ; ".join(f"{d} {canon[d]:,} ({100*(cb-canon[d])/max(cb,1):+.1f}%)" for d in nov_dbs))

    md = "\n".join(L) + "\n"
    Path(out).write_text(md)
    print(md); print(f"wrote {out}")


if __name__ == "__main__":
    main()
