#!/usr/bin/env python3
"""The UC-A headline: does the MODEL-selected ORF DB find novel peptides at a higher rate than the NULL
(all-candidate) DB? Class-specific FDR (rank-1, novel-class target-decoy) per db1_tradeoff/the "crazy FDR"
lesson: a NOVEL peptide is a PSM whose peptide maps ONLY to nuORF| ORFs (not shared with any canonical
protein); its FDR is estimated against REV_nuORF| decoys ONLY (the honest per-class denominator), NOT the
global decoy pool. Reports, per DB: novel peptides at 1% class FDR, novel target-DB size, discovery rate,
and canonical PSMs (PC-churn: does adding novel ORFs displace canonical IDs?). cas12a env (stdlib)."""
from __future__ import annotations

import argparse
import csv
import glob
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def load_rank1(db, glob_fmt):
    """rank-1 PSMs pooled over the fractions: (peptide, hyperscore, accs). glob_fmt must contain '{db}'."""
    out = []
    for t in sorted(glob.glob(glob_fmt.format(db=db))):
        for r in csv.DictReader(open(t), delimiter="\t"):
            if r.get("hit_rank") != "1":
                continue
            accs = [a for a in r["proteins"].split(",") if a]
            out.append((r["peptide"], float(r["hyperscore"]), tuple(accs)))
    return out


def klass(accs):
    """-> 'novel_t' | 'canon_t' | 'novel_d' | 'other_d'. novel = maps ONLY to nuORF| (target) or
    REV_nuORF| (decoy); canonical/other otherwise."""
    tgt = [a for a in accs if not a.startswith("REV_")]
    if tgt:                                           # TARGET PSM
        return "novel_t" if all(a.startswith("nuORF|") for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith("REV_nuORF|") for a in accs) else "other_d"


def class_fdr_novel(psms, fdr=0.01):
    """Novel peptides at class-specific FDR: rank novel_t + novel_d by hyperscore, cut where
    cum_decoy/cum_target > fdr, return the set of accepted unique novel peptide sequences."""
    sel = [(p, h, k) for (p, h, a) in psms for k in (klass(a),) if k in ("novel_t", "novel_d")]
    sel.sort(key=lambda x: -x[1])
    t = d = 0; keep = set(); last_ok_h = None
    for p, h, k in sel:
        if k == "novel_d":
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= fdr:
            last_ok_h = h
    # second pass: accept target PSMs with hyperscore >= the last score meeting FDR
    if last_ok_h is None:
        return set(), 0, 0
    t = d = 0
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
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--search-glob",
                    default=str(NEW / "proteogenomics/data/A549_pilot/search_{db}/*.tsv"),
                    help="glob for the per-DB MSFragger .tsv; MUST contain '{db}' (canonical/model/null)")
    ap.add_argument("--db-dir", default=str(NEW / "proteogenomics/data/A549_pilot/db"),
                    help="dir holding {model,null}_class_map.tsv")
    ap.add_argument("--title", default="A549 proteogenomics: model-selected vs null (3-frame) ORF DB")
    ap.add_argument("--out", default=str(NEW / "proteogenomics/data/A549_pilot/model_vs_null.md"))
    args = ap.parse_args()

    # novel-target DB sizes (unique protein seqs) from the class maps
    dbsize = {}
    for db in ("model", "null"):
        cm = Path(args.db_dir) / f"{db}_class_map.tsv"
        dbsize[db] = sum(1 for _ in open(cm)) - 1 if cm.exists() else 0

    lines = [f"# {args.title}\n",
             f"Novel peptides at {args.fdr:.0%} class-specific FDR (rank-1, novel-class target-decoy). "
             "Novel = peptide maps ONLY to model-called ORFs, not shared with any canonical protein.\n",
             "| DB | novel ORFs in DB | novel peptides @1% class FDR | discovery rate (pep/1k ORFs) | "
             "canonical target PSMs |", "|----|---:|---:|---:|---:|"]
    res = {}
    for db in ("canonical", "model", "null"):
        psms = load_rank1(db, args.search_glob)
        canon_t = sum(1 for (_, _, a) in psms if klass(a) == "canon_t")
        if db == "canonical":
            lines.append(f"| canonical (baseline) | 0 | -- | -- | {canon_t:,} |")
            res[db] = (0, 0, canon_t)
            continue
        keep, t, d = class_fdr_novel(psms, args.fdr)
        n_orf = dbsize[db]
        rate = 1000 * len(keep) / n_orf if n_orf else 0
        lines.append(f"| {db} | {n_orf:,} | **{len(keep):,}** | {rate:.2f} | {canon_t:,} |")
        res[db] = (n_orf, len(keep), canon_t)

    m_orf, m_pep, m_can = res["model"]; n_orf, n_pep, n_can = res["null"]
    lines.append("\n## Verdict\n")
    lines.append(f"- Model DB: {m_pep:,} novel peptides from {m_orf:,} ORFs "
                 f"({1000*m_pep/m_orf:.2f}/1k). Null DB: {n_pep:,} from {n_orf:,} ({1000*n_pep/n_orf:.2f}/1k).")
    if n_pep and m_pep:
        lines.append(f"- **Discovery-rate ratio (model/null): {(m_pep/m_orf)/(n_pep/n_orf):.2f}x** "
                     f"(>1 = the model's selection is enriched for detectable novel peptides).")
        lines.append(f"- Absolute novel peptides: model {m_pep:,} vs null {n_pep:,} "
                     f"({'model finds MORE despite a '+str(round(n_orf/m_orf,1))+'x smaller DB' if m_pep>=n_pep else 'null finds more (bigger net), but at '+str(round((n_orf/m_orf),1))+'x the DB'}).")
    lines.append(f"- PC-churn (canonical target PSMs): canonical-DB {res['canonical'][2]:,}, "
                 f"model {m_can:,}, null {n_can:,} (drop = novel ORFs stealing canonical IDs).")

    md = "\n".join(lines) + "\n"
    Path(args.out).write_text(md)
    print(md)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
