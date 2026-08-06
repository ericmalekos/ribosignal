#!/usr/bin/env python3
"""Aggregate the model-vs-null churn analysis across all 12 macrophage populations. Per population, per DB
(canonical/model/null): rank-1 PSMs pooled over the 18 fractions, novel peptides at 1% class-specific FDR
(novel = maps ONLY to nuORF| targets, scored vs REV_nuORF| decoys ONLY), and canonical target PSMs. Emits
the standing-rule DB-tradeoff table: novel discovery (model vs null) PLUS PC-churn (change in canonical PSMs)
and net total PSM -- not novel discovery alone. Reuses the exact FDR logic of compare_model_vs_null.py.
cas12a env (stdlib).

--tag namespaces the search/db trees per model checkpoint (matching macro_dumps.sbatch TAG). Because the
canonical DB carries zero model content, and db_null's SEQUENCE set is determined by the universe FASTA
rather than the checkpoint, those two searches are reusable across checkpoints: --reuse-tag names the tag
to fall back to when a (population, db) is absent under --tag. Every fallback is reported in the output so
a reused search can never be mistaken for a fresh one."""
import argparse
import csv
import glob
import sys
from pathlib import Path

PD = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
          "riboseq_signal_model/proteogenomics/data/macrophage_tissue")
FDR = 0.01
TAG = ""
REUSE_TAG = ""
REUSED = []


def search_dir(pop, db):
    """Resolve a (population, db) search dir, falling back to REUSE_TAG and recording the fallback."""
    d = PD / f"search{TAG}" / pop / db
    if TAG and not any(d.glob("*.tsv")):
        alt = PD / f"search{REUSE_TAG}" / pop / db
        if any(alt.glob("*.tsv")):
            REUSED.append((pop, db))
            return alt
    return d


def db_dir(pop):
    d = PD / f"db{TAG}" / pop
    if TAG and not d.is_dir():
        return PD / f"db{REUSE_TAG}" / pop
    return d


def klass(accs):
    tgt = [a for a in accs if not a.startswith("REV_")]
    if tgt:
        return "novel_t" if all(a.startswith("nuORF|") for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith("REV_nuORF|") for a in accs) else "other_d"


def load_rank1(pop, db):
    out = []
    for t in sorted(glob.glob(str(search_dir(pop, db) / "*.tsv"))):
        with open(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("hit_rank") != "1":
                    continue
                accs = tuple(a for a in r["proteins"].split(",") if a)
                out.append((r["peptide"], float(r["hyperscore"]), accs))
    return out


def _cut(ranked, tset, dset):
    """Descending-score target-decoy cut; returns the last score where decoy/target <= FDR."""
    t = d = 0
    last_ok = None
    for h, k in ranked:
        if k in dset:
            d += 1
        elif k in tset:
            t += 1
        if t > 0 and d / t <= FDR:
            last_ok = h
    return last_ok


def analyze(psms):
    """-> (novel peptides @class-specific FDR, canonical PSMs @GLOBAL FDR, net PSMs @GLOBAL FDR).

    TWO DIFFERENT FDRs ON PURPOSE (see Task 52 / the DB-tradeoff standing rule):

    * canonical + net use GLOBAL 1% FDR (all targets vs all decoys). This is what the search reports
      and it keeps the target-decoy competition intact. Class-specific FDR on the canonical side is
      confounded: a large novel space cannibalises canonical DECOYS faster than canonical targets
      (BMDM: other_d -29.5% vs canon_t -13.2%), which deflates the estimated canonical FDR and lets a
      junk database appear to GAIN canonical identifications.
    * novel uses CLASS-SPECIFIC FDR (novel targets vs REV_nuORF| decoys only) -- the honest
      denominator; a global FDR inflates non-canonical discovery ~10-13x.

    Every returned column is FDR-filtered. Earlier versions returned RAW rank-1 counts for canonical
    and net, which measured how many spectra a database ABSORBED rather than identified, and therefore
    scaled with database size: the 2.36M-sequence null "won" net PSM on 1,592,864 novel matches that
    yielded only 374 FDR-surviving peptides, while its apparent -1,480,032 canonical churn shrank to
    a 0.27% spread once filtered.
    """
    sel = [(p, h, klass(a)) for (p, h, a) in psms]

    # --- global 1% FDR for the canonical and net columns ---
    ranked_g = sorted([(h, k) for (_, h, k) in sel], key=lambda x: -x[0])
    g = _cut(ranked_g, {"canon_t", "novel_t"}, {"novel_d", "other_d"})
    if g is None:
        canon = total_t = 0
    else:
        canon = sum(1 for (_, h, k) in sel if k == "canon_t" and h >= g)
        total_t = sum(1 for (_, h, k) in sel if k in ("canon_t", "novel_t") and h >= g)
    nov = sorted([(h, k) for (_, h, k) in sel if k in ("novel_t", "novel_d")], key=lambda x: -x[0])
    t = d = 0
    last_ok = None
    for h, k in nov:
        if k == "novel_d":
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= FDR:
            last_ok = h
    keep = set()
    if last_ok is not None:
        for (p, h, a) in psms:
            if h >= last_ok and klass(a) == "novel_t":
                keep.add(p)
    return len(keep), canon, total_t


def dbsize(pop, db):
    cm = db_dir(pop) / f"{db}_class_map.tsv"
    return sum(1 for _ in open(cm)) - 1 if cm.exists() else 0


def main():
    global TAG, REUSE_TAG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="", help="model tag, e.g. _mamba4_union (default: original layout)")
    ap.add_argument("--reuse-tag", default="",
                    help="tag to fall back to for (pop, db) absent under --tag; canonical and null "
                         "searches are checkpoint-independent and reusable")
    ap.add_argument("--title", default=None, help="H1 title for the emitted markdown")
    ap.add_argument("--out", default=None, help="output .md path (default macro_model_vs_null<tag>.md)")
    a = ap.parse_args()
    TAG, REUSE_TAG = a.tag, a.reuse_tag

    root = PD / f"search{TAG}"
    if not root.is_dir():
        sys.exit(f"no search tree at {root}")
    pops = sorted(p.name for p in root.iterdir() if p.is_dir())
    rows = []
    for pop in pops:
        d = {db: analyze(load_rank1(pop, db)) for db in ("canonical", "model", "null")}
        rows.append((pop, d, dbsize(pop, "model"), dbsize(pop, "null")))
        print(f"  {pop} done", file=sys.stderr)

    title = a.title or ("# Macrophage proteogenomics: model-selected vs null ORF DB, per population"
                       + (f" ({TAG.lstrip('_')})" if TAG else ""))
    out = [title + "\n",
           "ALL columns are FDR-filtered (Task 52). Novel peptides use 1% CLASS-SPECIFIC FDR (novel = maps "
           "ONLY to model/null-called ORFs, scored against REV_nuORF| decoys only). Canonical and net PSMs use "
           "1% GLOBAL FDR (all targets vs all decoys), which keeps the target-decoy competition intact -- "
           "class-specific FDR on the canonical side is confounded by decoy cannibalisation, and RAW rank-1 "
           "counts measure spectra absorbed rather than identified and scale with DB size. Expect canonical to "
           "be near-FLAT across DBs (~0.3% spread); a large apparent churn means the filter was skipped.\n",
           "| population | model novel | null novel | model/null rate | canon base | canon model | "
           "ΔPC model | canon null | ΔPC null | net PSM model | net PSM null |",
           "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
    tot = {"mn": 0, "nn": 0, "cb": 0, "cm": 0, "cn": 0, "netm": 0, "netn": 0}
    for pop, d, msz, nsz in rows:
        mn, _, mnet = d["model"]
        nn, _, nnet = d["null"]
        cb = d["canonical"][1]
        cm = d["model"][1]
        cn = d["null"][1]
        rate = (mn / msz) / (nn / nsz) if (nn and msz and nsz and mn) else float("nan")
        out.append(f"| {pop} | {mn} | {nn} | {rate:.2f}x | {cb:,} | {cm:,} | {cm - cb:+,} | "
                   f"{cn:,} | {cn - cb:+,} | {mnet:,} | {nnet:,} |")
        tot["mn"] += mn; tot["nn"] += nn; tot["cb"] += cb; tot["cm"] += cm
        tot["cn"] += cn; tot["netm"] += mnet; tot["netn"] += nnet
    out.append(f"| **TOTAL (12 pop)** | **{tot['mn']}** | **{tot['nn']}** | | {tot['cb']:,} | "
               f"{tot['cm']:,} | {tot['cm'] - tot['cb']:+,} | {tot['cn']:,} | {tot['cn'] - tot['cb']:+,} | "
               f"{tot['netm']:,} | {tot['netn']:,} |")
    out += ["\n## Verdict\n",
            f"- Novel peptides (1% class FDR): model **{tot['mn']}** vs null **{tot['nn']}** across 12 populations.",
            f"- PC-churn: canonical PSMs {tot['cb']:,} (baseline) -> {tot['cm']:,} (model, {tot['cm']-tot['cb']:+,}) "
            f"-> {tot['cn']:,} (null, {tot['cn']-tot['cb']:+,}). The null DB displaces "
            f"{abs(tot['cn']-tot['cb'])-abs(tot['cm']-tot['cb']):,} MORE canonical PSMs than the model DB.",
            f"- Net PSM: model {tot['netm']:,} vs null {tot['netn']:,}."]
    if REUSED:
        by_db = {}
        for pop, db in REUSED:
            by_db.setdefault(db, []).append(pop)
        out.append("\n## Reused searches\n")
        for db, ps in sorted(by_db.items()):
            out.append(f"- `{db}` reused from tag `{REUSE_TAG or '<original>'}` for {len(ps)} population(s): "
                       f"{', '.join(sorted(ps))}. Valid because this DB's sequence content does not depend "
                       f"on the model checkpoint.")
    md = "\n".join(out) + "\n"
    dest = Path(a.out) if a.out else PD / f"macro_model_vs_null{TAG}.md"
    dest.write_text(md)
    print(md)
    print(f"wrote {dest}", file=sys.stderr)


if __name__ == "__main__":
    main()
