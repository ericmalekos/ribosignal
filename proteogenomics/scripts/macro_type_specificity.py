#!/usr/bin/env python3
"""Is the per-type macrophage design actually buying type specificity?

The pipeline already builds ONE database per macrophage population, from that population's own
RNA-seq: 12 universes (18,102-26,145 tx), 12 candidate sets, 12 model DBs, and each population's
mzML is searched only against its own DB. What has never been asked is whether that specificity is
REAL or merely nominal -- if the 12 universes were ~identical, per-type DBs would be bookkeeping.

Three questions, none of which the aggregate TOTAL row can answer:

  1. UNIVERSE divergence -- how much of each population's expressed transcriptome is shared with all
     12 (core) vs unique to it? Driven by the tissue RNA-seq, independent of the model.
  2. DATABASE divergence -- same for the model-selected novel ORF sequences. This is where the model
     and the tissue RNA-seq compound.
  3. PEPTIDE specificity -- of the novel peptides each population actually detects at 1% class FDR,
     how many are seen ONLY in that population? This is the payoff: a novel peptide found solely in
     Kupffer cells, from a transcript expressed in liver, is a tissue-resident microprotein
     candidate; one found in all 12 is a housekeeping ORF the per-type design did not need.

Note (3) is confounded by detection depth -- a peptide can be "unique" simply because the other
populations lacked the spectra. The report therefore gives, alongside uniqueness, how many of the
OTHER populations even carried that sequence in their search DB. A peptide unique to Kupffer but
present in 11 other DBs is a genuine biological signal; one unique to Kupffer because only Kupffer's
DB contained it is a database artefact.

cas12a env (stdlib).
"""
import argparse
import csv
import glob
import sys
from collections import Counter, defaultdict
from pathlib import Path

PD = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
          "riboseq_signal_model/proteogenomics/data/macrophage_tissue")
FDR = 0.01


def klass(accs):
    tgt = [a for a in accs if not a.startswith("REV_")]
    if tgt:
        return "novel_t" if all(a.startswith("nuORF|") for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith("REV_nuORF|") for a in accs) else "other_d"


def novel_peps(search_root, pop, db="model"):
    """Unique novel peptides at 1% CLASS-SPECIFIC FDR (novel vs REV_nuORF| only)."""
    psms = []
    for t in sorted(glob.glob(str(search_root / pop / db / "*.tsv"))):
        with open(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("hit_rank") != "1":
                    continue
                accs = tuple(a for a in r["proteins"].split(",") if a)
                psms.append((r["peptide"], float(r["hyperscore"]), klass(accs)))
    ranked = sorted([(h, k) for (_, h, k) in psms if k in ("novel_t", "novel_d")], key=lambda x: -x[0])
    t = d = 0
    last = None
    for h, k in ranked:
        if k == "novel_d":
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= FDR:
            last = h
    if last is None:
        return set()
    return {p for (p, h, k) in psms if k == "novel_t" and h >= last}


def db_seqs(db_root, pop):
    """Novel target SEQUENCES in this population's model DB (for the coverage control)."""
    out, name, buf = set(), None, []
    f = db_root / pop / "db_model.fasta"
    if not f.exists():
        return out
    for ln in open(f):
        if ln.startswith(">"):
            if name and name.startswith("nuORF|"):
                out.add("".join(buf))
            name = ln[1:].strip().split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if name and name.startswith("nuORF|"):
        out.add("".join(buf))
    return out


def jaccard_report(label, sets, pops):
    n = len(pops)
    core = set.intersection(*[sets[p] for p in pops]) if all(sets[p] for p in pops) else set()
    union = set.union(*[sets[p] for p in pops])
    print(f"\n  === {label} ===")
    print(f"  union across {n} populations: {len(union):,}   core (in all {n}): {len(core):,} "
          f"({len(core) / max(len(union), 1) * 100:.1f}% of union)")
    print(f"  {'population':<18}{'total':>9}{'core':>9}{'unique':>9}{'% unique':>10}")
    print("  " + "-" * 55)
    cnt = Counter()
    for p in pops:
        for x in sets[p]:
            cnt[x] += 1
    for p in pops:
        uniq = sum(1 for x in sets[p] if cnt[x] == 1)
        tot = len(sets[p])
        print(f"  {p:<18}{tot:>9,}{len(core):>9,}{uniq:>9,}{uniq / max(tot, 1) * 100:>9.1f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="_mamba4_union_cal", help="model tag for db/search trees")
    ap.add_argument("--skip-peptides", action="store_true", help="universe/DB divergence only (fast)")
    a = ap.parse_args()

    pops = sorted(p.name.replace("_universe_tx.txt", "")
                  for p in (PD / "universes").glob("*_universe_tx.txt"))
    print(f"macrophage type-specificity report   tag={a.tag}   populations={len(pops)}")

    # 1. universe divergence (tissue RNA-seq only, model-independent)
    uni = {p: {ln.strip() for ln in open(PD / "universes" / f"{p}_universe_tx.txt") if ln.strip()}
           for p in pops}
    jaccard_report("1. EXPRESSED UNIVERSE (from each type's own RNA-seq)", uni, pops)

    # 2. model-DB divergence
    dbs = {p: db_seqs(PD / f"db{a.tag}", p) for p in pops}
    if all(dbs[p] for p in pops):
        jaccard_report("2. MODEL-SELECTED NOVEL ORF SEQUENCES", dbs, pops)

    if a.skip_peptides:
        return

    # 3. detected-peptide specificity + the DB-coverage control
    root = PD / f"search{a.tag}"
    pep = {}
    for p in pops:
        pep[p] = novel_peps(root, p)
        print(f"    scored {p}: {len(pep[p])} novel peptides @1% class FDR", file=sys.stderr)
    jaccard_report("3. DETECTED NOVEL PEPTIDES (1% class-specific FDR)", pep, pops)

    cnt = Counter()
    for p in pops:
        for x in pep[p]:
            cnt[x] += 1
    print("\n  === 3b. Are 'unique' peptides biology or database artefact? ===")
    print("  For peptides detected in exactly ONE population: how many OTHER populations carried that")
    print("  sequence in their search DB but did not detect it? High = real biological specificity.")
    print(f"  {'population':<18}{'unique pep':>11}{'searchable elsewhere':>22}{'DB-limited':>12}")
    print("  " + "-" * 63)
    tot_bio = tot_art = 0
    for p in pops:
        u = [x for x in pep[p] if cnt[x] == 1]
        bio = art = 0
        for s in u:
            others = sum(1 for q in pops if q != p and s in dbs.get(q, set()))
            if others >= 1:
                bio += 1
            else:
                art += 1
        tot_bio += bio
        tot_art += art
        print(f"  {p:<18}{len(u):>11}{bio:>22}{art:>12}")
    print("  " + "-" * 63)
    print(f"  {'TOTAL':<18}{tot_bio + tot_art:>11}{tot_bio:>22}{tot_art:>12}")
    if tot_bio + tot_art:
        print(f"\n  {tot_bio / (tot_bio + tot_art) * 100:.0f}% of population-unique peptides were "
              f"searchable in at least one other population's DB and still not found there,\n"
              f"  i.e. genuine detection specificity rather than a database-coverage artefact.")


if __name__ == "__main__":
    main()
