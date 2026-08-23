#!/usr/bin/env python3
"""How many novel proteogenomic DISCOVERIES would a 30-AA ORF floor destroy?

THE INCONSISTENCY THIS MEASURES. The ORF-call track filters at 90 nt = exactly 30 aa
(`compare_dropin_calls.build_loader`, applied identically to model and real Ribo-seq calls, no
override anywhere). The proteogenomics track builds its search databases at `--min-aa 7`
("detectability"). So the two tracks report on DIFFERENT ORF populations, and nothing declares it:
22% of the model's macrophage database and ~65% of the naive nulls sit below the ORF-call floor.

Counting short DB ENTRIES is the wrong measure for deciding what to do, because short ORFs are also
the least likely to yield a detectable peptide. What matters is how many novel PSMs and peptides
that actually PASSED FDR come from them.

THE UNIT THAT DECIDES IT IS "SHORT-ONLY". A novel peptide can be supported by several ORFs. It is
only LOST under a 30-aa floor if EVERY supporting ORF is under 30 aa; if any one of them is >= 30 aa
the peptide survives the filter. Counting peptides that merely touch a short ORF would overstate the
cost, sometimes badly, since short ORFs are frequently nested inside longer ones.

FDR is the pipeline's own `pgx.report.novel_at_class_fdr` (1% class-specific, novel targets vs
`REV_nuORF|` decoys only), imported rather than reimplemented, and `by_pep` gives each surviving
peptide its supporting novel accessions directly.

  python3 proteogenomics/scripts/orf_length_audit.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pgx.report import load_rank1, novel_at_class_fdr  # noqa: E402

NEW = HERE.parents[1]
MIN_AA = 30          # the ORF-call track's floor: build_loader min_len=90 nt, ORF_length excludes stop
XS = NEW / "proteogenomics/data/macrophage_tissue/pgx_xsubtype"


NOVEL_PREFIX = "nuORF|"
DECOY_PREFIX = "REV_"


def fasta_lengths(fa):
    """-> {accession: protein length in aa} for NOVEL TARGET (nuORF|) entries only.

    The search fastas hold targets AND decoys AND the full canonical proteome -- 117,518 entries for
    1,533 novel ORFs. An unfiltered dict made the reported `db_short` fraction 2.0% where the true
    novel-target fraction is 22.0%, because ~115k canonical and decoy entries diluted it. The
    short-only peptide counts were never affected (they look up specific nuORF| accessions), but the
    two columns then disagreed inside one table, which is worse than either being wrong alone.
    """
    out, name, seq = {}, None, []

    def flush():
        if name and name.startswith(NOVEL_PREFIX):
            out[name] = len("".join(seq))

    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
    flush()
    return out


def audit(search_dir, db_fa, min_aa=MIN_AA, fdr=0.01):
    if not pathlib.Path(db_fa).exists():
        return None
    L = fasta_lengths(db_fa)
    n_psm, peps, accs, by_pep = novel_at_class_fdr(load_rank1(search_dir), fdr)
    if not peps:
        return dict(novel_psms=0, novel_peptides=0, short_only_peptides=0, kept_peptides=0,
                    unmapped_peptides=0, db_novel=len(L),
                    db_short=sum(1 for v in L.values() if v < min_aa))
    short_only = kept = unmapped = 0
    for p, a in by_pep.items():
        lens = [L[x] for x in a if x in L]
        if not lens:
            unmapped += 1                      # accession not in this DB fasta -- reported, not hidden
        elif max(lens) < min_aa:
            short_only += 1                    # every supporting ORF is short -> LOST at the floor
        else:
            kept += 1
    return dict(novel_psms=n_psm, novel_peptides=len(peps), short_only_peptides=short_only,
                kept_peptides=kept, unmapped_peptides=unmapped, db_novel=len(L),
                db_short=sum(1 for v in L.values() if v < min_aa))


def human(a):
    """The 5 human pgx datasets behind F2b / D12 / D13, same measurement, own data."""
    DS = ["A549", "HBL1", "SUDHL4", "DoHH2", "THP1"]
    ARMS = ["model_poisson", "model_standard", "null_atg", "null_nc", "cpat", "cpc2"]
    rows = []
    print(f"human pgx datasets (mamba4), floor = {a.min_aa} aa\n")
    print(f"  {'dataset':<10}{'arm':<18}{'PSMs':>6}{'pept':>6}{'SHORT-ONLY':>12}{'kept':>6}"
          f"{'unmapped':>10}{'DBnovel<30aa':>13}")
    for ds in DS:
        root = NEW / f"proteogenomics/data/{ds}_pilot/pgx_mamba4"
        for arm in ARMS:
            sd, db = root / "search" / arm, root / "db" / f"db_{arm}.fasta"
            if not sd.exists() or not db.exists():
                continue
            r = audit(sd, db, a.min_aa)
            if r is None:
                continue
            r.update(dataset=ds, arm=arm)
            rows.append(r)
            pct = 100 * r["short_only_peptides"] / max(1, r["novel_peptides"])
            print(f"  {ds:<10}{arm:<18}{r['novel_psms']:>6}{r['novel_peptides']:>6}"
                  f"{r['short_only_peptides']:>7} ({pct:>4.1f}%){r['kept_peptides']:>6}"
                  f"{r['unmapped_peptides']:>10}"
                  f"{100*r['db_short']/max(1,r['db_novel']):>12.1f}%")
    print()
    for arm in ARMS:
        sub = [r for r in rows if r["arm"] == arm]
        if not sub:
            continue
        tp = sum(r["novel_peptides"] for r in sub)
        ts = sum(r["short_only_peptides"] for r in sub)
        print(f"  TOTAL {arm:<18} peptides {tp:>5,}  short-only {ts:>5,} ({100*ts/max(1,tp):.1f}%)"
              f"  unmapped {sum(r['unmapped_peptides'] for r in sub)}")
    out = str(NEW / "proteogenomics/data/orf_length_audit_human.json")
    pathlib.Path(out).write_text(json.dumps({"min_aa": a.min_aa, "rows": rows}, indent=2))
    print(f"\n  wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-aa", type=int, default=MIN_AA)
    ap.add_argument("--out", default=str(XS / "orf_length_audit.json"))
    ap.add_argument("--human", action="store_true",
                    help="audit the 5 HUMAN pgx datasets (F2b / D12 / D13) instead of the mouse "
                         "macrophage sweep. Their conclusions must not rest on a macrophage "
                         "measurement extrapolated across species and assay.")
    a = ap.parse_args()

    if a.human:
        return human(a)
    arms = ["model_predicted", "model_ribocode_bmdm_nt", "null_atg"]
    rows = []
    pops = sorted(p.name for p in (XS / "frozen").iterdir() if p.is_dir())
    print(f"macrophage sweep: {len(pops)} populations, floor = {a.min_aa} aa\n")
    print(f"  {'population':<17}{'arm':<24}{'PSMs':>6}{'pept':>6}{'SHORT-ONLY':>12}"
          f"{'kept':>6}{'unmapped':>10}{'DBnovel<30aa':>13}")
    for pop in pops:
        for arm in arms:
            sd = XS / "frozen" / pop / arm
            db = XS / "db_percall" / pop / f"db_{arm}.fasta"
            if not db.exists():
                db = XS / "db" / f"db_{arm}.fasta"
            if not db.exists():
                db = XS / "db_null" / f"db_{arm}.fasta"
            r = audit(sd, db, a.min_aa)
            if r is None:
                print(f"  {pop:<17}{arm:<24}  (no DB fasta found -- skipped, not zeroed)")
                continue
            r.update(population=pop, arm=arm)
            rows.append(r)
            pct = 100 * r["short_only_peptides"] / max(1, r["novel_peptides"])
            print(f"  {pop:<17}{arm:<24}{r['novel_psms']:>6}{r['novel_peptides']:>6}"
                  f"{r['short_only_peptides']:>7} ({pct:>4.1f}%){r['kept_peptides']:>6}"
                  f"{r['unmapped_peptides']:>10}"
                  f"{100*r['db_short']/max(1,r['db_novel']):>12.1f}%")

    print()
    for arm in arms:
        sub = [r for r in rows if r["arm"] == arm]
        if not sub:
            continue
        tp = sum(r["novel_peptides"] for r in sub)
        ts = sum(r["short_only_peptides"] for r in sub)
        tu = sum(r["unmapped_peptides"] for r in sub)
        print(f"  TOTAL {arm:<24} peptides {tp:>5,}  short-only {ts:>5,} ({100*ts/max(1,tp):.1f}%)"
              f"  unmapped {tu}")
    pathlib.Path(a.out).write_text(json.dumps({"min_aa": a.min_aa, "rows": rows}, indent=2))
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
