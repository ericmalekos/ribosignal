#!/usr/bin/env python3
"""THE standard proteomics output table (standing rule, 2026-08-01).

Every proteogenomics DB-comparison table in this project reports, per population, EXACTLY these columns:

  novel        unique novel peptides from the MODEL DB, 1% CLASS-SPECIFIC FDR
  null         unique novel peptides from the NULL DB,  1% CLASS-SPECIFIC FDR
  DB seqs      model DB novel target sequences   (discovery is not comparable without the denominator)
  null seqs    null  DB novel target sequences
  canon base   canonical peptides/PSMs from the canonical-ONLY DB, 1% GLOBAL FDR (the baseline)
  dPC model    change in canonical vs that baseline when the model DB is used
  dPC null     change in canonical vs that baseline when the null  DB is used
  ncStart      model peptides absent from the NULL amino-acid space (see below)

WHY dPC FOR BOTH ARMS. Adding sequences to a search perturbs the target-decoy competition, so canonical
identifications can be LOST purely because the database grew -- independent of whether the added ORFs
are real. Reporting dPC for the model arm alone hides whether that cost is intrinsic to adding ORFs or
specific to adding BAD ones. Both arms must be shown so the reader can separate a database-size effect
from a database-quality effect. FDR control is the central methodological hazard in proteogenomics and
this table is where it is made visible rather than assumed away.

WHICH FDR WHERE (Task 52 -- do not mix filtered and unfiltered columns):
  canonical / dPC : GLOBAL 1% FDR (all targets vs all decoys) -- what the search reports, and it keeps
                    the target-decoy competition intact. Class-specific FDR on the canonical side is
                    confounded: a large novel space cannibalises canonical DECOYS faster than canonical
                    targets, deflating the estimate and letting a junk DB appear to GAIN canonical IDs.
  novel           : CLASS-SPECIFIC 1% FDR (novel targets vs REV_nuORF| decoys ONLY). A global FDR
                    inflates non-canonical discovery ~10-13x.
Raw rank-1 counts are never reported: they measure spectra ABSORBED, which scales with DB size.

THE ncStart COLUMN. Model peptides that pass FDR and are found NOWHERE in the null DB's amino-acid
space (substring scan against every null target sequence, since a peptide is a fragment of a protein).
These are discoveries a naive enumeration could not have made at any threshold. The expected mechanism
is a non-canonical start codon: the null enumerates AUG..stop only, so a CUG/GUG/ACG-initiated ORF
extends the protein N-terminally past the first in-frame AUG and yields tryptic peptides that exist in
no AUG-only database.

  *** As of 2026-08-01 this column is STRUCTURALLY ZERO and the table says so. ***
  `enumerate_score_orfs.candidate_orfs()` matches AUG only:
      is_atg = (c0 == A) & (c1 == T) & (c2 == G)
  verified empirically: 400,000/400,000 BMDM candidate ORFs begin ATG. Both DBs are built from that one
  enumeration and db_model is a strict SUBSET of db_null, so no model peptide can be absent from the
  null space. Making this column meaningful requires extending enumeration to non-AUG starts (the ORF
  track already SCORES them: CUG 0.5, GUG/ACG 0.35, UUG 0.3, AUA 0.25 ...), with only the model
  selecting them. The column is computed rather than hard-coded so it becomes live the moment that
  lands, and a non-zero value can never pass unnoticed.

cas12a env (stdlib).
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


def load_rank1(d):
    out = []
    for t in sorted(glob.glob(str(d / "*.tsv"))):
        with open(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("hit_rank") != "1":
                    continue
                out.append((r["peptide"], float(r["hyperscore"]),
                            klass(tuple(a for a in r["proteins"].split(",") if a))))
    return out


def cut(lab, tset, dset):
    ranked = sorted([(h, k) for (_, h, k) in lab if k in tset | dset], key=lambda x: -x[0])
    t = d = 0
    last = None
    for h, k in ranked:
        if k in dset:
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= FDR:
            last = h
    return last


def canon_at_global_fdr(lab):
    g = cut(lab, {"canon_t", "novel_t"}, {"novel_d", "other_d"})
    if g is None:
        return 0
    return len({p for (p, h, k) in lab if k == "canon_t" and h >= g})


def novel_at_class_fdr(lab):
    c = cut(lab, {"novel_t"}, {"novel_d"})
    if c is None:
        return set()
    return {p for (p, h, k) in lab if k == "novel_t" and h >= c}


def null_aa_space(fasta):
    """One concatenated string of every NULL target sequence, for O(1)-ish substring containment."""
    parts, name, buf = [], None, []
    for ln in open(fasta):
        if ln.startswith(">"):
            if name is not None and not name.startswith("REV_"):
                parts.append("".join(buf))
            name = ln[1:].strip().split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if name is not None and not name.startswith("REV_"):
        parts.append("".join(buf))
    return "|".join(parts)


def dbsize(class_map):
    return sum(1 for _ in open(class_map)) - 1 if Path(class_map).exists() else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="", help="model tag, e.g. _mamba4_union_cal")
    ap.add_argument("--reuse-tag", default="", help="fallback tag for canonical/null searches")
    ap.add_argument("--out", default=None, help="write markdown here as well as stdout")
    a = ap.parse_args()

    root = PD / f"search{a.tag}"
    alt = PD / f"search{a.reuse_tag}"
    dbr = PD / f"db{a.tag}"
    dba = PD / f"db{a.reuse_tag}"
    pops = sorted(p.name for p in root.iterdir() if p.is_dir())

    hdr = (f"| population | novel | null | DB seqs | null seqs | canon base | dPC model | "
           f"dPC null | ncStart |")
    sep = "|---|--:|--:|--:|--:|--:|--:|--:|--:|"
    lines = [hdr, sep]
    T = dict(nv=0, nl=0, ds=0, ns=0, cb=0, dm=0, dn=0, nc=0)
    for pop in pops:
        def sd(db):
            d = root / pop / db
            return d if any(d.glob("*.tsv")) else alt / pop / db
        lab_c = load_rank1(sd("canonical"))
        lab_m = load_rank1(sd("model"))
        lab_n = load_rank1(sd("null"))
        cb = canon_at_global_fdr(lab_c)
        cm = canon_at_global_fdr(lab_m)
        cn = canon_at_global_fdr(lab_n)
        pm = novel_at_class_fdr(lab_m)
        pn = novel_at_class_fdr(lab_n)
        dm = dbr / pop / "model_class_map.tsv"
        ds = dbsize(dm if dm.exists() else dba / pop / "model_class_map.tsv")
        nsz = dbsize(dba / pop / "null_class_map.tsv")

        # ncStart: FDR-passing model peptides absent from the NULL amino-acid space
        nullfa = dba / pop / "db_null.fasta"
        nc = 0
        if pm and nullfa.exists():
            space = null_aa_space(nullfa)
            nc = sum(1 for p in pm if p not in space)

        lines.append(f"| {pop} | {len(pm)} | {len(pn)} | {ds:,} | {nsz:,} | {cb:,} | "
                     f"{cm - cb:+,} | {cn - cb:+,} | {nc} |")
        T["nv"] += len(pm); T["nl"] += len(pn); T["ds"] += ds; T["ns"] += nsz
        T["cb"] += cb; T["dm"] += cm - cb; T["dn"] += cn - cb; T["nc"] += nc
        print(f"  {pop} done", file=sys.stderr)

    lines.append(f"| **TOTAL ({len(pops)} pop)** | **{T['nv']}** | **{T['nl']}** | {T['ds']:,} | "
                 f"{T['ns']:,} | {T['cb']:,} | **{T['dm']:+,}** | **{T['dn']:+,}** | **{T['nc']}** |")

    note = [
        "",
        "`novel` / `null` = unique novel peptides at 1% CLASS-SPECIFIC FDR (vs REV_nuORF| decoys only).",
        "`canon base` = canonical peptides from the canonical-ONLY DB at 1% GLOBAL FDR; `dPC` = change",
        "vs that baseline. Both arms' dPC are shown so a database-SIZE effect can be told apart from a",
        "database-QUALITY effect: canonical IDs can be lost purely because the search space grew.",
        "`ncStart` = model peptides passing FDR that appear NOWHERE in the null DB's amino-acid space",
        "(substring scan). Expected mechanism is a non-canonical start codon extending the ORF past the",
        "first in-frame AUG.",
        "",
        "**ncStart is structurally 0 until enumeration supports non-AUG starts.**",
        "`enumerate_score_orfs.candidate_orfs()` matches `is_atg = (c0==A)&(c1==T)&(c2==G)` only",
        "(verified: 400,000/400,000 BMDM candidates begin ATG), and db_model is a strict SUBSET of",
        "db_null, so no model peptide can be missing from the null space. The column is COMPUTED, not",
        "assumed, so it goes live automatically when non-AUG enumeration lands.",
    ]
    md = "\n".join(lines + note) + "\n"
    print("\n" + md)
    if a.out:
        Path(a.out).write_text(md)
        print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
