#!/usr/bin/env python3
"""FDR control and the standard output table (`pgx` step 8).

WHICH FDR WHERE (project standing rule; never mix filtered and unfiltered columns):
  canonical / dPSM : GLOBAL 1% FDR (all targets vs all decoys). This is what the search reports and
                     it keeps the target-decoy competition intact. A class-specific FDR on the
                     canonical side is confounded, because a large novel space cannibalises
                     canonical DECOYS faster than canonical targets, deflating the estimate until a
                     junk database appears to GAIN canonical identifications.
  novel            : CLASS-SPECIFIC 1% FDR (novel targets vs REV_nuORF| decoys ONLY). A global FDR
                     inflates non-canonical discovery roughly 10-13x.
Raw rank-1 counts are never reported: they measure spectra ABSORBED, which scales with DB size.

NO DOUBLE COUNTING AGAINST GENCODE. MSFragger reports every protein containing a peptide, so a
peptide that occurs in ANY GENCODE protein is classified canonical even when it also occurs in a
novel ORF. Those are real: 3,489 of 191,358 rank-1 PSMs in the BMDM model search map to both. The
protein list is `;`-delimited (NOT `,`); splitting on the wrong character silently merges the list
into one token, and the classification then survives only because MSFragger happens to list the
canonical protein first. This module splits on both so correctness does not depend on that order.

ncStart is the model-only discovery column: peptides passing FDR that appear NOWHERE in a null
database's amino-acid space (substring scan, since a peptide is a fragment of a protein). Reported
against BOTH nulls -- vs null_atg answers "would an AUG-only pipeline ever have found this", vs
null_nc answers "would even a naive near-cognate enumeration have found this".

  python -m pgx.report --search-root <dir> --db-dir <dir> --out table.md
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import json
import sys
from pathlib import Path

from .seqtools import iter_fasta

FDR = 0.01
NOVEL_PREFIX = "nuORF|"
DECOY_PREFIX = "REV_"


def split_accs(field):
    """MSFragger delimits the protein list with ';'. Tolerate ',' too."""
    return [a for a in field.replace(",", ";").split(";") if a]


def klass(accs):
    """canon_t / novel_t / novel_d / other_d. A peptide in ANY GENCODE protein is canonical."""
    tgt = [a for a in accs if not a.startswith(DECOY_PREFIX)]
    if tgt:
        return "novel_t" if all(a.startswith(NOVEL_PREFIX) for a in tgt) else "canon_t"
    return "novel_d" if all(a.startswith(DECOY_PREFIX + NOVEL_PREFIX) for a in accs) else "other_d"


def _rows(d):
    """Yield dict rows from the compact rank1.tsv.gz digest if present, else the raw *.tsv.

    pgx.digest reduces a search directory to rank-1 rows on the compute node, so only ~1% of the
    bytes reach the shared filesystem. Preferring it here means reports work identically whether a
    search kept its raw output or only the digest.
    """
    dz = Path(d) / "rank1.tsv.gz"
    if dz.exists():
        with gzip.open(dz, "rt") as fh:
            yield from csv.DictReader(fh, delimiter="\t")
        return
    for t in sorted(glob.glob(str(Path(d) / "*.tsv"))):
        with open(t) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("hit_rank") == "1":
                    yield r


def load_rank1(d):
    """[(peptide, hyperscore, klass, novel_accessions)] over every fraction in a search dir."""
    out = []
    for r in _rows(d):
        accs = split_accs(r.get("proteins", ""))
        k = klass(accs)
        nov = tuple(a for a in accs if a.startswith(NOVEL_PREFIX)) if k == "novel_t" else ()
        out.append((r["peptide"], float(r["hyperscore"]), k, nov))
    return out


def load_rank1_spec(d):
    """[(spec_id, peptide, hyperscore, klass)] -- spec_id = '<fraction stem>:<scannum>'.

    Keeps spectrum identity, which is what makes a genuinely paired comparison between two search
    arms possible: the SAME spectrum matched to the SAME peptide in both. Comparing
    max-score-per-peptide instead is biased by how many spectra each arm passes.
    """
    out = []
    for r in _rows(d):
        sid = r.get("spec_id")
        if sid is None:                       # raw *.tsv fallback has no spec_id column
            sid = f"{r.get('_stem', '')}:{r.get('scannum', '')}"
        out.append((sid, r["peptide"], float(r["hyperscore"]),
                    klass(split_accs(r.get("proteins", "")))))
    return out


def cut(rows, tset, dset, fdr=FDR):
    """Highest score threshold whose running decoy/target ratio is still <= fdr."""
    ranked = sorted(((h, k) for (_, h, k, _) in rows if k in tset | dset), key=lambda x: -x[0])
    t = d = 0
    last = None
    for h, k in ranked:
        if k in dset:
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= fdr:
            last = h
    return last


def canonical_at_global_fdr(rows, fdr=FDR):
    g = cut(rows, {"canon_t", "novel_t"}, {"novel_d", "other_d"}, fdr)
    if g is None:
        return 0, set()
    psms = [(p, h) for (p, h, k, _) in rows if k == "canon_t" and h >= g]
    return len(psms), {p for p, _ in psms}


def novel_at_class_fdr(rows, fdr=FDR):
    """-> (n_psms, peptides, accessions, peptide -> supporting novel accessions)."""
    c = cut(rows, {"novel_t"}, {"novel_d"}, fdr)
    if c is None:
        return 0, set(), set(), {}
    hit = [(p, a) for (p, h, k, a) in rows if k == "novel_t" and h >= c]
    peps = {p for p, _ in hit}
    accs = {x for _, a in hit for x in a}
    by_pep = {}
    for p, a in hit:
        by_pep.setdefault(p, set()).update(a)
    return len(hit), peps, accs, by_pep


def aa_space(fasta):
    """One string of every TARGET sequence, sentinel-joined, for substring containment."""
    return "|".join(s for n, s in iter_fasta(fasta) if not n.startswith(DECOY_PREFIX))


def db_novel_count(class_map):
    p = Path(class_map)
    return sum(1 for _ in open(p)) - 1 if p.exists() else 0


def load_class_map(class_map):
    """{accession: {class, start_codon, source}} for an arm's novel entries."""
    p = Path(class_map)
    if not p.exists():
        return {}
    with open(p) as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def nonaug_accessions(cmap):
    return {a for a, r in cmap.items() if (r.get("start_codon") or "ATG").upper() != "ATG"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search-root", required=True, help="dir containing <arm>/*.tsv")
    ap.add_argument("--db-dir", required=True, help="dir containing db_<arm>.fasta + class maps")
    ap.add_argument("--arms", default=None,
                    help="comma list; default = every subdirectory of --search-root")
    ap.add_argument("--baseline", default="gencode", help="arm used as the canonical baseline")
    ap.add_argument("--label", default="", help="row label (population / dataset)")
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--out", default=None, help="write markdown here as well as stdout")
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args()

    root, dbd = Path(a.search_root), Path(a.db_dir)
    arms = ([x.strip() for x in a.arms.split(",") if x.strip()] if a.arms
            else sorted(p.name for p in root.iterdir() if p.is_dir()))
    if a.baseline not in arms:
        raise SystemExit(f"baseline arm {a.baseline!r} not among {arms}")

    loaded = {}
    empty = []
    for arm in arms:
        rows = load_rank1(root / arm)
        if not rows:
            empty.append(arm)
            continue
        loaded[arm] = rows
        print(f"  loaded {arm}: {len(rows):,} rank-1 PSMs", file=sys.stderr)
    # Fail loudly rather than emit a table that silently omits a requested arm. A partial table
    # looks complete, and a report run before its searches finished has already produced one:
    # `model_predicted` was dropped from all 12 macrophage reports because the report jobs started
    # while those searches were still writing.
    if empty and not a.allow_missing:
        raise SystemExit(f"no search output for requested arm(s): {', '.join(empty)}\n"
                         f"  (searches still running? pass --allow-missing to report anyway)")
    if empty:
        print(f"[warn] omitting empty arm(s): {', '.join(empty)}", file=sys.stderr)
    if a.baseline not in loaded:
        raise SystemExit(f"baseline arm {a.baseline!r} has no search output")

    base_psm, base_pep = canonical_at_global_fdr(loaded[a.baseline], a.fdr)

    # null AA spaces, built once and reused for every arm's ncStart columns
    spaces = {}
    for null in ("null_atg", "null_nc"):
        fa = dbd / f"db_{null}.fasta"
        if fa.exists():
            spaces[null] = aa_space(fa)
            print(f"  {null} AA space: {len(spaces[null]) / 1e6:.1f} Mb", file=sys.stderr)

    res = {}
    for arm in arms:
        if arm not in loaded:
            continue
        rows = loaded[arm]
        c_psm, c_pep = canonical_at_global_fdr(rows, a.fdr)
        n_psm, n_pep, n_acc, by_pep = novel_at_class_fdr(rows, a.fdr)
        cm = dbd / f"db_{arm}_class_map.tsv"
        cmap = load_class_map(cm)
        nonaug = nonaug_accessions(cmap)
        r = {"novel_psms": n_psm, "novel_peptides": len(n_pep), "novel_seqs_with_psm": len(n_acc),
             "db_novel_seqs": db_novel_count(cm),
             "gencode_psms": c_psm, "gencode_peptides": len(c_pep),
             "d_gencode_psms": c_psm - base_psm, "d_gencode_peptides": len(c_pep) - len(base_pep)}
        for null, sp in spaces.items():
            if arm == null or not n_pep:
                r[f"absent_from_{null}"] = 0
                r[f"ncstart_vs_{null}"] = 0
                continue
            absent = {p for p in n_pep if p not in sp}
            r[f"absent_from_{null}"] = len(absent)
            # ncStart is absence attributable to the START CODON, which is what the column claims.
            # A peptide can also be absent because the two arms classify an ORF differently (the
            # null enumerator and RiboCode must agree, see seqtools.cds_relationship), so requiring
            # a non-AUG supporting entry makes the column immune to that and to any future
            # classifier drift.
            r[f"ncstart_vs_{null}"] = sum(1 for p in absent if by_pep.get(p, set()) & nonaug)
        res[arm] = r

    # Header wording is load-bearing: PSM columns are TOTAL rank-1 spectra (a peptide seen in five
    # spectra counts five times), peptide columns are DISTINCT sequences. Abbreviating either to
    # "pept" invites reading a spectrum count as a peptide count.
    cols = [("arm", "arm"),
            ("novel PSMs (total)", "novel_psms"),
            ("novel peptides (uniq)", "novel_peptides"),
            ("novel seqs w/PSM", "novel_seqs_with_psm"),
            ("DB novel seqs", "db_novel_seqs"),
            ("GENCODE PSMs (total)", "gencode_psms"),
            ("GENCODE peptides (uniq)", "gencode_peptides"),
            ("dPSM (total)", "d_gencode_psms"),
            ("dPeptide (uniq)", "d_gencode_peptides")]
    for null in spaces:
        short = null.replace("null_", "")
        cols.append((f"absent vs {short}", f"absent_from_{null}"))
        cols.append((f"ncStart vs {short}", f"ncstart_vs_{null}"))

    head = "| " + " | ".join(c[0] for c in cols) + " |"
    sep = "|---" + "|--:" * (len(cols) - 1) + "|"
    lines = [f"### {a.label}" if a.label else "", head, sep]
    for arm in arms:
        if arm not in res:
            continue
        cells = [arm]
        for _, key in cols[1:]:
            v = res[arm][key]
            cells.append(f"{v:+,}" if key.startswith("d_") else f"{v:,}")
        lines.append("| " + " | ".join(cells) + " |")

    lines += [
        "",
        "**PSM columns are TOTAL rank-1 spectra; peptide columns are DISTINCT sequences.** A peptide "
        "identified in five spectra contributes 5 to a PSM column and 1 to a peptide column, so the "
        "two never agree and the ratio (spectra per peptide) is itself informative -- identifications "
        "carried by a single spectrum are the weakest.",
        "`novel` columns: 1% CLASS-SPECIFIC FDR (novel targets vs `REV_nuORF|` decoys only).",
        f"`GENCODE` columns and `dPSM`/`dPept`: 1% GLOBAL FDR, relative to the `{a.baseline}` "
        "baseline, so a database-SIZE cost is separable from a database-QUALITY cost.",
        "A peptide occurring in ANY GENCODE protein is counted as GENCODE, never as novel.",
        "`absent vs <null>` = passing novel peptides found nowhere in that null's amino-acid space "
        "(substring scan). `ncStart` is the subset of those whose supporting entry opens at a "
        "NON-AUG codon, i.e. absence attributable to the start codon rather than to any "
        "classification difference between the arms. vs `atg` is the AUG-only naive pipeline, "
        "vs `nc` the near-cognate one.",
    ]
    md = "\n".join(x for x in lines if x is not None) + "\n"
    print("\n" + md)
    if a.out:
        Path(a.out).write_text(md)
        print(f"wrote {a.out}", file=sys.stderr)
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(
            {"label": a.label, "baseline": a.baseline, "fdr": a.fdr,
             "baseline_gencode_psms": base_psm, "arms": res}, indent=2) + "\n")
        print(f"wrote {a.json_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
