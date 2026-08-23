#!/usr/bin/env python3
"""Search-database construction (`pgx` step 6).

Every database shares ONE canonical base (the GENCODE proteome) plus inline reversed `REV_` decoys,
so the arms differ ONLY in their novel-ORF content and the comparison isolates SELECTION:

  db_gencode        GENCODE proteome only                        BASELINE detection
  db_model_<arm>    GENCODE + significant RiboCode calls on the model's predicted signal
                    + the N-terminal extension entries (a DISTINCT accession class)
  db_null_atg       GENCODE + every ATG..stop ORF on the expressed universe
  db_null_nc        GENCODE + every ATG-or-near-cognate ORF on the same universe

Both nulls use the SAME expressed universe as the model arm (TPM >= 1, protein_coding + lncRNA,
no chrM, length cap), so DB size differences reflect the selection rule and nothing else.
`db_null_atg` is a strict subset of `db_null_nc` by construction, which makes the pair directly
interpretable: the increment is exactly the cost of naive non-AUG enumeration.

NO DOUBLE COUNTING AGAINST GENCODE. A novel sequence identical to any GENCODE protein is dropped
here at build time. That is necessary but NOT sufficient, because a peptide is a FRAGMENT of a
protein: a novel ORF can be non-identical to every GENCODE protein yet still yield only peptides
that occur inside one. The complementary substring check runs at report time (pgx.report), which is
the only place peptides exist.

`internal` ORFs (CDS truncations and out-of-frame reads of a CDS) are excluded from every database,
matching build_a549_dbs / build_ribocode_db, so the classes stay comparable across DB builders.

  python -m pgx.build_dbs --species mouse --universe-fa u.fa --out <dir> \
      --calls <calls>/calls_poisson.tsv --extensions <ext>/extensions.tsv --arm poisson
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

from . import rc_io
from .refs import species_refs
from .seqtools import (candidate_orfs, cds_relationship, iter_fasta, orf_class, parse_starts,
                       translate)

NOVEL_CLASSES = {"uORF", "dORF", "lncRNA_orf", "nterm_ext"}
CLASS_MAP_COLS = ("accession", "class", "start_codon", "source", "score")


def load_canonical(path, min_aa):
    """{sequence: accession} for the GENCODE proteome, deduplicated by sequence."""
    canon = {}
    for name, seq in iter_fasta(path):
        seq = seq.replace("*", "").strip()
        if len(seq) >= min_aa and seq not in canon:
            canon[seq] = name.split("|")[0].split()[0]
    return canon


def _add(store, seq, acc, klass, codon, source, score, canon_seqs, min_aa):
    """Insert a novel protein unless too short or identical to a canonical protein."""
    if len(seq) < min_aa or seq in canon_seqs or "*" in seq:
        return False
    prev = store.get(seq)
    if prev is None or score > prev[3]:
        store[seq] = (acc, klass, codon, score, source)
    return True


def novel_from_calls(calls_tsv, canon_seqs, min_aa):
    """Model arm: significant RiboCode ORF calls."""
    store = {}
    with open(calls_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq = (r.get("aaseq") or "").replace("*", "").strip()
            acc = f"nuORF|{r['pgx_class']}|{r['tx']}_{r['start0']}_{r['end0']}"
            try:
                score = -float(r.get("qval", 1) or 1)
            except ValueError:
                score = -1.0
            _add(store, seq, acc, r["pgx_class"], r.get("start_codon", "ATG"),
                 "ribocode", score, canon_seqs, min_aa)
    return store


def novel_from_collapsed(paths, canon_seqs, min_aa, biotype, qcut=0.05):
    """Novel ORFs straight from RiboCode `*_collapsed.txt`, i.e. calls on MEASURED Ribo-seq.

    This is the experimental counterpart to novel_from_calls: same class mapping, same dedup, same
    GENCODE exclusion, but the evidence is a real ribosome profiling experiment rather than a
    predicted profile. Several files are pooled (e.g. an untreated and a stimulated condition of
    the same cell type), since a union over conditions is what "what Ribo-seq says is translated
    here" actually means.

    Extensions are NOT harvested here. RiboCode reports them inside its `annotated` rows and does
    not test whether the upstream start is used (see pgx.extensions), so including them would add
    unsupported entries; the model arms get their extensions from a test RiboCode cannot perform.
    Note the asymmetry when comparing arms.
    """
    store = {}
    for p in paths:
        for o in rc_io.read_collapsed(p):
            if o["qval"] > qcut:
                continue
            klass = rc_io.pgx_class(o, biotype.get(o["tx"]))
            if not klass or klass not in NOVEL_CLASSES:
                continue
            seq = (o["aaseq"] or "").replace("*", "").strip()
            acc = f"nuORF|{klass}|{o['tx']}_{o['start0']}_{o['end0']}"
            _add(store, seq, acc, klass, o.get("start_codon", "ATG"), "ribocode_measured",
                 -float(o["qval"]), canon_seqs, min_aa)
    return store


def novel_from_extensions(ext_tsv, canon_seqs, min_aa, qcut=None):
    """Model arm: N-terminal CDS extensions, kept as their OWN class/namespace (`nterm_ext`)."""
    store = {}
    with open(ext_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if qcut is not None and float(r.get("qval", 1)) > qcut:
                continue
            seq = (r.get("aaseq") or "").replace("*", "").strip()
            acc = f"nuORF|nterm_ext|{r['tx']}_{r['ext_start0']}_{r['end0']}"
            _add(store, seq, acc, "nterm_ext", r.get("start_codon", "ATG"),
                 "extension", -float(r.get("qval", 1) or 1), canon_seqs, min_aa)
    return store


def novel_enumerated(universe_fa, starts, canon_seqs, min_aa, tx2cds, biotype):
    """Null arms: every <start>..stop ORF on the expressed universe, naively enumerated."""
    store = {}
    min_nt = (min_aa + 1) * 3
    n_scanned = 0
    for tx, seq in iter_fasta(universe_fa):
        n_scanned += 1
        cds = tx2cds.get(tx)
        bt = biotype.get(tx, "unknown")
        for a, e, codon in candidate_orfs(seq, min_nt, starts):
            if e > len(seq):
                continue
            klass = orf_class(cds_relationship(a, e, cds), bt)
            if klass not in NOVEL_CLASSES:
                continue
            prot = translate(seq[a:e - 3].upper().replace("U", "T"))
            _add(store, prot, f"nuORF|{klass}|{tx}_{a}_{e}", klass, codon,
                 "enumerated", 0.0, canon_seqs, min_aa)
    return store, n_scanned


def write_db(path, canon, novel, class_map_path=None):
    """Canonical base + novel set, each followed by its inline reversed REV_ decoy."""
    with open(path, "w") as fh:
        for seq, acc in canon.items():
            fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
        for seq, (acc, klass, codon, score, source) in novel.items():
            fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
    if class_map_path:
        with open(class_map_path, "w") as fh:
            fh.write("\t".join(CLASS_MAP_COLS) + "\n")
            for _seq, (acc, klass, codon, score, source) in novel.items():
                fh.write(f"{acc}\t{klass}\t{codon}\t{source}\t{score:.6g}\n")
    return len(canon) + len(novel)


def summarize(novel):
    return {"n": len(novel),
            "by_class": dict(Counter(v[1] for v in novel.values()).most_common()),
            "by_start": dict(Counter(v[2] for v in novel.values()).most_common()),
            "by_source": dict(Counter(v[4] for v in novel.values()).most_common())}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--universe-fa", default=None, help="expressed universe FASTA (needed for nulls)")
    ap.add_argument("--calls", default=None, help="calls_<arm>.tsv from pgx.call_orfs")
    ap.add_argument("--ribocode-collapsed", nargs="*", default=None,
                    help="RiboCode *_collapsed.txt from MEASURED Ribo-seq; several are pooled. "
                         "Builds db_model_<arm> from real experimental calls instead of predicted "
                         "ones, so an experimental arm can be compared against the model arms on "
                         "the same spectra. Mutually exclusive with --calls.")
    ap.add_argument("--orf-qvalue", type=float, default=0.05,
                    help="q cutoff applied to --ribocode-collapsed calls")
    ap.add_argument("--extensions", default=None, help="extensions.tsv from pgx.extensions")
    ap.add_argument("--arm", default="poisson", help="label for the model DB (db_model_<arm>)")
    ap.add_argument("--dbs", default="gencode,model,null_atg,null_nc",
                    help="which databases to build")
    # STANDING RULE: the ORF-length floor is set BY ASSAY, not by pipeline default.
    #   tryptic whole-cell lysate -> 30 aa (the same floor the ORF-call track applies)
    #   MHC / HLA immunopeptidomics -> 7 aa (HLA-I peptides are 8-11 aa)
    # Leaving 7 aa on a tryptic search imports a large sub-30-aa population that the ORF-call track
    # excludes by rule. Measured cost of the floor (orf_length_audit.py): 0.0% of model discoveries
    # on mouse macrophage tryptic data, but 26.9-32.9% on human HLA-I -- so the two assays genuinely
    # need different floors, and neither number transfers to the other.
    ap.add_argument("--assay", choices=["tryptic", "mhc"], default=None,
                    help="assay type; SETS --min-aa (tryptic=30, mhc=7) and is the preferred way to "
                         "specify it. Passing --min-aa as well is allowed only if it agrees.")
    ap.add_argument("--min-aa", type=int, default=None,
                    help="minimum protein length. Prefer --assay. Defaults to 7 with a warning "
                         "when neither is given, for backward compatibility.")
    ap.add_argument("--null-starts", default="ATG", help="start codons for db_null_atg")
    ap.add_argument("--nc-starts", default="near_cognate", help="start codons for db_null_nc")
    ap.add_argument("--ext-qvalue", type=float, default=None,
                    help="extra q filter on extensions (they are already filtered by pgx.extensions)")
    ap.add_argument("--canonical", default=None, help="override the GENCODE proteome")
    ap.add_argument("--tx2cds", default=None)
    ap.add_argument("--tx2biotype", default=None)
    a = ap.parse_args()

    refs = species_refs(a.species, gencode_proteome=a.canonical, tx2cds=a.tx2cds,
                        tx2biotype=a.tx2biotype)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    want = {x.strip() for x in a.dbs.split(",") if x.strip()}

    # Resolve the floor from the assay, and refuse a contradiction rather than silently picking one.
    ASSAY_MIN_AA = {"tryptic": 30, "mhc": 7}
    if a.assay:
        # NOT `want` -- that name is already the SET of databases to build (line ~213), and
        # shadowing it with an int made `if "gencode" in want` raise
        # "TypeError: argument of type 'int' is not iterable" AFTER the assay had resolved
        # correctly, so the log showed a working floor and a crash together.
        want_min_aa = ASSAY_MIN_AA[a.assay]
        if a.min_aa is not None and a.min_aa != want_min_aa:
            sys.exit(f"ABORT: --assay {a.assay} requires --min-aa {want_min_aa}, got {a.min_aa}. "
                     "The floor is set by the assay (see docs/PIPELINE_POLICY.md).")
        a.min_aa = want_min_aa
        print(f"assay={a.assay} -> min_aa={a.min_aa}")
    elif a.min_aa is None:
        a.min_aa = 7
        print("WARNING: neither --assay nor --min-aa given; defaulting to min_aa=7 (MHC). "
              "Tryptic whole-cell-lysate databases MUST use --assay tryptic (30 aa).",
              file=sys.stderr)
    else:
        print(f"WARNING: --min-aa {a.min_aa} given without --assay; the floor should be set by "
              "assay (tryptic=30, mhc=7).", file=sys.stderr)

    canon = load_canonical(refs["gencode_proteome"], a.min_aa)
    canon_seqs = set(canon)
    print(f"GENCODE proteins (unique, >= {a.min_aa} aa): {len(canon):,}")
    # Merge into any existing summary: run.py invokes this once for the model-independent
    # databases and once per calling arm, so overwriting would drop the earlier entries.
    summary_path = out / "db_summary.json"
    report = {"species": a.species, "min_aa": a.min_aa, "canonical": len(canon), "dbs": {}}
    if summary_path.exists():
        try:
            prev = json.loads(summary_path.read_text())
            report["dbs"].update(prev.get("dbs", {}))
        except (ValueError, OSError):
            pass

    if "gencode" in want:
        n = write_db(out / "db_gencode.fasta", canon, {})
        print(f"db_gencode        : {n:,} target + {n:,} decoy")
        report["dbs"]["gencode"] = {"targets": n, "novel": 0}

    if "model" in want:
        if bool(a.calls) == bool(a.ribocode_collapsed):
            raise SystemExit("pass exactly one of --calls or --ribocode-collapsed")
        if a.ribocode_collapsed:
            bt = rc_io.load_biotype(refs["tx2biotype"])
            novel = novel_from_collapsed(a.ribocode_collapsed, canon_seqs, a.min_aa, bt,
                                         a.orf_qvalue)
        else:
            novel = novel_from_calls(a.calls, canon_seqs, a.min_aa)
        n_calls = len(novel)
        n_ext = 0
        if a.extensions and Path(a.extensions).exists():
            ext = novel_from_extensions(a.extensions, canon_seqs, a.min_aa, a.ext_qvalue)
            before = len(novel)
            for seq, v in ext.items():
                if seq not in novel:
                    novel[seq] = v
            n_ext = len(novel) - before
        tag = f"db_model_{a.arm}"
        n = write_db(out / f"{tag}.fasta", canon, novel, out / f"{tag}_class_map.tsv")
        s = summarize(novel)
        print(f"{tag:<18}: {n:,} target + {n:,} decoy   novel={s['n']:,} "
              f"(calls {n_calls:,} + extensions {n_ext:,})")
        print(f"   by class: {s['by_class']}")
        print(f"   by start: {s['by_start']}")
        report["dbs"][f"model_{a.arm}"] = {"targets": n, **s,
                                           "from_calls": n_calls, "from_extensions": n_ext}

    if want & {"null_atg", "null_nc"}:
        if not a.universe_fa:
            raise SystemExit("--universe-fa is required to build the null DBs")
        tx2cds = rc_io.load_tx2cds(refs["tx2cds"])
        biotype = rc_io.load_biotype(refs["tx2biotype"])
        for name, spec in (("null_atg", a.null_starts), ("null_nc", a.nc_starts)):
            if name not in want:
                continue
            starts = parse_starts(spec)
            novel, n_tx = novel_enumerated(a.universe_fa, starts, canon_seqs, a.min_aa,
                                           tx2cds, biotype)
            n = write_db(out / f"db_{name}.fasta", canon, novel, out / f"db_{name}_class_map.tsv")
            s = summarize(novel)
            print(f"db_{name:<14}: {n:,} target + {n:,} decoy   novel={s['n']:,} "
                  f"over {n_tx:,} tx, starts={','.join(starts)}")
            print(f"   by class: {s['by_class']}")
            report["dbs"][name] = {"targets": n, "universe_tx": n_tx,
                                   "starts": list(starts), **s}

    summary_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {out}/db_*.fasta + class maps + db_summary.json")


if __name__ == "__main__":
    main()
