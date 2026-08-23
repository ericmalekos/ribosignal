#!/usr/bin/env python3
"""Derive a min_aa=N search database from an existing one by length filtering. Exactly equivalent.

STANDING RULE (docs/PIPELINE_POLICY.md): the ORF-length floor is set BY ASSAY -- 30 aa for tryptic
whole-cell lysate, 7 aa for MHC/HLA immunopeptidomics. The macrophage databases were built at 7 aa
and are tryptic, so they need 30.

WHY FILTERING IS EQUIVALENT TO REBUILDING, AND WHY THAT MATTERS HERE. In `pgx.build_dbs`, `min_aa`
is a PURE LENGTH FILTER on the final protein sequences and nothing else:

    load_canonical: if len(seq) >= min_aa and seq not in canon
    _add:           if len(seq) < min_aa or seq in canon_seqs or "*" in seq: return

It does not change how ORFs are called, translated, deduplicated, or classified. So
`build_dbs --min-aa 30` produces exactly the >= 30 aa subset of `build_dbs --min-aa 7`. Filtering the
existing FASTA therefore reproduces the rebuild bit-for-bit while REUSING the original call sets --
which matters because `db_percall/` carries no `db_summary.json`, so a from-scratch rebuild would
risk silently substituting a different call set and confounding the floor change with a provenance
change.

One subtlety that also comes out identical: raising the floor shrinks the CANONICAL set too, so a
novel sequence identical to a sub-30-aa canonical protein would no longer be excluded by
`canon_seqs`. But such a sequence is itself sub-30-aa and is dropped by the length test regardless.

DECOYS FILTER CORRECTLY BY CONSTRUCTION. A `REV_` entry is its target reversed and therefore the same
length, so one length test keeps target and decoy in step. The script asserts the target and decoy
counts still match afterwards -- an imbalance would silently corrupt FDR.

  python3 -m pgx.filter_db_min_aa --in <db_dir> --out <db_dir> --min-aa 30
"""
from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import sys

NOVEL_PREFIX = "nuORF|"
DECOY_PREFIX = "REV_"


def read_fasta(p):
    name, seq = None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def filter_one(src, dst, min_aa):
    kept = dropped = 0
    n_t = n_d = 0
    acc_kept = set()
    with open(dst, "w") as out:
        for name, seq in read_fasta(src):
            if len(seq) < min_aa:
                dropped += 1
                continue
            kept += 1
            acc = name.split()[0]
            acc_kept.add(acc)
            if acc.startswith(DECOY_PREFIX):
                n_d += 1
            else:
                n_t += 1
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    # A target/decoy imbalance would silently corrupt every FDR computed downstream.
    assert n_t == n_d, f"{src.name}: target/decoy imbalance after filtering ({n_t} vs {n_d})"
    novel_t = sum(1 for a in acc_kept if a.startswith(NOVEL_PREFIX))
    return dict(kept=kept, dropped=dropped, targets=n_t, decoys=n_d, novel_targets=novel_t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, help="existing database directory")
    ap.add_argument("--out", required=True, help="output directory (must differ from --in)")
    ap.add_argument("--min-aa", type=int, required=True)
    ap.add_argument("--assay", choices=["tryptic", "mhc"], default=None,
                    help="if given, must agree with --min-aa (tryptic=30, mhc=7)")
    a = ap.parse_args()
    ASSAY_MIN_AA = {"tryptic": 30, "mhc": 7}
    if a.assay and ASSAY_MIN_AA[a.assay] != a.min_aa:
        sys.exit(f"ABORT: --assay {a.assay} requires --min-aa {ASSAY_MIN_AA[a.assay]}, got {a.min_aa}")

    src, dst = pathlib.Path(a.src), pathlib.Path(a.out)
    if src.resolve() == dst.resolve():
        sys.exit("ABORT: --out must differ from --in (this never edits in place)")
    dst.mkdir(parents=True, exist_ok=True)

    summary = {"derived_from": str(src), "min_aa": a.min_aa, "assay": a.assay,
               "method": "length filter of an existing min_aa=7 database; equivalent to "
                         "build_dbs --min-aa because min_aa is a pure length filter",
               "dbs": {}}
    fastas = sorted(src.glob("*.fasta"))
    if not fastas:
        sys.exit(f"ABORT: no *.fasta in {src}")
    for fa in fastas:
        r = filter_one(fa, dst / fa.name, a.min_aa)
        summary["dbs"][fa.stem] = r
        print(f"  {fa.name:<34} kept {r['kept']:>8,}  dropped {r['dropped']:>8,}  "
              f"novel targets {r['novel_targets']:>7,}")
    # Class maps are accession-keyed; carry them across, restricted to surviving accessions.
    for cm in sorted(src.glob("*_class_map.tsv")):
        fa = src / (cm.name.replace("_class_map.tsv", ".fasta"))
        keep = {n.split()[0] for n, s in read_fasta(fa) if len(s) >= a.min_aa} if fa.exists() else None
        lines = cm.read_text().splitlines()
        with open(dst / cm.name, "w") as out:
            out.write(lines[0] + "\n")
            n = 0
            for ln in lines[1:]:
                if keep is None or ln.split("\t")[0] in keep:
                    out.write(ln + "\n")
                    n += 1
        print(f"  {cm.name:<34} class-map rows kept {n:,}")
    (dst / "db_summary_filtered.json").write_text(json.dumps(summary, indent=2))
    print(f"\n  wrote {dst}")


if __name__ == "__main__":
    main()
