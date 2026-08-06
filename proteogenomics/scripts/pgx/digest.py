#!/usr/bin/env python3
"""Reduce a raw MSFragger search directory to the compact digest the report actually needs.

The search writes ~1-5 GB of .tsv and .pin per arm, but everything downstream (pgx.report,
pgx.replication, pgx.score_compare) reads only rank-1 rows and only four fields from them. Keeping
the raw output on the shared filesystem costs ~100x more space than the information is worth, and
the group quota is the binding constraint on this project.

So this runs ON THE COMPUTE NODE, against node-local scratch, and emits one gzipped file. Only that
file crosses to the group filesystem. `/data/tmp` is node-local, so the reduction MUST happen in
the same job that ran the search; a later job on another node cannot see the scratch.

Columns kept: spec_id ('<fraction>:<scannum>'), peptide, hyperscore, proteins.
That is sufficient for rank-1 FDR, class-specific FDR, per-fraction replication, and per-spectrum
score matching.

  python -m pgx.digest --search-dir <scratch> --out <group-fs>/rank1.tsv.gz
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import sys
from pathlib import Path

COLS = ("spec_id", "peptide", "hyperscore", "proteins")


def _msfragger(fh, stem, w):
    n_in = n_out = 0
    for r in csv.DictReader(fh, delimiter="\t"):
        n_in += 1
        if r.get("hit_rank") != "1":
            continue
        n_out += 1
        w.write(f"{stem}:{r.get('scannum', '')}\t{r['peptide']}\t"
                f"{r['hyperscore']}\t{r.get('proteins', '')}\n")
    return n_in, n_out


def _comet(fh, stem, w):
    """Comet .txt: a 1-line header banner, then a real header row.

    Score column is `xcorr` rather than `hyperscore`; both are 'higher is better', which is all the
    FDR code assumes. Comet lists additional matching proteins in `protein` as a comma/space list
    with the count in `duplicate_protein_count`, so the accession field is normalised to the same
    `;`-delimited form the MSFragger path produces.
    """
    first = fh.readline()
    if not first.lower().startswith("scan"):        # skip the "CometVersion ..." banner line
        pass
    else:
        fh.seek(0)
    n_in = n_out = 0
    for r in csv.DictReader(fh, delimiter="\t"):
        n_in += 1
        if str(r.get("num", r.get("rank", "1"))).strip() != "1":
            continue
        n_out += 1
        prot = (r.get("protein") or "").replace(", ", ";").replace(",", ";").strip()
        w.write(f"{stem}:{r.get('scan', '')}\t{r['plain_peptide'] if 'plain_peptide' in r else r['peptide']}\t"
                f"{r.get('xcorr', '0')}\t{prot}\n")
    return n_in, n_out


def digest(search_dir, out_path, engine="msfragger"):
    pat = "*.txt" if engine == "comet" else "*.tsv"
    files = sorted(glob.glob(str(Path(search_dir) / pat)))
    if not files:
        raise SystemExit(f"no {pat} under {search_dir}")
    parse = _comet if engine == "comet" else _msfragger
    n_in = n_out = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt", compresslevel=6) as w:
        w.write("\t".join(COLS) + "\n")
        for t in files:
            stem = Path(t).stem
            with open(t) as fh:
                a, b = parse(fh, stem, w)
                n_in += a
                n_out += b
    return len(files), n_in, n_out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search-dir", required=True, help="directory of MSFragger *.tsv (scratch)")
    ap.add_argument("--out", required=True, help="gzipped digest to write (group filesystem)")
    ap.add_argument("--engine", default="msfragger", choices=["msfragger", "comet"],
                    help="which search engine wrote the directory")
    a = ap.parse_args()
    n_f, n_in, n_out = digest(a.search_dir, a.out, a.engine)
    sz = Path(a.out).stat().st_size
    print(f"digest: {n_f} fractions, {n_in:,} rows -> {n_out:,} rank-1 -> "
          f"{sz / 1e6:.1f} MB gz  ({a.out})", file=sys.stderr)


if __name__ == "__main__":
    main()
