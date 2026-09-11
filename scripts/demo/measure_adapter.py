#!/usr/bin/env python3
"""Report which 3' adapter, if any, is actually present in a FASTQ.

The pipeline never defaults an adapter, and this is why: the public RNA-seq for this
demo (SRR15513269) turns out to be adapter-trimmed already -- Illumina TruSeq appears in
4 reads out of 400,000 -- so passing `-a AGATCGGAAGAGC` would be trimming something that
is not there, and `--discard-untrimmed`, which is correct for ribosome footprints, would
have thrown away essentially the whole library.

Prints the hit rate for the common adapters plus an unbiased scan of the most frequent
k-mers at read positions where read-through would appear. Decide from the numbers:

  a known adapter in a large fraction of reads   -> pass it with -a
  everything near zero, read lengths ragged      -> already trimmed, pass no -a

    python scripts/demo/measure_adapter.py reads.fastq.gz [--n 400000]
"""
from __future__ import annotations

import argparse
import collections
import gzip
import sys

KNOWN = {
    "Illumina TruSeq / universal": "AGATCGGAAGAGC",
    "Illumina small RNA 3'":       "TGGAATTCTCGGGTGCCAAGG",
    # The OLD small-RNA linker, standard in Ingolia-era Ribo-seq and still common in datasets
    # from ~2015-2020. It was missing here, and on GSE143393 that produced a confident
    # "no adapter detected -- the library is already trimmed" while 91.4% of reads carried it.
    # Only the unbiased k-mer listing caught it. Trusting the verdict would have sent reads to
    # STAR with the linker attached, which is the 0%-unique-mapping failure.
    "Illumina small RNA 3' (old linker)": "CTGTAGGCACCATCAAT",
    "Nextera":                     "CTGTCTCTTATACACATCT",
    "NEBNext small RNA":           "AGATCGGAAGAGCACACGTCT",
    "polyA":                       "AAAAAAAAAAAA",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq")
    ap.add_argument("--n", type=int, default=400_000, help="reads to sample")
    ap.add_argument("--from-pos", type=int, default=25,
                    help="unbiased k-mer scan starts here (default %(default)s)")
    a = ap.parse_args()

    op = gzip.open if a.fastq.endswith(".gz") else open
    hits = collections.Counter()
    kmers = collections.Counter()
    lens = collections.Counter()
    with op(a.fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i // 4 >= a.n:
                break
            if i % 4 != 1:
                continue
            s = line.strip()
            lens[len(s)] += 1
            for name, ad in KNOWN.items():
                if ad in s:
                    hits[name] += 1
            for j in range(a.from_pos, max(a.from_pos + 1, len(s) - 13)):
                kmers[s[j:j + 13]] += 1

    n = sum(lens.values())
    if not n:
        sys.exit(f"{a.fastq}: no reads")
    print(f"{a.fastq}: {n:,} reads sampled")
    print(f"read lengths (top 5): {dict(lens.most_common(5))}")
    print("\nknown adapters:")
    worst = 0.0
    for name, ad in KNOWN.items():
        pct = 100 * hits[name] / n
        worst = max(worst, pct if name != "polyA" else 0.0)
        print(f"  {hits[name]:>9,}  {pct:6.2f}%   {name:<28} {ad}")
    print(f"\nmost frequent 13-mers at position >= {a.from_pos} (unbiased):")
    for k, c in kmers.most_common(8):
        print(f"  {c:>10,}  {k}")
    print()
    if worst < 1.0:
        print("VERDICT: no adapter detected -- the library is already trimmed.")
        print("         Use `cutadapt --trim-n -m 20` with NO -a. Do not use --discard-untrimmed.")
    else:
        best = max((k for k in KNOWN if k != "polyA"), key=lambda k: hits[k])
        print(f"VERDICT: {best} present in {100*hits[best]/n:.1f}% of reads.")
        print(f"         Use `cutadapt -a {KNOWN[best]} ...`")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
