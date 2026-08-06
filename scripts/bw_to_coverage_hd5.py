#!/usr/bin/env python3
"""Convert a deepTools bamCoverage bigWig (per-nt, per-transcript) to the coverage hd5.

bamCoverage writes a per-transcript per-nucleotide read-depth track (binSize 1,
--filterRNAstrand forward = ISR sense, secondary alignments counted = posture A) as a
bigWig whose "chromosomes" are the transcripts. This reads it back into the same vlen hd5
schema as rnaseq_coverage.py so the two counters are interchangeable and poolable.

Usage: bw_to_coverage_hd5.py <cov.bw> <out_coverage.hd5> <SRR>

Output hd5: transcript_ids (vlen-str, bigWig chrom order = BAM @SQ order) +
coverage (vlen-int32, per-nt read depth). Attrs: srr, libtype, n_transcripts,
total_coverage, source.
"""
import sys

import h5py
import numpy as np
import pyBigWig


def main():
    bw_path, out_hd5, srr = sys.argv[1], sys.argv[2], sys.argv[3]
    bw = pyBigWig.open(bw_path)
    chroms = bw.chroms()                       # dict tx_id -> length (bigWig order)
    refs = list(chroms.keys())
    n = len(refs)

    str_dt = h5py.string_dtype("utf-8")
    obj = np.empty(n, dtype=object)
    total = 0
    covered = 0
    for i, tx in enumerate(refs):
        L = int(chroms[tx])
        vals = bw.values(tx, 0, L, numpy=True)         # float32, nan where no coverage
        a = np.rint(np.nan_to_num(vals, nan=0.0)).astype(np.int32)
        obj[i] = a
        s = int(a.sum(dtype=np.int64))
        total += s
        if s > 0:
            covered += 1
    bw.close()

    with h5py.File(out_hd5, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(refs, dtype=object), dtype=str_dt)
        o.create_dataset("coverage", data=obj, dtype=h5py.vlen_dtype(np.int32),
                         compression="gzip", compression_opts=4)
        o.attrs["srr"] = srr
        o.attrs["libtype"] = "ISR"
        o.attrs["n_transcripts"] = n
        o.attrs["total_coverage"] = int(total)
        o.attrs["source"] = "deeptools_bamCoverage_binSize1_filterRNAstrand_forward_postureA"

    print(f"{srr}: covered_tx={covered:,}/{n:,}  total_coverage={total:,}  wrote {out_hd5}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
