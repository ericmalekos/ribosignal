#!/usr/bin/env python3
# SUPERSEDED 2026-09-03: the real Hepatocytes RNA-seq BAMs exist in the warm archive
# (data/rnaseq_bam_mm1, "genomic + transcriptome", 165 files) and are used instead. This script is
# kept only for the case where a pack has no restorable RNA BAM. Do not use it when real reads are
# available: RiboTaper reads RNA for expression support, and a synthetic approximation there is an
# avoidable source of difference.
"""Synthesise a genome-coordinate RNA-seq BAM reproducing a pack's per-nt RNA coverage.

RiboTaper structurally requires an RNA bam (Ribotaper.sh argument 2; create_tracks.bash builds
RNA_tracks with `coverageBed -d -split -abam RNA_best.bam`). RiboCode and Ribo-TISH do not, so
this file exists only so RiboTaper can run at all. It is built ONCE from the pack's OBSERVED RNA
coverage and reused unchanged across every Ribo-seq variant, so it can never be the thing that
differs between the arms being compared.

SCHEME. Per transcript, N = round(scale * sum(c) / W) reads of length W are emitted with start
positions drawn from p(s) proportional to c[s]. Expected coverage at i is then
scale * c[i], so the shape is right in expectation at any depth and --scale sets the depth
without distorting it. --report_recovery measures the realised per-transcript error rather than
assuming it.

An earlier version emitted max(0, c[s]-c[s-1]) reads at each start, the positive first difference.
It was discarded on measurement, not on taste: at --scale 0.01 it produced only 23 of 200
transcripts and a mean absolute error of 0.90x mean coverage, because sub-1 coverage floors to
zero starts and RNA coverage fluctuates repeatedly inside one read length.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from synth_ribo_bam import cigar_from_blocks, load_exons, tx_to_genome_blocks  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", required=True, help="pack dir with coverage.npy/offsets/tx_order")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", default="")
    ap.add_argument("--read_len", type=int, default=50)
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply target coverage before differencing (0.01 keeps the BAM sane)")
    ap.add_argument("--report_recovery", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    pack = Path(a.pack)
    order = [l.strip() for l in (pack / "tx_order.txt").open()]
    offs = np.load(pack / "offsets.npy")
    cov = np.load(pack / "coverage.npy", mmap_mode="r")
    exons = load_exons(Path(a.gtf), set(order))
    chroms = set(x for x in a.chroms.split(",") if x)

    refs, lens = [], []
    for line in Path(a.fai).open():
        f = line.split("\t"); refs.append(f[0]); lens.append(int(f[1]))
    tid = {r: i for i, r in enumerate(refs)}
    hdr = {"HD": {"VN": "1.6", "SO": "unsorted"},
           "SQ": [{"SN": r, "LN": ln} for r, ln in zip(refs, lens)],
           "PG": [{"ID": "synth_rna_bam", "PN": "synth_rna_bam.py",
                   "CL": " ".join(sys.argv), "VN": "1"}]}

    rng = np.random.default_rng(a.seed)
    W = a.read_len
    n_reads = n_tx = 0
    errs = []
    with pysam.AlignmentFile(a.out, "wb", header=hdr) as bam:
        for j, t in enumerate(order):
            e = exons.get(t)
            if e is None:
                continue
            chrom, strand, blocks = e
            if (chroms and chrom not in chroms) or chrom not in tid:
                continue
            c = np.asarray(cov[offs[j]:offs[j + 1]], dtype=np.float64)
            L = c.size
            if sum(b - s for s, b in blocks) != L or L < W:
                continue
            tot = float(c.sum())
            if tot <= 0:
                continue
            n = int(round(a.scale * tot / W))
            if n <= 0:
                continue
            pstart = c[:L - W + 1]
            ps = pstart.sum()
            if ps <= 0:
                continue
            draws = rng.choice(L - W + 1, size=n, p=pstart / ps)
            starts = np.bincount(draws, minlength=L).astype(np.int64)
            nz = np.nonzero(starts)[0]
            if nz.size == 0:
                continue
            n_tx += 1
            if len(errs) < a.report_recovery:
                k = np.convolve(starts, np.ones(W, dtype=np.int64), mode="full")[:L]
                target = a.scale * c
                denom = max(1e-9, target.mean())
                c = target
                errs.append(float(np.abs(k - c).mean() / denom))
            for s in nz:
                s = int(s)
                if s + W > L:
                    continue
                bl = tx_to_genome_blocks(blocks, strand, s, W)
                if bl is None:
                    continue
                r = pysam.AlignedSegment()
                r.query_name = f"{t}_r{s}"
                r.flag = 0 if strand == "+" else 16
                r.reference_id = tid[chrom]
                r.reference_start = bl[0][0]
                r.mapping_quality = 255
                r.cigartuples = cigar_from_blocks(bl)
                r.query_sequence = "N" * W
                r.query_qualities = pysam.qualitystring_to_array("I" * W)
                r.set_tag("NH", 1)
                for _ in range(int(starts[s])):
                    bam.write(r)
                n_reads += int(starts[s])
    print(f"wrote {a.out}  reads={n_reads:,}  tx={n_tx:,}", file=sys.stderr)
    if errs:
        e = np.array(errs)
        print(f"coverage recovery over {len(e)} tx: mean abs error / mean coverage = "
              f"{e.mean():.4f}  median {np.median(e):.4f}  p95 {np.percentile(e,95):.4f}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
