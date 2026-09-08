#!/usr/bin/env python3
"""Per-nucleotide ribosome P-site counts from a Ribo-seq transcriptome BAM.

This is the missing half of the input pipeline. `rnaseq_coverage.py` turns an RNA-seq
transcriptome BAM into the model's coverage channel, and `build_pack.py --psites` takes a
P-site hd5 -- but nothing published here produced one, so the observed arm (`--variant
real`), the target for training, and any predicted-vs-observed comparison were all out of
reach from a public BAM.

    python scripts/ribo_psites.py --bam ribo.Aligned.toTranscriptome.out.bam \\
           --gtf annotation.gtf --out psites.hd5 --sample SRR15513208

The P-site is not the read. A ribosome footprint is ~28-30 nt and the P-site sits roughly
12 nt inside its 5' end, but the exact offset varies with read length and with the
nuclease digestion of the particular library, so it is CALIBRATED here rather than
assumed: for each read length, the offset is the one that piles the most 5' ends 12 nt
upstream of an annotated start codon, and the frame-0 fraction it achieves is reported as
the check. A library whose frame-0 fraction is near 1/3 has no periodicity and something
is wrong upstream -- wrong read-length window, RNA contamination, or not Ribo-seq at all.

Offsets are applied per read length, so `--lengths 25:35` keeps a window and each length
inside it gets its own offset. Lengths with too few reads to calibrate fall back to the
median of the calibrated ones and are named in the output.

Output matches rnaseq_coverage.py's schema exactly, so build_pack.py reads either:
  transcript_ids  vlen-str    (BAM @SQ order)
  psites          vlen-int32  (per-nt P-site counts, 0-based, length == tx length)
  attrs: sample, offsets_json, frame0_fraction, n_reads_counted, ...
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ATTR = re.compile(r'(\S+) "([^"]*)"')


def cds_starts_in_transcript_space(gtf_path):
    """(starts, annotated): {tx_id: 0-based transcript offset of the CDS first base},
    and the subset of those whose start is a real annotated `start_codon`.

    Needs the exon walk, not just a subtraction: the CDS start is a genome coordinate and
    the transcript is spliced, so the answer is the summed length of the exon pieces 5' of
    it. Strand-aware -- on the minus strand 'first base of the CDS' is the HIGHEST genome
    coordinate, and exons are walked from the top down.

    The `annotated` set matters for calibration and is not a formality. On GENCODE v49
    chr22, 272 of 5,754 CDS-bearing transcripts do not begin with ATG, and 261 of those are
    tagged mRNA_start_NF / cds_start_NF -- 5'-incomplete models whose CDS begins at
    transcript position 0 because the transcript itself is truncated, not because
    translation starts there. Calibrating a P-site offset against those would be
    calibrating against an annotation artifact. Restricted to `annotated`, the walk lands
    on ATG for 99.8% of transcripts (5,475 / 5,485), the remainder being GENCODE's
    genuinely non-AUG-initiated entries.
    """
    exons = defaultdict(list)      # tx -> [(start, end)] genome, 1-based inclusive
    cds = {}                       # tx -> (min_start, max_end)
    sc = {}                        # tx -> (min_start, max_end) of start_codon
    strand = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS", "start_codon"):
                continue
            a = dict(ATTR.findall(f[8]))
            tx = a.get("transcript_id")
            if not tx:
                continue
            s, e = int(f[3]), int(f[4])
            if f[2] == "exon":
                exons[tx].append((s, e))
                strand[tx] = f[6]
            else:
                d = cds if f[2] == "CDS" else sc
                lo, hi = d.get(tx, (s, e))
                d[tx] = (min(lo, s), max(hi, e))

    def walk(tx, g):
        ex = sorted(exons.get(tx, []))
        if not ex:
            return None
        st = strand[tx]
        pos = 0
        for a, b in (ex if st == "+" else [(y, x) for x, y in reversed(ex)]):
            lo, hi = (a, b) if st == "+" else (b, a)
            if lo <= g <= hi:
                return pos + (g - lo if st == "+" else hi - g)
            pos += hi - lo + 1
        return None

    out, annotated = {}, set()
    for tx, (clo, chi) in cds.items():
        if tx not in strand:
            continue
        # genome coordinate of the CDS's FIRST translated base
        found = walk(tx, clo if strand[tx] == "+" else chi)
        if found is not None:
            out[tx] = found
            if tx in sc:
                annotated.add(tx)
    return out, annotated


def calibrate(bam_path, cds0, lengths, max_reads, min_count):
    """Per-read-length 5'-end offset, and the frame-0 fraction each one achieves.

    Scores each candidate offset by how much of the CDS signal it puts in frame 0. The
    offset that maximises periodicity is the one that is actually placing the P-site: a
    wrong offset smears the reads across all three frames.
    """
    import pysam
    lo, hi = lengths
    # per length: histogram of (5' end - CDS start) for reads near an annotated start
    hist = {L: Counter() for L in range(lo, hi + 1)}
    n_used = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_reverse or r.is_supplementary:
            continue                      # transcriptome reads are sense-strand by construction
        L = r.query_length
        if L is None or not (lo <= L <= hi):
            continue
        c = cds0.get(bam.get_reference_name(r.reference_id))
        if c is None:
            continue                      # cds0 here is the ANNOTATED-start subset (see main)
        d = r.reference_start - c
        if -40 <= d <= 40:
            hist[L][d] += 1
        n_used += 1
        if max_reads and n_used >= max_reads:
            break
    bam.close()

    offsets, frames, weak = {}, {}, []
    for L in range(lo, hi + 1):
        h = hist[L]
        tot = sum(h.values())
        if tot < min_count:
            weak.append(L)
            continue
        best, best_score = None, -1.0
        for off in range(8, 20):
            # a read of length L with offset `off` has its P-site at d + off relative to
            # the CDS start, so frame 0 means (d + off) % 3 == 0
            f0 = sum(n for d, n in h.items() if (d + off) % 3 == 0)
            score = f0 / tot
            if score > best_score:
                best, best_score = off, score
        offsets[L], frames[L] = best, best_score
    if offsets:
        med = int(np.median(list(offsets.values())))
        for L in weak:
            offsets[L] = med
    return offsets, frames, weak


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="Ribo-seq STAR toTranscriptome BAM")
    ap.add_argument("--gtf", required=True, help="the annotation the BAM was aligned against")
    ap.add_argument("--out", required=True, help="output .hd5")
    ap.add_argument("--sample", required=True, help="sample id, stored as an attribute")
    ap.add_argument("--lengths", default="25:35",
                    help="read-length window LO:HI to keep (default %(default)s)")
    ap.add_argument("--offset", type=int, default=None,
                    help="skip calibration and use this 5'-end offset for every length. "
                         "For reproducing a fixed-offset run; normally leave unset.")
    ap.add_argument("--calibrate-reads", type=int, default=5_000_000,
                    help="reads to sample for calibration (default %(default)s; 0 = all)")
    ap.add_argument("--min-count", type=int, default=1000,
                    help="a read length needs this many start-proximal reads to get its own "
                         "offset; below it, the median of the calibrated lengths is used")
    ap.add_argument("--min-frame0", type=float, default=0.4,
                    help="fail if the pooled frame-0 fraction is below this. 1/3 is no "
                         "periodicity at all; a real Ribo-seq library reaches 0.5-0.7.")
    args = ap.parse_args()

    import h5py
    import pysam

    lo, hi = (int(x) for x in args.lengths.split(":"))
    print(f"reading CDS starts from {args.gtf} ...", file=sys.stderr)
    cds0, annotated = cds_starts_in_transcript_space(args.gtf)
    cal0 = {t: c for t, c in cds0.items() if t in annotated}
    print(f"  {len(cds0):,} transcripts with a CDS; {len(cal0):,} with an annotated "
          f"start codon (calibration uses only these)", file=sys.stderr)
    if not cal0:
        sys.exit(f"{args.gtf} has no start_codon features, so no offset can be calibrated. "
                 f"Pass --offset to set one explicitly.")

    if args.offset is not None:
        offsets = {L: args.offset for L in range(lo, hi + 1)}
        frames, weak = {}, []
        print(f"using a fixed offset of {args.offset} for lengths {lo}-{hi}", file=sys.stderr)
    else:
        print("calibrating P-site offsets against annotated start codons ...", file=sys.stderr)
        offsets, frames, weak = calibrate(args.bam, cal0, (lo, hi),
                                          args.calibrate_reads, args.min_count)
        if not offsets:
            sys.exit(f"no read length in {lo}-{hi} had {args.min_count} start-proximal reads. "
                     f"Widen --lengths, lower --min-count, or check the BAM is Ribo-seq "
                     f"aligned to the same annotation as --gtf.")
        print(f"  {'len':>4} {'offset':>7} {'frame0':>8}", file=sys.stderr)
        for L in sorted(offsets):
            f = f"{frames[L]:.3f}" if L in frames else "(median)"
            print(f"  {L:>4} {offsets[L]:>7} {f:>8}", file=sys.stderr)
        if weak:
            print(f"  lengths given the median offset for lack of reads: {weak}", file=sys.stderr)

    bam = pysam.AlignmentFile(args.bam, "rb")
    refs, lens = list(bam.references), list(bam.lengths)
    counts = {}
    n_seen = n_counted = 0
    f0 = f_all = 0
    for r in bam.fetch(until_eof=True):
        n_seen += 1
        if r.is_unmapped or r.is_reverse or r.is_supplementary:
            continue
        L = r.query_length
        if L is None or L not in offsets:
            continue
        p = r.reference_start + offsets[L]
        rid = r.reference_id
        if not (0 <= p < lens[rid]):
            continue                       # P-site pushed off the transcript end
        a = counts.get(rid)
        if a is None:
            a = counts[rid] = np.zeros(lens[rid], dtype=np.int32)
        a[p] += 1
        n_counted += 1
        c = cds0.get(refs[rid])
        if c is not None:
            f_all += 1
            f0 += (p - c) % 3 == 0
    bam.close()

    frac0 = (f0 / f_all) if f_all else 0.0
    print(f"\n{args.sample}: {n_seen:,} records, {n_counted:,} P-sites over "
          f"{len(counts):,} transcripts", file=sys.stderr)
    print(f"  pooled frame-0 fraction {frac0:.3f} over {f_all:,} CDS-transcript P-sites "
          f"(1/3 = no periodicity)", file=sys.stderr)
    if frac0 < args.min_frame0:
        sys.exit(f"frame-0 fraction {frac0:.3f} is below --min-frame0 {args.min_frame0}. "
                 f"This library shows no 3-nt periodicity, so the P-site assignment is not "
                 f"meaningful and the counts would be a misleading target. Check the read "
                 f"lengths kept, the adapter trimming, and that this is Ribo-seq.")

    str_dt = h5py.string_dtype("utf-8")
    obj = np.empty(len(refs), dtype=object)
    for i in range(len(refs)):
        obj[i] = counts.get(i, np.zeros(lens[i], dtype=np.int32))
    with h5py.File(args.out, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(refs, dtype=object), dtype=str_dt)
        o.create_dataset("psites", data=obj, dtype=h5py.vlen_dtype(np.int32),
                         compression="gzip", compression_opts=4)
        o.attrs["sample"] = args.sample
        o.attrs["offsets_json"] = json.dumps({str(k): v for k, v in sorted(offsets.items())})
        o.attrs["frame0_fraction"] = float(frac0)
        o.attrs["n_records_seen"] = int(n_seen)
        o.attrs["n_psites_counted"] = int(n_counted)
        o.attrs["length_window"] = args.lengths
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
