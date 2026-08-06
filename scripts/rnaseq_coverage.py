#!/usr/bin/env python3
"""Per-nucleotide RNAseq read coverage from a STAR toTranscriptome BAM (one sample).

Coverage is the model's second input feature (RNA abundance/context), built on the SAME
transcriptome-nucleotide axis as the Ribo-seq P-site target and the FM per-token
embeddings, so all three align 1:1 position-for-position.

Counting rules:
  - Strand: AUTO-DETECTED per library (a first pass samples mate orientations). On a
    transcriptome BAM the sense reads for a transcript are one consistent orientation; the
    antisense ~1% are contamination and must be dropped. "ISR-sense" = `is_read2 != is_reverse`
    (read2-forward / read1-reverse). If the sampled ISR-sense fraction is >= 0.7 the library
    is reverse-stranded (keep ISR-sense, = the original Chothani behavior); <= 0.3 it is
    forward-stranded (keep the complement, ISF); in between it is unstranded (keep all mapped).
    This generalizes the previously-hardcoded ISR rule so held-out libraries with a different
    prep (e.g. Ruiz-Orera forward-stranded, or a single-end mouse library) are counted on their
    OWN sense strand instead of silently discarding ~99% of reads. Reduces to ISR for Chothani,
    so the training coverage is unchanged. For single-end reads is_read2 is always False, so the
    same rule keys on is_reverse (reverse-stranded SE -> ISR branch, forward SE -> ISF branch).
  - Multimapper posture A (matches the P-site target): no NH filter; every alignment
    record (primary + secondary) contributes, so a fragment compatible with N isoforms
    covers all N. Skip unmapped + supplementary.
  - Depth semantics: each kept mate adds 1 to every transcript position it spans
    (reference_start .. reference_end, exclusive), like `samtools depth` per-mate.

Usage: rnaseq_coverage.py <toTranscriptome.bam> <out_coverage.hd5> <SRR>

Output hd5 (vlen schema, RiboCode-compatible so one reader serves inputs and target):
  transcript_ids  vlen-str    (STAR @SQ order; join to the P-site target by tx-id)
  coverage        vlen-int32  (per-nt read depth, 0-based, length == tx length)
  attrs: srr, libtype (ISR|ISF|U), strand_isr_frac, strand_sample_n,
         n_transcripts, records_seen, mapped_sense_counted
"""
import os
import sys

import h5py
import numpy as np
import pysam

SAMPLE_N = 2_000_000   # mapped records to sample for strand auto-detection


def isr_sense(r):
    """True if r is the ISR-sense orientation (read2-forward / read1-reverse). For single-end
    reads is_read2 is always False, so this reduces to `is_reverse` (reverse-strand read)."""
    return r.is_read2 != r.is_reverse


def detect_strand(bam_path):
    """Sample mapped records, return (libtype, isr_fraction, n_sampled)."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    n_map = n_isr = 0
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_supplementary:
            continue
        n_map += 1
        if isr_sense(r):
            n_isr += 1
        if n_map >= SAMPLE_N:
            break
    bam.close()
    frac = (n_isr / n_map) if n_map else 0.0
    libtype = "ISR" if frac >= 0.7 else "ISF" if frac <= 0.3 else "U"
    return libtype, frac, n_map


def main():
    bam_path, out_hd5, srr = sys.argv[1], sys.argv[2], sys.argv[3]
    libtype, frac, nsamp = detect_strand(bam_path)
    print(f"{srr}: strand auto-detect ISR-sense frac={frac:.3f} over {nsamp:,} sampled "
          f"-> libtype={libtype}", file=sys.stderr)

    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = list(bam.references)
    lengths = list(bam.lengths)
    n = len(refs)

    cov = {}                       # ref_id -> int32 per-nt array (lazy, covered tx only)
    seen = counted = 0
    for r in bam.fetch(until_eof=True):
        seen += 1
        if r.is_unmapped or r.is_supplementary:
            continue
        # drop the antisense strand for a stranded library; keep everything if unstranded
        if libtype == "ISR" and not isr_sense(r):
            continue
        if libtype == "ISF" and isr_sense(r):
            continue
        e = r.reference_end
        if e is None:
            continue
        rid = r.reference_id
        a = cov.get(rid)
        if a is None:
            a = np.zeros(lengths[rid], dtype=np.int32)
            cov[rid] = a
        a[r.reference_start:e] += 1
        counted += 1
    bam.close()

    # GUARD: fail loudly on near-empty coverage instead of silently writing a useless hd5 that
    # poisons the pack. covered_tx = # transcripts with any read; a healthy RNA-seq sample covers
    # tens of thousands (the corrupt+adapter Janich RNA gave ~1,530). Floor: $RNASEQ_MIN_COVERED_TX.
    covered = len(cov)
    floor = int(os.environ.get("RNASEQ_MIN_COVERED_TX", "5000"))
    if covered < floor:
        print(f"ERROR {srr}: only {covered:,}/{n:,} tx covered (records_seen={seen:,}, "
              f"sense_counted={counted:,}); below floor {floor:,}. Likely corrupt/adapter "
              f"fastq, wrong STAR index, or over-strict filtering. NOT writing {out_hd5}.",
              file=sys.stderr)
        sys.exit(1)

    str_dt = h5py.string_dtype("utf-8")
    obj = np.empty(n, dtype=object)
    for i in range(n):
        a = cov.get(i)
        obj[i] = a if a is not None else np.zeros(lengths[i], dtype=np.int32)
    with h5py.File(out_hd5, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(refs, dtype=object), dtype=str_dt)
        o.create_dataset("coverage", data=obj, dtype=h5py.vlen_dtype(np.int32),
                         compression="gzip", compression_opts=4)
        o.attrs["srr"] = srr
        o.attrs["libtype"] = libtype
        o.attrs["strand_isr_frac"] = float(frac)
        o.attrs["strand_sample_n"] = int(nsamp)
        o.attrs["n_transcripts"] = n
        o.attrs["records_seen"] = int(seen)
        o.attrs["mapped_sense_counted"] = int(counted)

    print(f"{srr}: libtype={libtype} records_seen={seen:,}  sense_counted={counted:,}  "
          f"covered_tx={len(cov):,}/{n:,}  wrote {out_hd5}", file=sys.stderr)


if __name__ == "__main__":
    main()
