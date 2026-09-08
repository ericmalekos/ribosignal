#!/usr/bin/env python3
"""Synthesise a genome-coordinate Ribo-seq BAM whose P-sites reproduce a per-nt transcript
density array, so RiboTaper can be run UNMODIFIED on the model's predicted profiles.

WHY A BAM AND NOT A DROP-IN. Ribo-TISH exposes `--inprofile` and takes a transcript profile
directly. RiboTaper has no such hook: `Ribotaper.sh` derives P_sites_all and Centered_RNA from
BAMs, and `create_tracks.bash` additionally builds RIBO_tracks / RNA_tracks with
`coverageBed -d -split -abam`. Two of the four tracks the ORF finder reads therefore come from a
BAM. Hand-writing all four track files would mean reverse-engineering RiboTaper's on-disk format
and would fail silently if a field were wrong, so the reads are synthesised instead and every
downstream stage runs as it does on real data.

THE P-SITE CONTRACT. Reads are emitted at a single length R with a single offset O, and
Ribotaper.sh is then called with <read_lengths>=R <cutoffs>=O. RiboTaper recovers the P-site as
"the base O nt into the read along the transcript", which by construction is the position the
density array put it at. Emitting one length removes any dependence on RiboTaper's per-length
offset handling.

Reads are spliced when the R-nt window crosses an exon junction; the CIGAR carries the N blocks,
which is what `coverageBed -split` and RiboTaper's bed12 P-site awk both expect.

VERIFY, DO NOT ASSUME: --selfcheck re-derives the P-site from each emitted read exactly as
RiboTaper's awk does and asserts it lands on the intended genome base. A silent off-by-one here
would shift every ORF's frame.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# pysam is imported inside main(), not at module scope, so the pure coordinate routines
# (tx_to_genome_blocks, cigar_from_blocks, project, build_cigar) stay importable in an
# interpreter without pysam. tests/test_orfcallers.py depends on that.

ATTR = re.compile(r'(\S+) "([^"]*)"')


def load_exons(gtf: Path, want: set[str]):
    """tx -> (chrom, strand, [(g_start0, g_end) ...] in TRANSCRIPT order)."""
    ex = defaultdict(list)
    info = {}
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            t = a.get("transcript_id")
            if t is None or (want and t not in want):
                continue
            ex[t].append((int(f[3]) - 1, int(f[4])))
            info[t] = (f[0], f[6])
    out = {}
    for t, blocks in ex.items():
        chrom, strand = info[t]
        blocks.sort()
        if strand == "-":
            blocks = blocks[::-1]
        out[t] = (chrom, strand, blocks)
    return out


def tx_to_genome_blocks(blocks, strand, t_start, t_len):
    """Map transcript interval [t_start, t_start+t_len) to sorted genome blocks."""
    got, need, skip = [], t_len, t_start
    for gs, ge in blocks:
        blen = ge - gs
        if skip >= blen:
            skip -= blen
            continue
        take = min(blen - skip, need)
        if strand == "+":
            s = gs + skip
            got.append((s, s + take))
        else:
            e = ge - skip
            got.append((e - take, e))
        need -= take
        skip = 0
        if need == 0:
            break
    if need:
        return None
    got.sort()
    merged = [list(got[0])]
    for s, e in got[1:]:
        if s == merged[-1][1]:
            merged[-1][1] = e
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def psite_window_start(i, read_len, offset, strand):
    """Transcript start of the read whose P-site lands on transcript position `i`.

    The P-site offset is measured from the read's 5' END on BOTH strands, so this is i - offset
    regardless of strand. RiboTaper agrees: for '-' it sets the genome start to
    leftmost + (read_len - offset - 1), and the genome-right end IS the 5' end of a minus-strand
    read (libexec/P_sites_RNA_sites_calc.bash).

    `strand` is accepted and unused on purpose. An earlier version branched on it and used
    i - (read_len - 1 - offset) for '-', which put 36,198 of 77,682 emitted reads on the wrong
    genome base while producing a completely valid bam.
    """
    del strand
    return i - offset


def cigar_from_blocks(bl):
    c = []
    for i, (s, e) in enumerate(bl):
        if i:
            c.append((3, s - bl[i - 1][1]))   # N
        c.append((0, e - s))                  # M
    return c


def main() -> int:
    import pysam
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True, help="genome .fai, for the BAM header")
    ap.add_argument("--variant", required=True,
                    choices=["real", "pred_obsdepth", "pred_preddepth"])
    ap.add_argument("--out", required=True, help="output BAM (unsorted; sort/index after)")
    ap.add_argument("--chroms", default="", help="comma-separated chromosome restriction")
    ap.add_argument("--read_len", type=int, default=29)
    ap.add_argument("--offset", type=int, default=12)
    ap.add_argument("--pred_scale", type=float, default=1.0)
    ap.add_argument("--pred_poisson", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ref-seq", default=None,
                    help="genome fasta. With it, each emitted read carries the REFERENCE sequence "
                         "at its projected position instead of a run of N. PRICE needs this: it "
                         "estimates its cleavage model from read mismatches, and an all-N read "
                         "makes every base a mismatch, which aborts the estimator.")
    ap.add_argument("--selfcheck", type=int, default=20000,
                    help="re-derive the P-site for this many emitted reads and assert it matches")
    a = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from ribocode_dropin import build_density

    fa = pysam.FastaFile(a.ref_seq) if a.ref_seq else None
    z = np.load(a.profiles, allow_pickle=True)
    ids = list(map(str, z["tx_ids"])); lens = z["lengths"]
    off = np.concatenate([[0], np.cumsum(lens)])
    pf, ob, pt = z["pred_flat"], z["obs_flat"], z["pred_total"]
    rng = np.random.default_rng(a.seed) if a.pred_poisson else None

    keep = set(ids)
    exons = load_exons(Path(a.gtf), keep)
    chroms = set(x for x in a.chroms.split(",") if x)
    refs, lengths = [], []
    for line in Path(a.fai).open():
        f = line.split("\t")
        refs.append(f[0]); lengths.append(int(f[1]))
    tid = {r: i for i, r in enumerate(refs)}
    hdr = {"HD": {"VN": "1.6", "SO": "unsorted"},
           "SQ": [{"SN": r, "LN": ln} for r, ln in zip(refs, lengths)],
           "PG": [{"ID": "synth_ribo_bam", "PN": "synth_ribo_bam.py",
                   "CL": " ".join(sys.argv), "VN": "1"}]}

    R, O = a.read_len, a.offset
    n_reads = n_tx = n_skip_edge = n_skip_map = n_skip_tx = 0
    checked = bad = 0
    with pysam.AlignmentFile(a.out, "wb", header=hdr) as bam:
        for j, t in enumerate(ids):
            e = exons.get(t)
            if e is None:
                n_skip_tx += 1
                continue
            chrom, strand, blocks = e
            if chroms and chrom not in chroms:
                continue
            if chrom not in tid:
                n_skip_tx += 1
                continue
            L = int(lens[j])
            if sum(b - s for s, b in blocks) != L:
                n_skip_tx += 1
                continue
            aa, bb = off[j], off[j + 1]
            dens = build_density(a.variant, pf[aa:bb], ob[aa:bb], pt[j],
                                 pred_scale=a.pred_scale, rng=rng, poisson=a.pred_poisson)
            nz = np.nonzero(dens)[0]
            if nz.size == 0:
                continue
            n_tx += 1
            for i in nz:
                i = int(i)
                # The P-site offset is measured from the read's 5' END on both strands, so in
                # TRANSCRIPT coordinates the window always starts at i - O. RiboTaper's awk
                # agrees: for '-' it sets start = leftmost + (readlen - offset - 1), which is the
                # base O in from the genome-RIGHT end, and the genome-right end is the 5' end of
                # a minus-strand read. An earlier version used i - (R-1-O) here and the
                # --selfcheck caught it: 36,198 of 77,682 reads landed on the wrong base.
                ts = psite_window_start(i, R, O, strand)
                if ts < 0 or ts + R > L:
                    n_skip_edge += 1
                    continue
                bl = tx_to_genome_blocks(blocks, strand, ts, R)
                if bl is None:
                    n_skip_map += 1
                    continue
                cig = cigar_from_blocks(bl)
                if checked < a.selfcheck:
                    # RiboTaper walks O bases into the read ALONG THE TRANSCRIPT
                    want = tx_to_genome_blocks(blocks, strand, i, 1)[0][0]
                    # genome-leftmost + O on '+', + (R-1-O) on '-', matching
                    # libexec/P_sites_RNA_sites_calc.bash
                    walk, got = (O if strand == "+" else (R - 1 - O)), None
                    for s, en in bl:
                        if walk < en - s:
                            got = s + walk; break
                        walk -= en - s
                    checked += 1
                    if got != want:
                        bad += 1
                r = pysam.AlignedSegment()
                r.query_name = f"{t}_{i}"
                r.flag = 0 if strand == "+" else 16
                r.reference_id = tid[chrom]
                r.reference_start = bl[0][0]
                r.mapping_quality = 255
                r.cigartuples = cig
                if fa is None:
                    r.query_sequence = "N" * R
                else:
                    # SAM stores SEQ in REFERENCE-FORWARD orientation whatever the strand flag,
                    # so this must NOT be reverse-complemented: doing so would make every base a
                    # mismatch and reproduce the exact PRICE failure this option exists to fix.
                    sq = "".join(fa.fetch(chrom, s0, e0) for s0, e0 in bl).upper()
                    r.query_sequence = sq if len(sq) == R else "N" * R
                r.query_qualities = pysam.qualitystring_to_array("I" * R)
                r.set_tag("NH", 1)
                c = int(dens[i])
                for _ in range(c):
                    bam.write(r)
                n_reads += c
    print(f"wrote {a.out}  reads={n_reads:,}  tx={n_tx:,}  "
          f"skipped: tx={n_skip_tx:,} edge={n_skip_edge:,} unmappable={n_skip_map:,}",
          file=sys.stderr)
    print(f"selfcheck: {checked:,} reads re-derived, {bad} P-site mismatches", file=sys.stderr)
    if bad:
        print("FATAL: P-site round-trip failed; the profile would be frame-shifted.",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
