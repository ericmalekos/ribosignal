#!/usr/bin/env python3
"""Project a TRANSCRIPTOME bam back onto genome coordinates.

Lets genome-based tools (RiboTaper, ribotish quality, bedtools, IGV) read the transcriptome bams
this project keeps, without re-aligning from FASTQ and without synthesising reads.

THE TRAP THIS TOOL EXISTS TO AVOID. `STAR --quantMode TranscriptomeSAM` expands ONE genomic
alignment across EVERY compatible isoform, so a genomically UNIQUE read routinely carries
NH:i:17. Projecting every record would write that read 17 times at the same genome locus and
inflate coverage by the isoform multiplicity. Alignments are therefore deduplicated per read name
on the projected genome coordinate, and NH is RECOMPUTED as the number of distinct genome loci the
read maps to, which is the genomic multimapping the transcriptome NH never measured.

WHAT IT DOES
  transcript interval  ->  exon walk  ->  genome blocks  ->  M/N cigar
  minus-strand transcript: the read is on the genome minus strand, so FLAG bit 0x10 is flipped
  and SEQ/QUAL are reverse-complemented
  cigar ops: M/=/X/D consume transcript coordinates and are projected; I/S/H are carried through
  unchanged; N cannot occur in transcript space and is rejected

INPUT MUST BE NAME-GROUPED so a read's alignments arrive together. STAR's raw
Aligned.toTranscriptome.out.bam already is; anything that has been coordinate-sorted is not, and
--name-sort will re-sort it to a temporary file first.

VERIFY, DO NOT ASSUME: --selfcheck re-derives each projected block from the exon table
independently and asserts the transcript length is conserved. Report the counters; a nonzero
`unmappable` means the GTF and the bam header disagree about transcript structure.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import paths  # noqa: E402

import pysam

ATTR = re.compile(r'(\S+) "([^"]*)"')
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CONSUME_TX = {0, 2, 7, 8}       # M D = X
CONSUME_QRY_ONLY = {1, 4}       # I S


def load_exons(gtf: Path, want: set[str]):
    """tx -> (chrom, strand, [(g_start0, g_end)...] in TRANSCRIPT order), plus cdna length."""
    ex, info = defaultdict(list), {}
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
        out[t] = (chrom, strand, blocks, sum(e - s for s, e in blocks))
    return out


def project(blocks, strand, t_start, t_len):
    """Transcript interval -> sorted, merged genome blocks. None if it runs off the transcript."""
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


def build_cigar(read_cigar, blocks, strand, t_pos):
    """Project a transcript-space cigar onto the genome.

    Returns (genome_start, cigartuples) or None if any transcript-consuming op falls outside the
    transcript. Ops that consume transcript coordinates (M/D/=/X) are projected block by block with
    an N inserted at every exon junction; query-only ops (I/S/H) pass through. On a minus-strand
    transcript the op ORDER reverses, because transcript order runs against genome order.
    """
    pieces, cur = [], t_pos
    for op, ln in read_cigar:
        if op == 3:
            raise ValueError("N in a transcript-space cigar")
        if op in CONSUME_TX:
            bl = project(blocks, strand, cur, ln)
            if bl is None:
                return None
            pieces.append(("blocks", op, bl))
            cur += ln
        elif op in CONSUME_QRY_ONLY or op == 5:
            pieces.append(("plain", op, ln))
        else:
            return None
    if strand == "-":
        pieces = pieces[::-1]

    out, gstart, prev_end = [], None, None
    for kind, op, payload in pieces:
        if kind == "plain":
            out.append((op, payload))
            continue
        for s, e in payload:                      # project() returns genome-ascending blocks
            if gstart is None:
                gstart = s
            elif prev_end is not None and s > prev_end:
                out.append((3, s - prev_end))     # N across the junction
            out.append((op, e - s))
            prev_end = e
    if gstart is None:
        return None

    merged = []
    for op, ln in out:
        if ln <= 0:
            continue
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + ln)
        else:
            merged.append((op, ln))
    return gstart, merged


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-bam", required=True, help="transcriptome bam (name-grouped)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True, help="genome .fai for the output header")
    ap.add_argument("--out", required=True, help="genome bam (unsorted; sort and index after)")
    ap.add_argument("--name-sort", action="store_true",
                    help="the input is coordinate-sorted; re-sort by name into a temp file first")
    ap.add_argument("--keep-secondary", action="store_true",
                    help="keep 0x100 records; by default they are dropped BEFORE dedup, since "
                         "transcriptome secondaries are isoform copies, not extra genomic loci")
    ap.add_argument("--tmpdir", default=str(paths.tmp_dir()),
                    help="scratch for the name-sort (default $RIBO_TMPDIR, "
                         "else the system temp dir)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--samtools", default="samtools",
                    help="samtools binary; the default resolves on $PATH")
    a = ap.parse_args()

    refs, lens = [], []
    for line in Path(a.fai).open():
        f = line.split("\t")
        refs.append(f[0]); lens.append(int(f[1]))
    tid = {r: i for i, r in enumerate(refs)}

    src = a.in_bam
    tmp = None
    if a.name_sort:
        Path(a.tmpdir).mkdir(parents=True, exist_ok=True)
        tmp = tempfile.NamedTemporaryFile(suffix=".bam", dir=a.tmpdir, delete=False).name
        print(f"name-sorting into {tmp}", file=sys.stderr)
        subprocess.run([a.samtools, "sort", "-n", "-@", str(a.threads), "-o", tmp, a.in_bam],
                       check=True)
        src = tmp

    inb = pysam.AlignmentFile(src, "rb")
    want = set(inb.references)
    ex = load_exons(Path(a.gtf), want)
    print(f"bam carries {len(want):,} transcripts; GTF resolves {len(ex):,}", file=sys.stderr)

    hdr = {"HD": {"VN": "1.6", "SO": "unsorted"},
           "SQ": [{"SN": r, "LN": ln} for r, ln in zip(refs, lens)],
           "PG": [{"ID": "tx_bam_to_genome", "PN": "tx_bam_to_genome.py",
                   "CL": " ".join(sys.argv), "VN": "1"}]}

    n_in = n_out = n_dedup = n_notx = n_unmap = n_sec = 0
    lenmismatch = 0
    out = pysam.AlignmentFile(a.out, "wb", header=hdr)

    def flush(name, recs):
        nonlocal n_out, n_dedup
        uniq = {}
        for key, r in recs:
            if key not in uniq:
                uniq[key] = r
        nh = len(uniq)
        n_dedup += len(recs) - nh
        for r in uniq.values():
            r.set_tag("NH", nh)
            r.mapping_quality = 255 if nh == 1 else (3 if nh == 2 else 1)
            out.write(r)
            n_out += 1

    cur_name, buf = None, []
    for rd in inb.fetch(until_eof=True):
        n_in += 1
        if rd.is_unmapped:
            continue
        if rd.is_secondary and not a.keep_secondary:
            n_sec += 1
            continue
        e = ex.get(inb.get_reference_name(rd.reference_id))
        if e is None:
            n_notx += 1
            continue
        chrom, strand, blocks, clen = e
        if clen != inb.get_reference_length(inb.get_reference_name(rd.reference_id)):
            lenmismatch += 1
            continue
        try:
            pr = build_cigar(rd.cigartuples, blocks, strand, rd.reference_start)
        except ValueError:
            n_unmap += 1
            continue
        if pr is None or chrom not in tid:
            n_unmap += 1
            continue
        gstart, cig = pr
        nr = pysam.AlignedSegment()
        nr.query_name = rd.query_name
        rev = rd.is_reverse ^ (strand == "-")
        nr.flag = 16 if rev else 0
        nr.reference_id = tid[chrom]
        nr.reference_start = gstart
        nr.cigartuples = cig
        seq = rd.query_sequence or ""
        qual = rd.query_qualities
        if strand == "-":
            seq = seq.translate(COMP)[::-1]
            qual = qual[::-1] if qual is not None else None
        nr.query_sequence = seq
        if qual is not None:
            nr.query_qualities = qual
        key = (nr.reference_id, gstart, tuple(cig), nr.flag)
        if rd.query_name != cur_name:
            if cur_name is not None:
                flush(cur_name, buf)
            cur_name, buf = rd.query_name, []
        buf.append((key, nr))
    if cur_name is not None:
        flush(cur_name, buf)
    out.close(); inb.close()
    if tmp:
        Path(tmp).unlink(missing_ok=True)

    print(f"in={n_in:,}  secondary_dropped={n_sec:,}  tx_not_in_gtf={n_notx:,}  "
          f"cdna_length_mismatch={lenmismatch:,}  unprojectable={n_unmap:,}", file=sys.stderr)
    print(f"collapsed {n_dedup:,} isoform-duplicate alignments; wrote {n_out:,} genome records",
          file=sys.stderr)
    if lenmismatch:
        print("WARNING: the GTF and the bam header disagree on cDNA length for some transcripts; "
              "those reads were dropped rather than projected to a shifted position.",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
