#!/usr/bin/env python3
"""Regression tests for the ORF-caller tooling. Each one checks a bug that has already happened.

Every case here produced a plausible-looking result at the time. The minus-strand offset bug wrote
a complete, valid, sorted bam; the read-loss bug produced 105 million perfectly good genome records
while silently discarding 89 million more. Nothing looked wrong until it was measured.

Synthetic fixtures only. No bam, no genome, no SLURM, no real data. Runs in well under a second.

    python tests/test_orfcallers.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts" / "orfcallers"))

# imported without pysam/numpy: these are the pure-python coordinate routines
from synth_ribo_bam import (  # noqa: E402
    cigar_from_blocks,
    psite_window_start,
    tx_to_genome_blocks,
)

FAILURES = []


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


# A two-exon transcript. Genome blocks are always stored genome-ascending; TRANSCRIPT order is
# reversed for a minus-strand gene, which is where the sign errors live.
PLUS = [(1000, 1030), (2000, 2030)]          # 30 + 30 nt, transcript 5' at 1000
MINUS = [(2000, 2030), (1000, 1030)]         # transcript order: 5' at 2029, running down


def test_tx_to_genome_plus():
    check("plus strand: first base maps to the first exon start",
          tx_to_genome_blocks(PLUS, "+", 0, 1) == [(1000, 1001)])
    check("plus strand: last base of exon 1",
          tx_to_genome_blocks(PLUS, "+", 29, 1) == [(1029, 1030)])
    check("plus strand: first base of exon 2 is transcript position 30",
          tx_to_genome_blocks(PLUS, "+", 30, 1) == [(2000, 2001)])
    check("plus strand: a window spanning the junction returns TWO blocks",
          tx_to_genome_blocks(PLUS, "+", 28, 4) == [(1028, 1030), (2000, 2002)])


def test_tx_to_genome_minus():
    check("minus strand: transcript base 0 is the HIGHEST genome coordinate",
          tx_to_genome_blocks(MINUS, "-", 0, 1) == [(2029, 2030)])
    check("minus strand: transcript advances DOWN the genome",
          tx_to_genome_blocks(MINUS, "-", 1, 1) == [(2028, 2029)])
    check("minus strand: crossing into exon 2 lands at the top of the lower exon",
          tx_to_genome_blocks(MINUS, "-", 30, 1) == [(1029, 1030)])
    check("minus strand: a junction-spanning window returns genome-ASCENDING blocks",
          tx_to_genome_blocks(MINUS, "-", 28, 4) == [(1028, 1030), (2000, 2002)])


def test_runs_off_the_end_returns_none():
    """A read window past the transcript end must be REJECTED, not silently truncated.
    Truncating would place a short read at a real coordinate and it would look fine."""
    check("window past the 3' end returns None",
          tx_to_genome_blocks(PLUS, "+", 58, 4) is None)
    check("window exactly filling the transcript is accepted",
          tx_to_genome_blocks(PLUS, "+", 56, 4) == [(2026, 2030)])


def test_cigar_inserts_N_at_the_junction():
    check("contiguous blocks give a single M",
          cigar_from_blocks([(1000, 1029)]) == [(0, 29)])
    check("split blocks give M / N / M with the intron length",
          cigar_from_blocks([(1028, 1030), (2000, 2002)]) == [(0, 2), (3, 970), (0, 2)])


def test_psite_offset_is_measured_from_the_read_5_prime_end():
    """THE BUG: the minus-strand read window was computed as i-(R-1-O) instead of i-O, and 36,198
    of 77,682 emitted reads landed on the wrong genome base. The emitter and the P-site recovery
    used different conventions, so nothing failed; the reads were simply in the wrong place.

    Contract: the window starts at i-O in TRANSCRIPT coordinates on BOTH strands, and RiboTaper
    recovers the P-site at genome-leftmost + O on '+' and + (R-1-O) on '-'
    (libexec/P_sites_RNA_sites_calc.bash).
    """
    R, O = 29, 12
    for strand, blocks in (("+", PLUS), ("-", MINUS)):
        for i in (14, 20, 33, 41):                       # includes junction-spanning windows
            ts = psite_window_start(i, R, O, strand)      # the REAL emitter path
            bl = tx_to_genome_blocks(blocks, strand, ts, R)
            if bl is None:
                continue
            want = tx_to_genome_blocks(blocks, strand, i, 1)[0][0]
            walk = O if strand == "+" else (R - 1 - O)
            got = None
            for s, e in bl:
                if walk < e - s:
                    got = s + walk
                    break
                walk -= e - s
            check(f"P-site round-trips on strand {strand} at transcript position {i}",
                  got == want, f"got {got} want {want}")


def test_secondary_alignments_must_not_be_dropped_before_the_gtf_check():
    """THE BUG: tx_bam_to_genome dropped 0x100 records BEFORE testing whether the transcript was in
    the restricted GTF. STAR picks the primary arbitrarily among isoform copies, so a read whose
    primary happened to sit on a non-universe isoform was lost entirely even though a secondary sat
    on a universe transcript at the same locus. On real data this discarded 89,034,820 primaries,
    45.8% of all primaries, and a controlled test recovered 78.6% of them from their secondaries.

    This asserts the ORDER in the source, since the behaviour needs a bam to exercise.
    """
    src = (REPO / "scripts/orfcallers/tx_bam_to_genome.py").read_text()
    i_sec = src.index("if rd.is_secondary")
    i_gtf = src.index("e = ex.get(inb.get_reference_name")
    check("secondaries are only skipped under an explicit --primary-only opt-in",
          "a.primary_only" in src and "keep_secondary" not in src)
    check("the default path reaches the GTF membership test for secondaries too",
          i_sec < i_gtf and "not a.primary_only" not in src)


def test_class_map_is_resolved_against_the_caller_that_owns_it():
    """THE BUG: canon_of scanned every caller's map and returned the first hit, ignoring which
    caller it was asked about. `uORF` and `dORF` appear in both the RiboCode and RiboTaper maps.
    They agree today, so nothing was wrong yet, which is exactly why it needed a test."""
    from compare_three_callers import CANON, canon_of
    check("canon_of takes the caller as its first argument",
          canon_of("RiboTaper", "ORFs_ccds") == "canonical")
    check("a name absent from THAT caller's map resolves to None, not another caller's answer",
          canon_of("RiboCode", "ORFs_ccds") is None)
    shared = [k for k in CANON["RiboCode"] if k in CANON["RiboTaper"]]
    check("shared category names still agree across callers",
          all(CANON["RiboCode"][k] == CANON["RiboTaper"][k] for k in shared),
          f"shared: {shared}")


def main():
    for fn in (test_tx_to_genome_plus, test_tx_to_genome_minus,
               test_runs_off_the_end_returns_none, test_cigar_inserts_N_at_the_junction,
               test_psite_offset_is_measured_from_the_read_5_prime_end,
               test_secondary_alignments_must_not_be_dropped_before_the_gtf_check,
               test_class_map_is_resolved_against_the_caller_that_owns_it):
        print(f"\n{fn.__name__}")
        fn()
    print(f"\n{'FAILED: ' + ', '.join(FAILURES) if FAILURES else 'all passed'}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
