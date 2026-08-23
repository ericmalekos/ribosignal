#!/usr/bin/env python3
"""Shared BAM/alignment primitives. USE THESE -- do not hand-roll NH parsing or record counting.

This module exists because the same two mistakes were made repeatedly in this project, each time
producing a plausible-looking number that was wrong:

  1. Treating NH in a --quantMode TranscriptomeSAM BAM as GENOMIC multimapping. It is not. STAR
     expands ONE genomic alignment across every compatible isoform, so a genomically-UNIQUE read
     routinely carries NH:i:17. Splitting on NH measures ISOFORM multiplicity (~99% under posture A)
     and says nothing about what --outFilterMultimapNmax removed. Cost: a whole analysis reporting
     "99.8% of reads discarded" against STAR's own logs saying 30%.
  2. Comparing counts across BAMs produced with different multimap caps WITHOUT holding isoform
     expansion constant.

The sound comparison of two multimap postures is TOTAL RECORDS PER TRANSCRIPT, mm_a vs mm_b:
isoform expansion inflates both sides identically and cancels. `multimap_cost()` does that.
"""
from __future__ import annotations

import gzip
import subprocess
from collections import defaultdict
from pathlib import Path

SAMTOOLS = "/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools"


def is_transcriptome_bam(bam: str | Path) -> bool:
    """True if the @SQ names look like transcripts (ENST/ENSMUST), not chromosomes."""
    h = subprocess.run([SAMTOOLS, "view", "-H", str(bam)], capture_output=True, text=True).stdout
    for line in h.splitlines():
        if line.startswith("@SQ"):
            return "ENST" in line or "ENSMUST" in line
    return False


def quickcheck(bam: str | Path) -> bool:
    """EOF-block validity. A truncated BAM keeps a VALID HEADER, so a header test is NOT a
    validity test -- that is how job 36834132 died after 92 minutes (2026-08-21)."""
    return subprocess.run([SAMTOOLS, "quickcheck", str(bam)]).returncode == 0


def records_per_transcript(bam: str | Path, threads: int = 4) -> dict[str, int]:
    """tx_id -> TOTAL alignment records. The correct unit for comparing multimap postures.

    Deliberately does NOT split on NH: see the module docstring. If you think you want an
    NH split on a transcriptome BAM, you almost certainly want `multimap_cost` instead.
    """
    if not quickcheck(bam):
        raise RuntimeError(f"{bam} fails samtools quickcheck (truncated?)")
    p = subprocess.Popen([SAMTOOLS, "view", "-@", str(threads), str(bam)],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    out: dict[str, int] = defaultdict(int)
    for line in p.stdout:
        i = line.find("\t")
        j = line.find("\t", i + 1)
        k = line.find("\t", j + 1)
        out[line[j + 1:k]] += 1
    p.wait()
    if p.returncode:
        raise RuntimeError(f"samtools view failed on {bam}")
    return dict(out)


def load_counts_tsv(path: str | Path, value_cols: tuple[int, ...] = (1,)) -> dict[str, int]:
    """Read a gzipped `tx_id \\t v1 [\\t v2 ...]` count table, summing the requested columns."""
    out: dict[str, int] = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(value_cols):
                continue
            out[f[0]] = sum(int(f[c]) for c in value_cols)
    return out


def multimap_cost(counts_strict: dict[str, int], counts_loose: dict[str, int]) -> dict[str, float]:
    """Per transcript: fraction of its LOOSE-posture signal that the STRICT posture discarded.

        cost = 1 - strict/loose        in [0, 1]

    `counts_strict` = e.g. mm1 records per transcript; `counts_loose` = e.g. mm25. Both must be
    TOTAL record counts from the SAME transcriptome, so isoform expansion cancels. A transcript
    absent from the strict set has cost 1.0 (all of its signal was filtered away).
    """
    out = {}
    for t, loose in counts_loose.items():
        if loose <= 0:
            continue
        strict = counts_strict.get(t, 0)
        out[t] = max(0.0, 1.0 - strict / loose)
    return out
