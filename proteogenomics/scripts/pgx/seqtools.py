#!/usr/bin/env python3
"""Sequence helpers shared across `pgx` steps: FASTA iteration, translation, codon sets."""
from __future__ import annotations

import gzip

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
for _i, _c1 in enumerate(_BASES):
    for _j, _c2 in enumerate(_BASES):
        for _k, _c3 in enumerate(_BASES):
            CODON[_c1 + _c2 + _c3] = _AAS[_i * 16 + _j * 4 + _k]

STOPS = ("TAA", "TAG", "TGA")

# The nine single-nucleotide variants of ATG, i.e. the set with documented mammalian initiation
# activity. The model's ORF track already scores these (CUG 0.5, GUG/ACG 0.35, UUG 0.3, AUA 0.25,
# AUU/AUC 0.2, AAG/AGG 0.15).
NEAR_COGNATE = ("CTG", "GTG", "TTG", "ACG", "AGG", "AAG", "ATA", "ATT", "ATC")
ALL_STARTS = ("ATG",) + NEAR_COGNATE


def parse_starts(spec):
    """'ATG' | 'near_cognate' | 'ATG,CTG,GTG' -> tuple of validated DNA codons."""
    s = (spec or "").strip().lower()
    if s in ("near_cognate", "all"):
        return ALL_STARTS
    out = tuple(c.strip().upper().replace("U", "T") for c in str(spec).split(",") if c.strip())
    bad = [c for c in out if len(c) != 3 or set(c) - set("ACGT")]
    if bad:
        raise SystemExit(f"bad start codon(s): {bad}")
    return out


def norm(seq):
    return seq.upper().replace("U", "T")


def translate(s):
    """DNA -> protein. Unknown codons become X; the caller decides how to treat '*'."""
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def iter_fasta(path):
    """Yield (name, sequence); name is the token before the first '|' or whitespace."""
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = ln[1:].split("|")[0].split()[0].strip(); buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            yield name, "".join(buf)


def read_fasta(path, keep=None):
    """{name: sequence}, optionally restricted to the `keep` set."""
    return {n: s for n, s in iter_fasta(path) if keep is None or n in keep}


def candidate_orfs(seq, min_nt, starts=("ATG",)):
    """Every <start>..in-frame-stop ORF -> [(start0, end0_excl, start_codon)].

    One ORF per (start position, next in-frame stop), so nested ORFs sharing a stop but opening at
    different starts are ALL emitted. That is required for the null databases: a naive enumeration
    is exactly this set, and the model's value is the selection made from it.
    """
    import numpy as np
    s = norm(seq)
    L = len(s)
    if L < 6:
        return []
    b = np.frombuffer(s.encode("ascii", "replace"), dtype=np.uint8)
    T_, A_, G_ = 84, 65, 71
    out = []
    for f in range(3):
        n = (L - f) // 3
        if n < 2:
            continue
        idx = f + 3 * np.arange(n)
        c0, c1, c2 = b[idx], b[idx + 1], b[idx + 2]
        is_stop = (((c0 == T_) & (c1 == A_) & (c2 == A_)) | ((c0 == T_) & (c1 == A_) & (c2 == G_))
                   | ((c0 == T_) & (c1 == G_) & (c2 == A_)))
        stop_pos = idx[is_stop]
        if stop_pos.size == 0:
            continue
        for sc in starts:
            k0, k1, k2 = ord(sc[0]), ord(sc[1]), ord(sc[2])
            m = (c0 == k0) & (c1 == k1) & (c2 == k2)
            for a in idx[m]:
                j = int(np.searchsorted(stop_pos, a, side="right"))
                if j >= stop_pos.size:
                    continue
                e = int(stop_pos[j]) + 3
                if e - a >= min_nt:
                    out.append((int(a), e, sc))
    return out


def cds_relationship(a, e, cds):
    """Position of ORF [a, e) relative to the annotated CDS (cds_start0, cds_end0_excl).

    Out-of-frame overlaps are split three ways to MATCH RiboCode's `classfy_orf`, which treats
    Overlap_uORF and Overlap_dORF as novel and only a fully-contained out-of-frame ORF as
    `internal`. Collapsing all out-of-frame overlaps to `internal` (as an earlier version did)
    left the null databases missing an entire class that the RiboCode-derived model arm contains,
    which makes the model-vs-null comparison unmatched: measured on BMDM, 117 of 1,533 model
    sequences were absent from the near-cognate null for that reason alone, all of them ATG.
    """
    if cds is None:
        return "no_cds"
    cs, ce = cds
    inframe = (a - cs) % 3 == 0
    if inframe and a >= ce:
        return "dorf_inframe"
    if inframe and a > cs:
        return "cds_inframe"
    if inframe and a == cs:
        return "cds_canonical"
    if max(a, cs) < min(e, ce):           # overlaps the CDS
        if inframe:                       # same frame, opens 5' of the CDS start
            return "ext_inframe"
        if a < cs:
            return "overlap_uorf"         # RiboCode Overlap_uORF
        if e > ce:
            return "overlap_dorf"         # RiboCode Overlap_dORF
        return "cds_offframe"             # contained and out of frame -> RiboCode `internal`
    if e <= cs:
        return "utr5"
    if a >= ce:
        return "utr3"
    return "other"


def orf_class(crel, biotype):
    """(cds_relationship, transcript biotype) -> pgx ORF class.

    Mirrors RiboCode's ORF_type -> pgx mapping (rc_io.TYPE_MAP) so the model arm and the naively
    enumerated null arms are classified by ONE rule and the databases stay comparable.
    `canonical` and `internal` are not novel and are dropped by the DB builder.
    """
    if crel == "cds_canonical":
        return "canonical"
    if biotype == "lncRNA" or crel == "no_cds":
        return "lncRNA_orf"
    if crel in ("utr5", "overlap_uorf"):
        return "uORF"
    if crel in ("utr3", "dorf_inframe", "overlap_dorf"):
        return "dORF"
    if crel == "ext_inframe":
        return "nterm_ext"
    if crel in ("cds_offframe", "cds_inframe"):
        return "internal"
    return "other"


def bh(pvals):
    """Benjamini-Hochberg q-values (monotone), same procedure RiboCode uses for adjusted_pval.

    NaN-SAFE, and it must be. A single NaN p-value used to null EVERY q-value: numpy sorts NaN
    last, and the reverse `np.minimum.accumulate` therefore starts on that NaN, which propagates
    through the whole array because np.minimum(NaN, x) is NaN. One degenerate test (a zero-variance
    Wilcoxon on a flat extension region) silently zeroed entire arms -- measured: B721/mamba4,
    DoHH2/attn and SU-DHL-4/attn standard extension arms each reported 0 passing from 4,473 / 5,261
    / 7,757 tested, while 4,434 / 5,179 / 7,678 of those would pass the whole-ORF f0 criterion.
    It looks like a real biological result and is arithmetic.

    NaN p-values are excluded from the correction (they are not tests) and get NaN q-values back,
    so they can never pass a q threshold but also never poison their neighbours.
    """
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.size, np.nan, dtype=float)
    if p.size == 0:
        return out
    finite = np.isfinite(p)
    n = int(finite.sum())
    if n == 0:
        return out
    idx = np.flatnonzero(finite)
    pf = p[idx]
    order = np.argsort(pf)
    q = pf[order] * n / (np.arange(n) + 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    qf = np.empty(n, dtype=float)
    qf[order] = np.clip(q, 0.0, 1.0)
    out[idx] = qf
    return out
