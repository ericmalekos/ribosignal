#!/usr/bin/env python3
"""Readers for RiboCode output, shared by the `pgx` steps.

Only `<outname>_collapsed.txt` is used: `detectORF.write_result` is called with
`report_adjusted_pval=True` for the collapsed file ONLY (detectORF.py:472), so the plain
`<outname>.txt` has no BH-adjusted q-value and cannot be significance-filtered.

Coordinate convention (detectORF.py:354-364), converted here to 0-based half-open:
    ORF_tstart      = orf_iv.start + 1        -> start0     = ORF_tstart - 1
    ORF_tstop       = orf_iv.end + 3          -> end0_excl  = ORF_tstop      (stop codon included)
    annotated_tstart= startcodon.start + 1    -> cds_start0 = annotated_tstart - 1

When alternative start codons are enabled, RiboCode inserts an extra `start_codon` column
(detectORF.py:152). Reading by header name rather than index absorbs that automatically.

RiboCode ORF_type -> pgx class. `annotated` and `internal` are NOT novel: `annotated` is the
canonical protein (or an N-terminal extension of it, see is_nterm_extension) and `internal` is a
CDS truncation. Matches build_ribocode_db.py's SKIP_TYPES / TYPE_MAP.
"""
from __future__ import annotations

import csv

SKIP_TYPES = {"annotated", "internal"}
TYPE_MAP = {"uORF": "uORF", "Overlap_uORF": "uORF", "dORF": "dORF", "Overlap_dORF": "dORF"}


def _f(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _i(v, default=None):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def read_collapsed(path):
    """Yield one dict per called ORF with typed, 0-based-converted fields."""
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tstart, tstop = _i(r.get("ORF_tstart")), _i(r.get("ORF_tstop"))
            ann = _i(r.get("annotated_tstart"))
            yield {
                "orf_id": r.get("ORF_ID", ""),
                "orf_type": r.get("ORF_type", ""),
                "tx": r.get("transcript_id", ""),
                "tx_type": r.get("transcript_type", ""),
                "gene_id": r.get("gene_id", ""),
                "gene_name": r.get("gene_name", ""),
                "gene_type": r.get("gene_type", ""),
                "start0": None if tstart is None else tstart - 1,
                "end0": tstop,
                "cds_start0": None if ann is None else ann - 1,
                "aa_len": _i(r.get("ORF_length"), 0),
                "start_codon": (r.get("start_codon") or "ATG").upper(),
                "psites_f0": _f(r.get("Psites_sum_frame0"), 0.0),
                "pval": _f(r.get("pval_combined"), 1.0),
                "qval": _f(r.get("adjusted_pval"), _f(r.get("pval_combined"), 1.0)),
                "aaseq": (r.get("AAseq") or "").replace("*", "").strip(),
            }


def pgx_class(orf, gene_type=None):
    """RiboCode ORF_type + transcript biotype -> pgx novel class, or None if not novel.

    lncRNA transcripts get `lncRNA_orf` regardless of the positional call, mirroring
    build_a549_dbs / build_ribocode_db so classes are comparable across DB builders.
    """
    t = orf["orf_type"]
    if t in SKIP_TYPES:
        return None
    bt = gene_type if gene_type is not None else orf.get("gene_type", "")
    if t == "novel" or bt == "lncRNA" or orf.get("tx_type") == "lncRNA":
        return "lncRNA_orf"
    return TYPE_MAP.get(t)


def is_nterm_extension(orf):
    """True if this `annotated` call starts UPSTREAM of the annotated CDS start (same stop).

    RiboCode's classfy_orf returns "annotated" whenever orf.end == cds.end, and start_check with
    only_longest_orf=True takes the most 5' in-frame start, so an ATG N-terminal extension is
    already emitted here -- just labelled `annotated`. This is the ONLY route by which RiboCode
    reports an extension: orf_finder.orf_find sets alt_flag=0 as soon as any in-frame ATG exists
    before the stop, so a near-cognate extension of an annotated CDS is unreachable at any
    parameter setting and needs the separate scan in pgx/extensions.py.
    """
    return (orf["orf_type"] == "annotated" and orf["cds_start0"] is not None
            and orf["start0"] is not None and orf["start0"] < orf["cds_start0"])


def load_tx2cds(path):
    """{tx_id: (cds_start0, cds_end0_excl)} for transcripts with an annotated CDS."""
    out = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len")}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            u5, cl = int(f[ix["utr5_len"]]), int(f[ix["cds_len"]])
            if cl > 0:
                out[f[ix["tx_id"]]] = (u5, u5 + cl)
    return out


def load_biotype(path):
    """{tx_id: gene_type}."""
    out = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ti, gi = h.index("tx_id"), h.index("gene_type")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            out[f[ti]] = f[gi]
    return out
