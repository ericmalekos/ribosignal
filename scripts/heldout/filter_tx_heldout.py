#!/usr/bin/env python3
"""Parameterized transcriptome-BAM cleaner for held-out Ribo-seq (posture A).

Byte-for-byte the same logic as
expression_context_human/scripts/phase30/filter_transcriptome_bam.py, but with the
ncRNA-tx list and tx->gene map passed as arguments so it works for either species.
Drops records pointing at ncRNA transcripts, then per read: keep all records if they
map within a single gene (isoform multimapping is legitimate), drop the whole read if
it spans multiple genes (cross-gene paralog ambiguity). Sorts + indexes + overwrites
the BAM in place. Relies on STAR's TranscriptomeSAM being grouped by query name.

Usage: filter_tx_heldout.py <bam> <ncrna_tx.txt> <tx_to_gene.tsv>
"""
import os
import shutil
import sys
from pathlib import Path

import pysam


def main(bam_path, ncrna_file, t2g_file):
    bam_in_path = Path(bam_path)
    if not bam_in_path.exists():
        sys.exit(f"missing input BAM: {bam_in_path}")
    ncrna_tx = {ln.strip() for ln in open(ncrna_file) if ln.strip()}
    tx_to_gene = {}
    with open(t2g_file) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                tx_to_gene[p[0]] = p[1]
    print(f"loaded {len(ncrna_tx):,} ncRNA tx + {len(tx_to_gene):,} tx->gene", flush=True)

    scratch = Path(f"/data/tmp/emalekos/heldout_filter/{bam_in_path.stem}_{os.getpid()}")
    scratch.mkdir(parents=True, exist_ok=True)
    out_path = scratch / "clean.bam"

    n_total = n_ncrna = n_xgene = n_uniq = n_wg = 0
    n_reads_kept = n_reads_xgene = 0
    # HARD GUARD: this filter groups records by query name to detect cross-gene reads, so the input
    # MUST be query-grouped (raw STAR TranscriptomeSAM output). On a COORDINATE-SORTED bam the
    # grouping silently collapses -- every record looks like a singleton, cross-gene detection never
    # fires, and the filter reports a plausible ~0.3% drop instead of the true 10-22%. That exact
    # failure produced a wrong "the filter is a no-op" conclusion on 2026-08-12. Fail loudly instead.
    _probe = pysam.AlignmentFile(str(bam_in_path), "rb")
    _hdr = _probe.header.to_dict().get("HD", {})
    if _hdr.get("SO") == "coordinate":
        _probe.close()
        sys.exit(f"REFUSING coordinate-sorted input: {bam_in_path}\n"
                 "  filter_tx_heldout requires query-name-grouped records (raw STAR "
                 "--quantMode TranscriptomeSAM output).\n"
                 "  On a coordinate-sorted bam cross-gene detection silently fails and the drop rate "
                 "is meaningless.\n"
                 "  Run this BEFORE sorting; it sorts and indexes the output itself.")
    _probe.close()

    bam_in = pysam.AlignmentFile(str(bam_in_path), "rb")
    bam_out = pysam.AlignmentFile(str(out_path), "wb", template=bam_in)

    def flush(buf):
        nonlocal n_xgene, n_uniq, n_wg, n_reads_kept, n_reads_xgene
        if not buf:
            return
        genes = {tx_to_gene.get(r.reference_name) for r in buf}
        genes.discard(None)
        if len(genes) <= 1:
            for r in buf:
                bam_out.write(r)
            n_reads_kept += 1
            if len(buf) == 1:
                n_uniq += 1
            else:
                n_wg += 1
        else:
            n_xgene += len(buf)
            n_reads_xgene += 1

    cur = None
    buf = []
    for r in bam_in.fetch(until_eof=True):
        n_total += 1
        if r.is_unmapped:
            continue
        if r.reference_name in ncrna_tx:
            n_ncrna += 1
            continue
        if cur is None or r.query_name == cur:
            buf.append(r)
            cur = r.query_name
        else:
            flush(buf)
            cur = r.query_name
            buf = [r]
    flush(buf)
    bam_in.close()
    bam_out.close()
    print(f"records total={n_total:,} drop_ncRNA={n_ncrna:,} drop_xgene={n_xgene:,} "
          f"kept_uniq={n_uniq:,} kept_within_gene={n_wg:,}; "
          f"reads kept={n_reads_kept:,} dropped_xgene={n_reads_xgene:,}", flush=True)

    sp = scratch / "sorted.bam"
    pysam.sort("-@", "8", "-o", str(sp), str(out_path))
    pysam.index("-@", "8", str(sp))
    shutil.move(str(sp), str(bam_in_path))
    shutil.move(str(sp) + ".bai", str(bam_in_path) + ".bai")
    shutil.rmtree(scratch)
    print(f"DONE -> {bam_in_path}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: filter_tx_heldout.py <bam> <ncrna_tx.txt> <tx_to_gene.tsv>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
