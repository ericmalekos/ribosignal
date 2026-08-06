#!/usr/bin/env python3
"""Clean a per-sample transcriptome-coord Ribo-seq BAM, identical logic to the training pipeline's
expression_context_human/scripts/phase30/filter_transcriptome_bam.py, but with paths as args so it runs on
the HBL-1 Ribo-seq. STAR --quantMode TranscriptomeSAM writes one record per (read, transcript) hit:
  1. drop records on ncRNA (rRNA/tRNA/miRNA/snoRNA/...) transcripts (ncrna_filter_tx.txt).
  2. group surviving records by read name: keep if all map to the SAME gene_id (within-gene isoform
     multimapping = posture A, retained + each counted); drop the read entirely if it spans >1 gene_id
     (cross-gene paralog ambiguity). STAR's transcriptome BAM is grouped by query name, so stream it.
  3. write cleaned BAM, sort, index, overwrite input. riboseq/ribocode env (pysam).

Usage: filter_tx_bam_generic.py <bam> <ncrna_filter_tx.txt> <tx_to_gene.tsv>"""
import os
import shutil
import sys
from pathlib import Path

import pysam


def main(bam_in_path, ncrna_path, tx2gene_path):
    bam_in_path = Path(bam_in_path)
    if not bam_in_path.exists():
        sys.exit(f"missing input BAM: {bam_in_path}")
    ncrna_tx = {ln.strip() for ln in open(ncrna_path) if ln.strip()}
    tx_to_gene = {}
    with open(tx2gene_path) as f:
        for ln in f:
            a, b = ln.rstrip("\n").split("\t")
            tx_to_gene[a] = b
    print(f"[{bam_in_path.name}] {len(ncrna_tx):,} ncRNA tx + {len(tx_to_gene):,} tx->gene", flush=True)

    scratch = Path(f"/data/tmp/emalekos/rc_clean/{bam_in_path.stem}_{os.getpid()}")
    scratch.mkdir(parents=True, exist_ok=True)
    out_path = scratch / "clean.bam"

    n_total = n_dropped_ncrna = n_dropped_xgene = n_unique = n_within_gene = 0
    n_reads_kept = n_reads_dropped_xgene = 0
    bam_in = pysam.AlignmentFile(str(bam_in_path), "rb")
    bam_out = pysam.AlignmentFile(str(out_path), "wb", template=bam_in)

    def flush_read(buf):
        nonlocal n_dropped_xgene, n_unique, n_within_gene, n_reads_kept, n_reads_dropped_xgene
        if not buf:
            return
        genes = {tx_to_gene.get(rec.reference_name) for rec in buf}
        genes.discard(None)
        if len(genes) <= 1:
            for rec in buf:
                bam_out.write(rec)
            n_reads_kept += 1
            if len(buf) == 1:
                n_unique += 1
            else:
                n_within_gene += 1
        else:
            n_dropped_xgene += len(buf)
            n_reads_dropped_xgene += 1

    cur_query, buf = None, []
    for rec in bam_in.fetch(until_eof=True):
        n_total += 1
        if rec.is_unmapped:
            continue
        if rec.reference_name in ncrna_tx:
            n_dropped_ncrna += 1
            continue
        if cur_query is None or rec.query_name == cur_query:
            buf.append(rec); cur_query = rec.query_name
        else:
            flush_read(buf); cur_query = rec.query_name; buf = [rec]
    flush_read(buf)
    bam_in.close(); bam_out.close()

    print(f"[{bam_in_path.name}] records: total={n_total:,} drop_ncRNA={n_dropped_ncrna:,} "
          f"drop_xgene={n_dropped_xgene:,} kept_unique={n_unique:,} kept_within_gene={n_within_gene:,}",
          flush=True)
    print(f"[{bam_in_path.name}] reads: kept={n_reads_kept:,} dropped_xgene={n_reads_dropped_xgene:,}",
          flush=True)

    sorted_path = scratch / "sorted.bam"
    pysam.sort("-@", "8", "-o", str(sorted_path), str(out_path))
    pysam.index("-@", "8", str(sorted_path))
    shutil.move(str(sorted_path), str(bam_in_path))
    shutil.move(str(sorted_path) + ".bai", str(bam_in_path) + ".bai")
    shutil.rmtree(scratch)
    print(f"[{bam_in_path.name}] DONE -> {bam_in_path}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: filter_tx_bam_generic.py <bam> <ncrna_filter_tx.txt> <tx_to_gene.tsv>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
