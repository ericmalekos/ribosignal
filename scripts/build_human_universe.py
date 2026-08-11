#!/usr/bin/env python3
"""Build an expression-restricted human transcript universe from salmon quants.

DEFAULT: one universe per dataset, built from that dataset's own RNA-seq. Borrowing another
dataset's quant imports its expression bias, which is how the Soellner artefact happened.

Borrowing is ALLOWED but must be declared with --borrowed-from, which stamps the donor into the
provenance and the TSV header so no downstream table can quote the numbers without the caveat
travelling with them. Use it only when the recipient has no RNA-seq of its own and the biological
system matches (`feedback_never_reuse_universe_across_datasets`, downgraded to a warning
2026-08-10).

Universe = protein_coding + lncRNA transcripts, mature length <= --max-len, with mean salmon TPM
across the dataset's RNA samples >= --min-tpm.

Two exclusions are applied by default and both are project rules rather than taste:

  * chrM. The 13 mitochondrial protein-coding mRNAs are dropped project-wide -- mitoribosome, a
    different genetic code, and leaderless transcripts make them a different problem
    (`feedback_riboseq_exclude_mito_genes`).
  * --exclude-genes. A file of gene IDs to drop. Needed for GSE304796 (CAR-T), whose RNA is
    Ribo-Zero TOTAL RNA: 43.1% of its TPM sits in `ncRNA_host` genes against 15.3% for the poly(A)
    control. Those hosts are typed lncRNA so no biotype filter reaches them, and the reads are NOT
    snoRNA reads (read-level snoRNA depletion moved the number by 0.0 points), so the exclusion has
    to happen here, at quantification, rather than upstream.

Writes <out>_tx.txt (one versioned ENST per line), <out>_universe.fa (bare-header, one line per
sequence) and <out>_universe.tsv (tx_id, length, biotype, mean_tpm).
"""
import argparse
import pathlib
import statistics
import sys

ANN = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations")


def iter_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = ln[1:].split("|")[0].split()[0]
                buf = []
            else:
                buf.append(ln.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salmon", nargs="+", required=True, help="quant.sf files (one per RNA sample)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--tx2biotype", required=True)
    ap.add_argument("--fastas", nargs="+",
                    default=[str(ANN / "gencode.v49.pc_transcripts.fa"),
                             str(ANN / "gencode.v49.lncRNA_transcripts.fa")])
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--max-len", type=int, default=10000)
    ap.add_argument("--exclude-genes", help="file of gene IDs (version-insensitive) to drop")
    ap.add_argument("--borrowed-from",
                    help="donor dataset name when --salmon is NOT this dataset's own RNA; recorded "
                         "in the outputs so the caveat cannot be lost downstream")
    ap.add_argument("--keep-chrm", action="store_true",
                    help="override the project-wide chrM exclusion (do not use without a reason)")
    a = ap.parse_args()

    # mean TPM across samples, keyed by versioned tx id
    per_tx = {}
    for q in a.salmon:
        p = pathlib.Path(q)
        if not p.exists():
            print(f"missing {p}", file=sys.stderr)
            return 1
        with open(p) as fh:
            fh.readline()
            for line in fh:
                f = line.split("\t")
                per_tx.setdefault(f[0].split("|")[0], []).append(float(f[3]))
    mean_tpm = {t: statistics.mean(v) for t, v in per_tx.items()}
    print(f"  salmon: {len(a.salmon)} sample(s), {len(mean_tpm):,} transcripts quantified")

    # biotype / chrom
    bt = {}
    with open(a.tx2biotype) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ti, bi, ci, gi = (hdr.index(x) for x in ("tx_id", "transcript_type", "chrom", "gene_id"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ti, bi, ci, gi):
                bt[f[ti]] = (f[bi], f[ci], f[gi])

    drop_genes = set()
    if a.exclude_genes:
        drop_genes = {g.split(".")[0] for g in open(a.exclude_genes).read().split() if g}
        print(f"  excluding {len(drop_genes):,} gene ids from {a.exclude_genes}")

    keep, rows, n_chrm, n_gene, n_tpm, n_bt = set(), [], 0, 0, 0, 0
    for tx, tpm in mean_tpm.items():
        info = bt.get(tx)
        if info is None:
            continue
        biotype, chrom, gene = info
        if biotype not in ("protein_coding", "lncRNA"):
            n_bt += 1
            continue
        if not a.keep_chrm and chrom in ("chrM", "MT", "chrMT"):
            n_chrm += 1
            continue
        if gene.split(".")[0] in drop_genes:
            n_gene += 1
            continue
        if tpm < a.min_tpm:
            n_tpm += 1
            continue
        keep.add(tx)
        rows.append((tx, biotype, tpm))
    print(f"  filtered out: biotype {n_bt:,} | chrM {n_chrm:,} | excluded genes {n_gene:,} | "
          f"TPM<{a.min_tpm} {n_tpm:,}")

    # sequences, with the length cap applied here so the FASTA and the list cannot disagree
    out = pathlib.Path(a.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    written, lens = set(), {}
    with open(f"{out}_universe.fa", "w") as fh:
        for f in a.fastas:
            for tx, seq in iter_fasta(f):
                if tx in keep and tx not in written and len(seq) <= a.max_len:
                    fh.write(f">{tx}\n{seq}\n")
                    written.add(tx)
                    lens[tx] = len(seq)
    n_long = len(keep) - len(written)
    print(f"  dropped {n_long:,} over --max-len {a.max_len} or absent from the source FASTAs")

    with open(f"{out}_tx.txt", "w") as fh:
        for tx in sorted(written):
            fh.write(tx + "\n")
    with open(f"{out}_universe.tsv", "w") as fh:
        if a.borrowed_from:
            fh.write(f"# BORROWED UNIVERSE: expression comes from {a.borrowed_from}, not from this "
                     f"dataset's own RNA-seq. Label this in every table that uses it.\n")
        fh.write("tx_id\tlength\tbiotype\tmean_tpm\n")
        for tx, biotype, tpm in sorted(rows):
            if tx in written:
                fh.write(f"{tx}\t{lens[tx]}\t{biotype}\t{tpm:.4f}\n")

    pc = sum(1 for t, b, _ in rows if t in written and b == "protein_coding")
    if a.borrowed_from:
        print(f"  WARNING: borrowed universe -- expression from {a.borrowed_from}, "
              f"not this dataset's own RNA")
    print(f"  UNIVERSE: {len(written):,} tx  ({pc:,} protein_coding, {len(written)-pc:,} lncRNA)")
    print(f"  -> {out}_tx.txt / _universe.fa / _universe.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
