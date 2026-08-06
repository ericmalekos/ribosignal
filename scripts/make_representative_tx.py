#!/usr/bin/env python3
"""Posture-B (multimap sensitivity check, results.md Task 21): pick the highest-expressed isoform per
gene, among the packed Fibroblast universe, from the matched Fibroblast RNA-seq salmon quant.

For the representative (highest-TPM) isoform of a gene, its posture-A per-nt Ribo-seq counts already
EQUAL its posture-B (one-isoform-commit) counts -- a footprint on that isoform lands in an exon it
contains, so committing it there leaves the representative's counts unchanged; only low-expressed
siblings lose inherited signal, and those drop below min_signal anyway. So restricting the universe to
representatives yields the posture-B target with NO target re-derivation: dataset.representative_tx()
reads this list (via env RIBO_REPRESENTATIVE_TX) and every split keeps only these tx.

Tie-break for equal TPM: longest transcript, then lexicographic tx_id (deterministic). Genes with no
packed isoform contribute nothing. Writes data/fibroblast_representative_tx.txt (one versioned tx_id per
line) + prints a composition report. cas12a env (stdlib only)."""
from __future__ import annotations

import argparse
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpm", default=str(NEW / "data/fibroblast_salmon_mean_tpm.tsv"))
    ap.add_argument("--tx2biotype", default=str(NEW / "data/tx2biotype.tsv"))
    ap.add_argument("--universe", default=str(NEW / "data/packed/tx_order.txt"),
                    help="restrict to isoforms present in this packed universe")
    ap.add_argument("--out", default=str(NEW / "data/fibroblast_representative_tx.txt"))
    args = ap.parse_args()

    packed = {ln.split()[0] for ln in Path(args.universe).read_text().splitlines() if ln.strip()}

    tpm = {}
    with open(args.tpm) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tpm[f[ci["tx_id"]]] = float(f[ci["mean_tpm"]])

    # tx -> (gene_id, gene_type, length)
    meta = {}
    with open(args.tx2biotype) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tx = f[ci["tx_id"]]
            if tx in packed:
                meta[tx] = (f[ci["gene_id"]], f[ci["gene_type"]], int(f[ci["length"]]))

    # group packed tx by gene, pick argmax TPM (tie: longest, then tx_id)
    by_gene = {}
    for tx, (gid, gtype, length) in meta.items():
        by_gene.setdefault(gid, []).append(tx)
    reps = []
    for gid, txs in by_gene.items():
        best = max(txs, key=lambda t: (tpm.get(t, 0.0), meta[t][2], t))
        reps.append(best)
    reps.sort()

    Path(args.out).write_text("\n".join(reps) + "\n")

    # report
    n_packed = len(packed)
    n_meta = len(meta)
    n_genes = len(by_gene)
    from collections import Counter
    bt = Counter(meta[t][1] for t in reps)
    multi = sum(1 for txs in by_gene.values() if len(txs) > 1)
    print(f"packed universe tx           : {n_packed:,}")
    print(f"  with gene mapping (kept)   : {n_meta:,}  ({n_packed - n_meta:,} unmapped, dropped)")
    print(f"genes in packed universe     : {n_genes:,}  ({multi:,} multi-isoform)")
    print(f"representatives written      : {len(reps):,}  -> {args.out}")
    print(f"  = {100*len(reps)/n_packed:.1f}% of packed tx kept "
          f"({n_packed - len(reps):,} sibling isoforms dropped)")
    print("representative biotype breakdown:")
    for b, n in bt.most_common(8):
        print(f"  {b:<22} {n:,}")


if __name__ == "__main__":
    main()
