#!/usr/bin/env python3
"""Held-out-chromosome, gene-disjoint splits for the Fibroblast universe.

Grouping transcripts by chromosome makes the split gene-disjoint automatically (a gene
lives on one chromosome, so holding out a chromosome holds out all its genes and all
their isoforms). Chromosomes are assigned to K folds greedily to balance transcript
count, so each fold is a chromosome group with roughly equal size.

Input:  data/fibroblast_universe.tsv (tx_id, ..., chrom, ...)
Output: data/splits/fibroblast_chrom_kfold.json
          {scheme, n_folds, chrom_to_fold, folds:[{fold, test_chroms, n_train_tx,
           n_test_tx, train:[tx...], test:[tx...]}]}
        Fold 0 is also the canonical single held-out-chromosome split for the first check.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
UNIV = NEW / "data" / "fibroblast_universe.tsv"
OUT = NEW / "data" / "splits" / "fibroblast_chrom_kfold.json"
K = 5


def main():
    by_chrom = defaultdict(list)
    with UNIV.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = header.index("chrom")
        ti = header.index("tx_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            by_chrom[f[ci]].append(f[ti])

    # greedy: assign chromosomes (largest first) to the currently-smallest fold
    chroms = sorted(by_chrom, key=lambda c: -len(by_chrom[c]))
    fold_tx = [[] for _ in range(K)]
    fold_chroms = [[] for _ in range(K)]
    chrom_to_fold = {}
    for c in chroms:
        j = min(range(K), key=lambda k: len(fold_tx[k]))
        fold_tx[j].extend(by_chrom[c])
        fold_chroms[j].append(c)
        chrom_to_fold[c] = j

    all_tx = [t for c in by_chrom for t in by_chrom[c]]
    all_set = set(all_tx)
    folds = []
    for k in range(K):
        test = set(fold_tx[k])
        train = [t for t in all_tx if t not in test]
        folds.append({
            "fold": k,
            "test_chroms": sorted(fold_chroms[k]),
            "n_train_tx": len(train),
            "n_test_tx": len(test),
            "train": train,
            "test": sorted(test),
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as o:
        json.dump({
            "scheme": "chrom_kfold_gene_disjoint",
            "n_folds": K,
            "n_total_tx": len(all_set),
            "chrom_to_fold": chrom_to_fold,
            "folds": folds,
        }, o)

    print(f"universe tx: {len(all_set):,}  chromosomes: {len(by_chrom)}", file=sys.stderr)
    for k in range(K):
        print(f"  fold {k}: {folds[k]['n_test_tx']:>6,} test tx  "
              f"chroms={','.join(folds[k]['test_chroms'])}", file=sys.stderr)
    print(f"wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
