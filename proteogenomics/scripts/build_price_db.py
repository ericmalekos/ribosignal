#!/usr/bin/env python3
"""Build db_price.fasta: the authors' PRICE Ribo-seq ORF calls wrapped in the SAME canonical
GENCODE v49 base + REV_ decoy scheme as db_{canonical,model,null}, so the immunopeptidome search
isolates ORF-SET provenance (model-selected vs 3-frame-null vs authors' experimental PRICE calls).
PRICE novel ORFs are tagged nuORF|price|<id> so compare_model_vs_null's class-specific FDR treats
them as the novel class (REV_nuORF|price| decoys as the honest per-class denominator). cas12a env."""
from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def iter_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = ln[1:].rstrip("\n"); buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--price", default=str(NEW / "proteogenomics/data/HBL1_pilot/ref/HBL1_Price_orfs.fasta"))
    ap.add_argument("--canonical", default="/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
                                           "annotations/gencode_proteins/gencode.v49.pc_translations.fa.gz")
    ap.add_argument("--min_aa", type=int, default=7)
    ap.add_argument("--out_dir", default=str(NEW / "proteogenomics/data/HBL1_pilot/db"))
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    canon = {}
    for name, seq in iter_fasta(args.canonical):
        seq = seq.replace("*", "").strip()
        if len(seq) >= args.min_aa and seq not in canon:
            canon[seq] = name.split("|")[0].split()[0]
    canon_seqs = set(canon)
    print(f"canonical proteins (unique, >= {args.min_aa} aa): {len(canon):,}")

    novel = {}                       # seq -> price_id
    n_read = same_canon = short = 0
    for name, seq in iter_fasta(args.price):
        n_read += 1
        seq = seq.replace("*", "").strip()
        if len(seq) < args.min_aa:
            short += 1; continue
        if seq in canon_seqs:
            same_canon += 1; continue
        pid = name.split()[0].replace("|", "_")
        if seq not in novel:
            novel[seq] = pid
    print(f"PRICE ORFs read: {n_read:,}; == canonical (dropped): {same_canon:,}; "
          f"< {args.min_aa} aa: {short:,}; unique novel: {len(novel):,}")

    cm = open(out / "price_class_map.tsv", "w"); cm.write("accession\tclass\tpred_frame0\n")
    with open(out / "db_price.fasta", "w") as fh:
        for seq, acc in canon.items():
            fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
        for seq, pid in novel.items():
            acc = f"nuORF|price|{pid}"
            fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
            cm.write(f"{acc}\tprice\t1.0000\n")
    cm.close()
    n = len(canon) + len(novel)
    print(f"wrote {out}/db_price.fasta: {n:,} target + {n:,} decoy = {2*n:,} seqs "
          f"({len(novel):,} novel PRICE)")
    print(f"wrote {out}/price_class_map.tsv")


if __name__ == "__main__":
    main()
