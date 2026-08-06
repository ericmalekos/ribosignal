#!/usr/bin/env python3
"""Build db_ribocode.fasta from RiboCode's OWN ORF calls on the HBL-1 DMSO Ribo-seq, wrapped in the SAME
canonical GENCODE v49 base + REV_ decoy scheme as db_{canonical,model,null,price}. This is the PRIMARY
experimental baseline (replacing the PRICE processed output): RiboCode is the same tool lineage as the
model's training P-sites, so model-vs-RiboCode isolates predicted-vs-measured translation rather than
caller discordance. RiboCode ORF_type -> novel classes: annotated = canonical (SKIP); internal = CDS
truncation (SKIP, matches build_a549_dbs excluding 'internal'); uORF/Overlap_uORF -> uORF;
dORF/Overlap_dORF -> dORF; novel -> lncRNA_orf if on a lncRNA transcript else novel. Novel proteins are
deduped by AAseq, dropped if == canonical, tagged nuORF|ribocode|<ORF_ID> for class-specific FDR. cas12a env."""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
SKIP_TYPES = {"annotated", "internal"}            # canonical + CDS-truncation, excluded (parallels db_model)
TYPE_MAP = {"uORF": "uORF", "Overlap_uORF": "uORF", "dORF": "dORF", "Overlap_dORF": "dORF"}


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
    ap.add_argument("--collapsed", default=str(NEW / "proteogenomics/data/HBL1_pilot/ribocode/HBL1_DMSO_collapsed.txt"))
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

    novel = {}                        # seq -> (orf_id, klass)
    n_read = skipped = same_canon = short = 0
    with open(args.collapsed) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            n_read += 1
            otype = r["ORF_type"]
            if otype in SKIP_TYPES:
                skipped += 1; continue
            seq = (r.get("AAseq") or "").replace("*", "").strip()
            if len(seq) < args.min_aa:
                short += 1; continue
            if seq in canon_seqs:
                same_canon += 1; continue
            if otype in TYPE_MAP:
                klass = TYPE_MAP[otype]
            elif otype == "novel":
                klass = "lncRNA_orf" if (r.get("transcript_type") == "lncRNA") else "novel"
            else:
                klass = otype
            oid = r["ORF_ID"]
            if seq not in novel:
                novel[seq] = (oid, klass)
    print(f"RiboCode ORFs read: {n_read:,}; skipped(annotated/internal): {skipped:,}; "
          f"== canonical: {same_canon:,}; < {args.min_aa} aa: {short:,}; unique novel: {len(novel):,}")
    print("  novel by class:", dict(Counter(v[1] for v in novel.values())))

    cm = open(out / "ribocode_class_map.tsv", "w"); cm.write("accession\tclass\tpred_frame0\n")
    with open(out / "db_ribocode.fasta", "w") as fh:
        for seq, acc in canon.items():
            fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
        for seq, (oid, klass) in novel.items():
            acc = f"nuORF|ribocode|{oid.replace('|', '_')}"
            fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
            cm.write(f"{acc}\t{klass}\t1.0000\n")
    cm.close()
    n = len(canon) + len(novel)
    print(f"wrote {out}/db_ribocode.fasta: {n:,} target + {n:,} decoy = {2*n:,} seqs ({len(novel):,} novel)")
    print(f"wrote {out}/ribocode_class_map.tsv")


if __name__ == "__main__":
    main()
