#!/usr/bin/env python3
"""Build the three A549 proteogenomics search DBs from the model-scored candidate ORFs (candidates.faa)
plus the canonical GENCODE v49 proteome. All three share the identical canonical base + REV_ decoys; they
differ ONLY in the novel-ORF set, so the search isolates the model's SELECTION value:

  db_canonical : canonical proteome only                         (baseline)
  db_model     : canonical + model-SELECTED novel ORFs           (pred_frame0 >= --thresh)
  db_null      : canonical + ALL novel-ORF candidates            (3-frame-ish naive background)

START-CODON ASYMMETRY (2026-08-01). --null_starts restricts the NULL to ORFs opening at those codons
(default ATG), while the MODEL may select any start present in candidates.faa. This is deliberate and
is the point of the ncStart column: the null represents what a naive pipeline builds -- every AUG..stop
ORF in the expressed transcriptome -- so a CUG/GUG/ACG-initiated ORF that the model selects yields
peptides that exist in NO AUG-only database, at any threshold. Letting non-AUG ORFs into the null would
both destroy that contrast and inflate the null ~10x. State the asymmetry explicitly in any writeup.

Novel classes = uORF / dORF / lncRNA_orf / nterm_ext (internal in-frame = CDS truncations, excluded).
Proteins are deduplicated by sequence and dropped if identical to a canonical protein (not novel). Decoys
are inline reversed sequences with a REV_ accession prefix (MSFragger/Comet decoy_prefix=REV_). A class map
(accession -> class) is written for class-specific FDR. cas12a env (stdlib + gzip)."""
from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
NOVEL_CLASSES = {"uORF", "dORF", "lncRNA_orf", "nterm_ext"}


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
    ap.add_argument("--candidates", default=str(NEW / "proteogenomics/data/A549_pilot/orfs/candidates.faa"))
    ap.add_argument("--canonical", default="/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
                                           "annotations/gencode_proteins/gencode.v49.pc_translations.fa.gz")
    ap.add_argument("--null_starts", default="ATG",
                    help="start codons the NULL db may contain (default ATG = the naive baseline a lab "
                         "would actually build). The MODEL db is unrestricted, so a model-selected "
                         "non-AUG ORF yields peptides present in NO AUG-only database at any threshold. "
                         "See the start-codon asymmetry note in the module docstring.")
    ap.add_argument("--thresh", type=float, default=0.5, help="pred_frame0 cut for model-selected novel")
    ap.add_argument("--min_aa", type=int, default=7, help="min protein length (tryptic detectability)")
    ap.add_argument("--out_dir", default=str(NEW / "proteogenomics/data/A549_pilot/db"))
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    # ---- canonical proteome: dedup by sequence, >= min_aa ----
    canon = {}                       # seq -> accession
    for name, seq in iter_fasta(args.canonical):
        seq = seq.replace("*", "").strip()
        if len(seq) >= args.min_aa and seq not in canon:
            canon[seq] = name.split("|")[0].split()[0]
    canon_seqs = set(canon)
    print(f"canonical proteins (unique, >= {args.min_aa} aa): {len(canon):,}")

    # ---- novel candidates: parse class + f0 from candidates.faa, dedup, drop == canonical ----
    null_starts = {c.strip().upper() for c in args.null_starts.split(",") if c.strip()}
    novel_all = {}                   # seq -> (accession, class, best_f0, start) -- NULL db
    novel_sel = {}                   # model-selected subset                     -- MODEL db
    start_all = Counter(); start_sel = Counter()
    n_read = 0
    for name, seq in iter_fasta(args.candidates):
        n_read += 1
        seq = seq.replace("*", "").strip()
        if len(seq) < args.min_aa or seq in canon_seqs:
            continue
        parts = name.split()
        oid = parts[0]               # tx|a|e|class
        klass = oid.split("|")[-1]
        if klass not in NOVEL_CLASSES:
            continue
        f0 = 0.0
        start = "ATG"                # pre-2026-08-01 candidates.faa carry no start= token
        for tok in parts[1:]:
            if tok.startswith("f0="):
                f0 = float(tok[3:])
            elif tok.startswith("start="):
                start = tok[6:].upper()
        # NULL: naive baseline, restricted to --null_starts (default ATG only)
        if start in null_starts:
            prev = novel_all.get(seq)
            if prev is None or f0 > prev[2]:
                novel_all[seq] = (oid, klass, f0, start)
                start_all[start] += 1
        # MODEL: any start, gated on the score
        if f0 >= args.thresh:
            prev = novel_sel.get(seq)
            if prev is None or f0 > prev[2]:
                novel_sel[seq] = (oid, klass, f0, start)
                start_sel[start] += 1
    print(f"candidate ORFs read: {n_read:,}; unique novel proteins (all): {len(novel_all):,}; "
          f"model-selected (f0>={args.thresh}): {len(novel_sel):,}")
    print("  all-novel by class:", dict(Counter(v[1] for v in novel_all.values())))
    print("  selected  by class:", dict(Counter(v[1] for v in novel_sel.values())))
    print(f"  NULL  start codons (restricted to {sorted(null_starts)}):",
          dict(Counter(v[3] for v in novel_all.values()).most_common()))
    print("  MODEL start codons:", dict(Counter(v[3] for v in novel_sel.values()).most_common()))
    n_nonaug = sum(1 for v in novel_sel.values() if v[3] != "ATG")
    print(f"  MODEL non-AUG selected: {n_nonaug:,} / {len(novel_sel):,} "
          f"({n_nonaug / max(len(novel_sel), 1) * 100:.1f}%) -- these can yield ncStart peptides")

    def acc_for(seq, info):
        # unique, class-encoding accession; keep short + FDR-parseable
        oid, klass, f0 = info[0], info[1], info[2]
        return f"nuORF|{klass}|{oid.replace('|', '_')}"

    def write_db(path, novel_dict, class_map_fh):
        with open(path, "w") as fh:
            for seq, acc in canon.items():
                fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
            for seq, info in novel_dict.items():
                acc = acc_for(seq, info)
                fh.write(f">{acc}\n{seq}\n>REV_{acc}\n{seq[::-1]}\n")
                if class_map_fh is not None:
                    class_map_fh.write(f"{acc}\t{info[1]}\t{info[2]:.4f}\t{info[3]}\n")
        n = len(canon) + len(novel_dict)
        print(f"  wrote {path.name}: {n:,} target + {n:,} decoy = {2*n:,} seqs")

    write_db(out / "db_canonical.fasta", {}, None)
    cm_model = open(out / "model_class_map.tsv", "w"); cm_model.write("accession\tclass\tpred_frame0\tstart_codon\n")
    write_db(out / "db_model.fasta", novel_sel, cm_model); cm_model.close()
    cm_null = open(out / "null_class_map.tsv", "w"); cm_null.write("accession\tclass\tpred_frame0\tstart_codon\n")
    write_db(out / "db_null.fasta", novel_all, cm_null); cm_null.close()
    print(f"\nDBs + class maps in {out}")


if __name__ == "__main__":
    main()
