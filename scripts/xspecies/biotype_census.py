#!/usr/bin/env python3
"""Gene and transcript counts by biotype class across all nine species.

Reads the normalized `tx2biotype.tsv` that scripts/normalize_annotation.py writes,
so every species is counted in the SAME (GENCODE) vocabulary. Counting the raw GTFs
instead would compare `mRNA` against `protein_coding` against a bare feature-type
column and be meaningless.

Human and mouse come from the pre-existing GENCODE tables and are included as the
reference points the model was actually trained on.

Category assignment is deliberately explicit, and two choices are worth stating:

  * `rRNA_pseudogene` and `tRNA_pseudogene` are counted under dropped_ncRNA, not
    PSEUDOGENE. They are pseudogenes by name but they are in the ncRNA drop set,
    which is what matters here; putting them in both buckets would double-count and
    putting them under pseudogene would make the drop column understate the drop set.
  * dropped_ncRNA is exactly the extended drop set (INCLUDING miRNA), so its
    transcript total equals `ncrna_tx.txt` for that species. That is asserted, not
    assumed, and the assertion caught miRNA being filed in the wrong bucket.

Usage: python3 scripts/xspecies/biotype_census.py [--tsv <out.tsv>]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")

# species -> tx2biotype path. Human and mouse predate data/annot/ and keep their own.
SOURCES = [
    ("human",     NEW / "data" / "tx2biotype.tsv",        NEW / "data" / "ncrna_tx_human_v49.txt"),
    ("mouse",     NEW / "data" / "tx2biotype_mouse.tsv",  NEW / "data" / "ncrna_tx_mouse_vM38.txt"),
    ("gorilla",   NEW / "data" / "annot" / "gorilla"   / "tx2biotype.tsv", None),
    ("chimp",     NEW / "data" / "annot" / "chimp"     / "tx2biotype.tsv", None),
    ("macaque",   NEW / "data" / "annot" / "macaque"   / "tx2biotype.tsv", None),
    ("zebrafish", NEW / "data" / "annot" / "zebrafish" / "tx2biotype.tsv", None),
    ("celegans",  NEW / "data" / "annot" / "celegans"  / "tx2biotype.tsv", None),
    ("fly",       NEW / "data" / "annot" / "fly"       / "tx2biotype.tsv", None),
    ("yeast",     NEW / "data" / "annot" / "yeast"     / "tx2biotype.tsv", None),
    # Same assemblies as the GENCODE human/mouse above, annotated by RefSeq instead.
    # These exist so the nine-species comparison can be made in ONE annotation source.
    ("human_refseq", NEW / "data" / "annot" / "human_refseq" / "tx2biotype.tsv", None),
    ("mouse_refseq", NEW / "data" / "annot" / "mouse_refseq" / "tx2biotype.tsv", None),
]

# EXACTLY the extended drop set from normalize_annotation.py. miRNA belongs here and
# not under small_ncRNA: it is dropped before any ORF caller sees the BAM, and that is
# what this column is for. Leaving it out was the first version of this script, and the
# match assertion below caught it -- human came out 568 against a 2,447-line drop list,
# the difference being exactly the 1,879 miRNA transcripts.
CONTAMINANT = {"rRNA", "Mt_rRNA", "rRNA_pseudogene", "miRNA", "tRNA", "Mt_tRNA",
               "tRNA_pseudogene", "SRP_RNA", "RNase_P_RNA", "RNase_MRP_RNA"}
SMALL_NC = {"snRNA", "snoRNA", "scaRNA", "misc_RNA", "sRNA", "scRNA",
            "vault_RNA", "ribozyme", "piRNA"}
LNC = {"lncRNA", "lnc_RNA", "processed_transcript", "TEC"}
CODING = {"protein_coding", "nonsense_mediated_decay", "non_stop_decay",
          "protein_coding_CDS_not_defined", "protein_coding_LoF", "retained_intron"}

ORDER = ["protein_coding", "lncRNA", "small_ncRNA", "dropped_ncRNA", "pseudogene", "other"]


def classify(bt: str) -> str:
    if bt in CONTAMINANT:
        return "dropped_ncRNA"
    if bt in CODING:
        return "protein_coding"
    if bt in LNC:
        return "lncRNA"
    if bt in SMALL_NC:
        return "small_ncRNA"
    if "pseudogene" in bt:
        return "pseudogene"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", type=Path,
                    default=NEW / "results" / "xspecies_biotype_census.tsv")
    args = ap.parse_args()

    rows = []
    for sp, path, ncfile in SOURCES:
        if not path.exists():
            print(f"  {sp}: MISSING {path}", file=sys.stderr)
            continue
        tx = defaultdict(int)
        genes: dict[str, set] = defaultdict(set)
        other_terms: dict[str, int] = defaultdict(int)
        n_tx = 0
        with path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n_tx += 1
                cls = classify(r["transcript_type"])
                tx[cls] += 1
                genes[classify(r["gene_type"])].add(r["gene_id"])
                if cls == "other":
                    other_terms[r["transcript_type"]] += 1
        rec = {"species": sp, "n_tx": n_tx,
               "n_genes": len({g for s in genes.values() for g in s})}
        for c in ORDER:
            rec[f"tx_{c}"] = tx.get(c, 0)
            rec[f"gene_{c}"] = len(genes.get(c, set()))
        rec["other_terms"] = ", ".join(f"{k} {v}" for k, v in
                                       sorted(other_terms.items(), key=lambda x: -x[1])[:4])
        # The dropped_ncRNA class must equal the ncRNA drop list exactly.
        if ncfile and ncfile.exists():
            n_drop = sum(1 for _ in ncfile.open())
        else:
            f = path.parent / "ncrna_tx.txt"
            n_drop = sum(1 for _ in f.open()) if f.exists() else -1
        rec["ncrna_tx_file"] = n_drop
        rec["drop_matches"] = "yes" if n_drop == rec["tx_dropped_ncRNA"] else f"NO ({n_drop})"
        rows.append(rec)

    w = 10
    print("\nTRANSCRIPTS by class\n")
    hdr = f"  {'species':<10}{'total':>9}" + "".join(f"{c.replace('protein_','p_'):>14}" for c in ORDER)
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['species']:<10}{r['n_tx']:>9,}" +
              "".join(f"{r['tx_'+c]:>14,}" for c in ORDER))

    print("\nGENES by class\n")
    hdr = f"  {'species':<10}{'total':>9}" + "".join(f"{c.replace('protein_','p_'):>14}" for c in ORDER)
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['species']:<10}{r['n_genes']:>9,}" +
              "".join(f"{r['gene_'+c]:>14,}" for c in ORDER))

    print("\ndropped_ncRNA class vs the ncRNA drop list actually written\n")
    for r in rows:
        print(f"  {r['species']:<10} dropped_ncRNA tx {r['tx_dropped_ncRNA']:>7,}   "
              f"ncrna_tx.txt {r['ncrna_tx_file']:>7,}   match: {r['drop_matches']}")

    print("\n'other' (no GENCODE equivalent, deliberately passed through)\n")
    for r in rows:
        if r["other_terms"]:
            print(f"  {r['species']:<10} {r['other_terms']}")

    args.tsv.parent.mkdir(parents=True, exist_ok=True)
    cols = ["species", "n_genes", "n_tx"] + \
           [f"gene_{c}" for c in ORDER] + [f"tx_{c}" for c in ORDER] + \
           ["ncrna_tx_file", "drop_matches", "other_terms"]
    with args.tsv.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                            extrasaction="ignore", lineterminator="\n")
        wr.writeheader(); wr.writerows(rows)
    print(f"\nwrote {args.tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
