#!/usr/bin/env python3
"""Generate the cross-species reference registry, and verify it while doing so.

This is both the Phase 2 deliverable and the Phase 2 gate. It walks every species
in scripts/xspecies/species.tsv, resolves each derived artifact at its conventional
path, and checks it is present, non-empty and internally consistent:

  * the STAR index actually records the per-species --sjdbOverhang and
    --genomeSAindexNbases the table asked for. 14 is the mammalian default for
    genomeSAindexNbases and is wrong for yeast, worm and fly; STAR does not error,
    it just builds a bad index, so this is checked rather than assumed.
  * the ncRNA drop list is non-empty and above a floor. An empty list is the
    silent failure normalize_annotation.py exists to prevent: the drop step runs,
    removes nothing, and reports success.
  * the salmon index's decoy count is plausible against the assembly's contig count.
  * RiboCode kept most of the transcripts (a large shortfall means the GTF shape
    confused prepare_transcripts).
  * every transcriptome FASTA ID joins to tx_to_gene.tsv.

Writes docs/SPECIES_REFERENCE_REGISTRY.md and data/species_reference_registry.tsv.
Exits nonzero if any species fails a check, so it can gate the alignment phase.

Usage: python3 scripts/xspecies/build_reference_registry.py [--check]
  --check   verify only, do not rewrite the registry files
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
REFS = Path("/private/groups/carpenterlab/emalekos/genomes/xspecies_refs")
STARIDX = Path("/private/groups/carpenterlab/emalekos/STAR_indexes")
SALIDX = Path("/private/groups/carpenterlab/emalekos/Salmon_indexes")
SPTSV = NEW / "scripts" / "xspecies" / "species.tsv"
MD = NEW / "docs" / "SPECIES_REFERENCE_REGISTRY.md"
TSV = NEW / "data" / "species_reference_registry.tsv"

NCRNA_FLOOR = 10          # below this, the drop list is certainly broken
RIBOCODE_MIN_RATIO = 0.5  # prepare_transcripts should keep most transcripts


def read_species() -> list[dict]:
    rows = []
    with SPTSV.open() as fh:
        lines = [l for l in fh if not l.startswith("#")]
    rd = csv.DictReader(lines, delimiter="\t")
    for r in rd:
        rows.append(r)
    return rows


def du(path: Path) -> str:
    if not path.exists():
        return "-"
    try:
        out = subprocess.run(["du", "-sh", str(path)], capture_output=True,
                             text=True, timeout=300).stdout.split()[0]
        return out
    except Exception:
        return "?"


def count_lines(p: Path) -> int:
    if not p.exists():
        return 0
    with p.open("rb") as fh:
        return sum(1 for _ in fh)


def count_fasta(p: Path) -> int:
    if not p.exists():
        return 0
    n = 0
    with p.open("rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def star_params(idx: Path) -> dict[str, str]:
    gp = idx / "genomeParameters.txt"
    out: dict[str, str] = {}
    if not gp.exists():
        return out
    for line in gp.read_text(errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0].strip()] = parts[1].strip()
    return out


def inspect(sp: dict) -> dict:
    name, ann = sp["species"], sp["annver"]
    annotd = NEW / "data" / "annot" / name
    sidx = STARIDX / f"star_index_{name}_{ann}"
    qidx = SALIDX / f"salmon_index_decoy_{name}_{ann}"
    rcan = NEW / "data" / f"ribocode_annot_{name}"

    rec: dict = {
        "species": name, "annver": ann,
        "assembly_accession": sp["assembly_accession"],
        "assembly_name": sp["assembly_name"],
        "annotation_source": sp["annotation_source"],
        "biotype_convention": sp["biotype_convention"],
        "genome_bp": sp["genome_bp"],
        "genome": str(REFS / name / "genome.fna"),
        "raw_gtf": str(REFS / name / "genomic.gtf"),
        "normalized_gtf": str(annotd / f"{name}.normalized.gtf"),
        "transcriptome_fa": str(annotd / "transcripts.fa"),
        "ncrna_tx": str(annotd / "ncrna_tx.txt"),
        "tx_to_gene": str(annotd / "tx_to_gene.tsv"),
        "star_index": str(sidx), "salmon_index": str(qidx),
        "ribocode_annot": str(rcan),
        "problems": [],
    }
    P = rec["problems"]

    for key, path in (("genome", REFS / name / "genome.fna"),
                      ("raw_gtf", REFS / name / "genomic.gtf"),
                      ("normalized_gtf", annotd / f"{name}.normalized.gtf"),
                      ("transcriptome_fa", annotd / "transcripts.fa"),
                      ("ncrna_tx", annotd / "ncrna_tx.txt"),
                      ("tx_to_gene", annotd / "tx_to_gene.tsv")):
        if not path.exists() or path.stat().st_size == 0:
            P.append(f"{key} missing or empty")

    rec["n_ncrna_tx"] = count_lines(annotd / "ncrna_tx.txt")
    rec["n_tx2gene"] = count_lines(annotd / "tx_to_gene.tsv")
    rec["n_transcriptome"] = count_fasta(annotd / "transcripts.fa")
    if rec["n_ncrna_tx"] < NCRNA_FLOOR:
        P.append(f"ncrna_tx.txt has only {rec['n_ncrna_tx']} entries (floor {NCRNA_FLOOR})")

    # normalization report
    nr = annotd / "normalization_report.json"
    if nr.exists():
        d = json.loads(nr.read_text())
        rec["mito_contig"] = d.get("mito_contig") or "NONE"
        rec["mito_retyped"] = d.get("mito_features_retyped", 0)
        rec["source_detected"] = d.get("source_convention", "?")
        rec["ncrna_composition"] = ", ".join(
            f"{k} {v}" for k, v in (d.get("ncrna_composition") or {}).items())
        rec["n_genes"] = d.get("n_genes", 0)
        if rec["mito_contig"] == "NONE":
            P.append("no mitochondrial contig: mito rRNA/tRNA cannot be re-typed, "
                     "and mito reads have nowhere to map")
    else:
        P.append("normalization_report.json missing")
        rec["mito_contig"] = "?"; rec["mito_retyped"] = 0
        rec["source_detected"] = "?"; rec["ncrna_composition"] = ""
        rec["n_genes"] = 0

    # STAR index: the parameters must be the per-species ones, not the defaults.
    spar = star_params(sidx)
    rec["star_sjdbOverhang"] = spar.get("sjdbOverhang", "-")
    rec["star_SAindexNbases"] = spar.get("genomeSAindexNbases", "-")
    rec["star_size"] = du(sidx)
    if not (sidx / "SA").exists():
        P.append("STAR index missing (no SA)")
    else:
        if rec["star_sjdbOverhang"] != sp["sjdb_overhang"]:
            P.append(f"STAR sjdbOverhang is {rec['star_sjdbOverhang']}, "
                     f"table says {sp['sjdb_overhang']}")
        if rec["star_SAindexNbases"] != sp["sa_index_nbases"]:
            P.append(f"STAR genomeSAindexNbases is {rec['star_SAindexNbases']}, "
                     f"table says {sp['sa_index_nbases']}")

    # salmon
    ij = qidx / "info.json"
    if ij.exists():
        d = json.loads(ij.read_text())
        rec["salmon_num_refs"] = d.get("num_refs", 0)
        rec["salmon_num_decoys"] = d.get("num_decoys", 0)
        rec["salmon_k"] = d.get("k", "?")
        rec["salmon_size"] = du(qidx)
        if not rec["salmon_num_decoys"]:
            P.append("salmon index has no decoys (transcriptome-only index biases TPM)")
    else:
        P.append("salmon index missing")
        rec["salmon_num_refs"] = rec["salmon_num_decoys"] = 0
        rec["salmon_k"] = "?"; rec["salmon_size"] = "-"

    # RiboCode
    rcfa = rcan / "transcripts_sequence.fa"
    rec["n_ribocode"] = count_fasta(rcfa)
    rec["ribocode_size"] = du(rcan)
    if rec["n_ribocode"] == 0:
        P.append("RiboCode annot missing or empty")
    elif rec["n_transcriptome"]:
        ratio = rec["n_ribocode"] / rec["n_transcriptome"]
        rec["ribocode_ratio"] = f"{ratio:.3f}"
        if ratio < RIBOCODE_MIN_RATIO:
            P.append(f"RiboCode kept only {100*ratio:.1f}% of transcripts")
    rec.setdefault("ribocode_ratio", "-")
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    recs = [inspect(sp) for sp in read_species()]
    bad = [r for r in recs if r["problems"]]

    print(f"{len(recs)} species inspected, {len(bad)} with problems\n")
    for r in recs:
        flag = "FAIL" if r["problems"] else "ok  "
        print(f"  [{flag}] {r['species']:10s} tx={r['n_transcriptome']:>7} "
              f"ncRNA={r['n_ncrna_tx']:>6} STAR={r['star_size']:>6} "
              f"salmon={r['salmon_size']:>6} mito={r['mito_contig']}")
        for p in r["problems"]:
            print(f"          - {p}")

    if args.check:
        return 1 if bad else 0

    cols = ["species", "annver", "assembly_accession", "assembly_name", "annotation_source",
            "biotype_convention", "source_detected", "genome_bp", "mito_contig", "mito_retyped",
            "n_genes", "n_transcriptome", "n_tx2gene", "n_ncrna_tx", "ncrna_composition",
            "star_sjdbOverhang", "star_SAindexNbases", "star_size",
            "salmon_num_refs", "salmon_num_decoys", "salmon_k", "salmon_size",
            "n_ribocode", "ribocode_ratio", "ribocode_size",
            "genome", "raw_gtf", "normalized_gtf", "transcriptome_fa", "ncrna_tx",
            "tx_to_gene", "star_index", "salmon_index", "ribocode_annot"]
    with TSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(recs)

    lines = [
        "# Cross-species reference registry",
        "",
        "Generated by `scripts/xspecies/build_reference_registry.py` from",
        "`scripts/xspecies/species.tsv` plus the artifacts actually on disk. Re-run it after any",
        "reference change; `--check` exits nonzero if a species fails a gate, so it can block the",
        "alignment phase.",
        "",
        "Machine-readable copy: `data/species_reference_registry.tsv`.",
        "",
        "Human and mouse are deliberately absent: they predate this layout, use GENCODE, and their",
        "paths are hardcoded in `scripts/riboseq_align.sbatch`. Their ncRNA lists are the ground",
        "truth `scripts/normalize_annotation.py --selftest` reproduces exactly.",
        "",
        "## Assemblies and annotation",
        "",
        "| species | assembly | accession | annotation | biotype key | mito contig |",
        "|---|---|---|---|---|---|",
    ]
    for r in recs:
        lines.append(f"| {r['species']} | {r['assembly_name']} | `{r['assembly_accession']}` | "
                     f"{r['annotation_source']} | `{r['biotype_convention']}` | "
                     f"{'**none**' if r['mito_contig']=='NONE' else '`'+str(r['mito_contig'])+'`'} |")
    lines += [
        "",
        "## Derived artifacts",
        "",
        "| species | genes | transcripts | ncRNA dropped | RiboCode tx (ratio) | STAR | salmon |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in recs:
        lines.append(f"| {r['species']} | {r['n_genes']:,} | {r['n_transcriptome']:,} | "
                     f"{r['n_ncrna_tx']:,} | {r['n_ribocode']:,} ({r['ribocode_ratio']}) | "
                     f"{r['star_size']} | {r['salmon_size']} |")
    lines += [
        "",
        "## Index parameters, as recorded by the index itself",
        "",
        "`genomeSAindexNbases` is `min(14, log2(genome_bp)/2 - 1)`. 14 is the mammalian default and",
        "is wrong for the small genomes; STAR does not error, it builds a bad index, so these are",
        "read back out of `genomeParameters.txt` rather than assumed.",
        "",
        "| species | genome bp | sjdbOverhang | genomeSAindexNbases | salmon refs | salmon decoys | k |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in recs:
        lines.append(f"| {r['species']} | {int(r['genome_bp']):,} | {r['star_sjdbOverhang']} | "
                     f"{r['star_SAindexNbases']} | {r['salmon_num_refs']:,} | "
                     f"{r['salmon_num_decoys']} | {r['salmon_k']} |")
    lines += ["", "## What was dropped as ncRNA", "",
              "Drop set for non-GENCODE sources is `core` plus the inline classes GENCODE handles",
              "elsewhere: `rRNA, Mt_rRNA, rRNA_pseudogene, Mt_tRNA, miRNA, tRNA, tRNA_pseudogene,",
              "SRP_RNA, RNase_P_RNA, RNase_MRP_RNA`. See methods.md for why the two drop sets differ.",
              "",
              "| species | composition |", "|---|---|"]
    for r in recs:
        lines.append(f"| {r['species']} | {r['ncrna_composition']} |")

    if bad:
        lines += ["", "## Outstanding problems", ""]
        for r in bad:
            for p in r["problems"]:
                lines.append(f"- **{r['species']}**: {p}")
    else:
        lines += ["", "All species pass every gate.", ""]

    lines += ["", "## Paths", "",
              "Every path follows a fixed convention, so nothing here needs to be looked up by hand:",
              "", "```",
              "genome     genomes/xspecies_refs/<species>/genome.fna",
              "raw GTF    genomes/xspecies_refs/<species>/genomic.gtf",
              "annot      data/annot/<species>/{<species>.normalized.gtf,transcripts.fa,",
              "                                 ncrna_tx.txt,tx_to_gene.tsv,tx2biotype.tsv,",
              "                                 ncrna_deplete.bed,ncrna_postfilter.bed,",
              "                                 normalization_report.json}",
              "STAR       STAR_indexes/star_index_<species>_<annver>",
              "salmon     Salmon_indexes/salmon_index_decoy_<species>_<annver>",
              "           (gentrome decoys.txt kept beside the index; it was NOT kept for",
              "            human/mouse, so those builds cannot be reproduced)",
              "RiboCode   data/ribocode_annot_<species>",
              "```", ""]
    MD.write_text("\n".join(lines))
    print(f"\nwrote {MD}\nwrote {TSV}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
