#!/usr/bin/env python3
"""Supplement a reference assembly with a mitochondrial genome from a separate
RefSeq nucleotide record, emitting GTF in the same shape RefSeq uses for the
assemblies that already include one.

WHY THIS EXISTS
---------------
`GCF_049350105.2` (T2T-MMU8v2.0, rhesus macaque) ships 22 contigs -- chr1-20, X,
Y -- and NO mitochondrial sequence. Every other species in the cross-species
expansion carries the standard 22 Mt_tRNA + 2 Mt_rRNA; macaque has zero. The
consequence is not cosmetic: with no chrM in the index, mitochondrial footprints
(a real and often large fraction of a Ribo-seq library) have nowhere correct to
map, and the likely destination is a NUMT, which a T2T assembly resolves in
unusual completeness. That is silent signal corruption on one of the four species
in the iPSC-CM comparison.

THE INCOMPLETE STOP CODON, WHICH DRIVES THE DESIGN
--------------------------------------------------
Metazoan mitochondrial mRNAs routinely end in T or TA, with the stop codon
completed by polyadenylation. NCBI annotates this explicitly: 6 of the 13 macaque
mito CDS carry `transl_except=(pos:N,aa:TERM)` and a Note reading "TAA stop codon
is completed by the addition of 3' A residues to the mRNA", and their CDS lengths
are 1 mod 3.

So a stop codon is emitted ONLY for the CDS that actually have one in the genome
(length divisible by 3, stop taken as the final 3 bases and excluded from the CDS
feature, which is the RefSeq GTF convention). For the incomplete ones the CDS is
emitted whole and no stop_codon is written, because writing one would place it
outside the transcript's exon -- exactly the condition that makes RiboCode's
prepare_transcripts die with "Can't transform the genomic interval", as it did on
the FlyBase mt:ND2/CoII/ND4/ND5 transcripts.

Usage:
    add_mito_contig.py --species macaque --accession NC_005943.1 [--out-dir <dir>]
                       [--dry-run]

Idempotent: refuses to append a contig the FASTA already contains.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REFS = Path("/private/groups/carpenterlab/emalekos/genomes/xspecies_refs")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
SVIEWER = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"


def fetch(url: str, params: dict) -> str:
    full = f"{url}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(full, timeout=180) as fh:
        return fh.read().decode()


def gff_attrs(field: str) -> dict[str, str]:
    out = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = urllib.parse.unquote(v)
    return out


def gtf_attrs(pairs: list[tuple[str, str]]) -> str:
    return " ".join(f'{k} "{v}";' for k, v in pairs)


def convert(gff: str, contig: str) -> tuple[list[str], dict]:
    """GFF3 -> RefSeq-shaped GTF. Returns (lines, stats)."""
    genes: dict[str, dict] = {}
    rnas: dict[str, dict] = {}
    cds: dict[str, list] = {}
    exons: dict[str, list] = {}

    for line in gff.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9 or f[0] != contig:
            continue
        a = gff_attrs(f[8])
        feat, start, end, strand = f[2], int(f[3]), int(f[4]), f[6]
        if feat == "gene":
            genes[a["ID"]] = {"start": start, "end": end, "strand": strand, "attrs": a}
        elif feat in ("mRNA", "tRNA", "rRNA", "ncRNA", "transcript"):
            rnas[a["ID"]] = {"start": start, "end": end, "strand": strand,
                             "attrs": a, "parent": a.get("Parent", ""), "kind": feat}
        elif feat == "CDS":
            cds.setdefault(a.get("Parent", ""), []).append((start, end, f[7], a))
        elif feat == "exon":
            exons.setdefault(a.get("Parent", ""), []).append((start, end))

    lines: list[str] = []
    n_tx = 0
    stats = {"gene": 0, "transcript": 0, "exon": 0, "CDS": 0,
             "start_codon": 0, "stop_codon": 0, "incomplete_stop": 0}

    def emit(feat, s, e, strand, frame, attrs_pairs):
        lines.append("\t".join([contig, "RefSeq", feat, str(s), str(e), ".",
                                strand, frame, gtf_attrs(attrs_pairs) + " "]))
        stats[feat] = stats.get(feat, 0) + 1

    for gid, g in sorted(genes.items(), key=lambda kv: kv[1]["start"]):
        ga = g["attrs"]
        locus = ga.get("locus_tag", ga.get("Name", gid))
        biotype = ga.get("gene_biotype", "other")
        gene_name = ga.get("gene", "")
        base = [("gene_id", locus)]
        gene_pairs = base + [("transcript_id", "")]
        if "Dbxref" in ga:
            gene_pairs += [("db_xref", ga["Dbxref"].split(",")[0])]
        gene_pairs += [("gbkey", "Gene")]
        if gene_name:
            gene_pairs += [("gene", gene_name)]
        gene_pairs += [("gene_biotype", biotype), ("locus_tag", locus)]
        emit("gene", g["start"], g["end"], g["strand"], ".", gene_pairs)

        # The RNA child, if the GFF3 has one; protein-coding mito genes often have
        # only a CDS, so fall back to the gene span.
        child = next((r for r in rnas.values() if r["parent"] == gid), None)
        gene_cds = cds.get(f"cds-{ga.get('Name','')}", [])
        if not gene_cds:
            gene_cds = [c for parent, cl in cds.items() if parent == gid for c in cl]
        if not gene_cds and child is None and biotype == "protein_coding":
            gene_cds = [c for parent, cl in cds.items()
                        for c in cl if c[3].get("locus_tag") == locus]

        n_tx += 1
        tx_id = f"mito_transcript_{n_tx}"
        tx_kind = {"protein_coding": "mRNA"}.get(biotype, biotype)
        tstart = child["start"] if child else g["start"]
        tend = child["end"] if child else g["end"]
        product = (child["attrs"].get("product") if child else None) or \
                  (gene_cds[0][3].get("product") if gene_cds else None)

        tx_pairs = base + [("transcript_id", tx_id), ("gbkey", tx_kind)]
        if gene_name:
            tx_pairs += [("gene", gene_name)]
        tx_pairs += [("locus_tag", locus)]
        if product:
            tx_pairs += [("product", product)]
        tx_pairs += [("transcript_biotype", tx_kind)]
        emit("transcript", tstart, tend, g["strand"], ".", tx_pairs)

        child_id = child["attrs"]["ID"] if child else ""
        ex = exons.get(child_id, []) or exons.get(f"rna-{locus}", []) or [(tstart, tend)]
        for i, (es, ee) in enumerate(sorted(ex), 1):
            ex_pairs = base + [("transcript_id", tx_id)]
            if gene_name:
                ex_pairs += [("gene", gene_name)]
            ex_pairs += [("locus_tag", locus)]
            if product:
                ex_pairs += [("product", product)]
            ex_pairs += [("transcript_biotype", tx_kind), ("exon_number", str(i))]
            emit("exon", es, ee, g["strand"], ".", ex_pairs)

        if not gene_cds:
            continue
        cs, ce, frame, ca = sorted(gene_cds)[0]
        length = ce - cs + 1
        complete = (length % 3 == 0) and ("transl_except" not in ca)
        c_pairs = base + [("transcript_id", tx_id), ("gbkey", "CDS")]
        if gene_name:
            c_pairs += [("gene", gene_name)]
        c_pairs += [("locus_tag", locus)]
        if ca.get("product"):
            c_pairs += [("product", ca["product"])]
        if ca.get("protein_id"):
            c_pairs += [("protein_id", ca["protein_id"])]
        c_pairs += [("transl_table", ca.get("transl_table", "2")), ("exon_number", "1")]

        if complete:
            # RefSeq GTF convention: the CDS feature EXCLUDES the stop codon, which
            # is emitted separately.
            if g["strand"] == "+":
                emit("CDS", cs, ce - 3, g["strand"], frame, c_pairs)
                emit("stop_codon", ce - 2, ce, g["strand"], "0", c_pairs)
                emit("start_codon", cs, cs + 2, g["strand"], "0", c_pairs)
            else:
                emit("CDS", cs + 3, ce, g["strand"], frame, c_pairs)
                emit("stop_codon", cs, cs + 2, g["strand"], "0", c_pairs)
                emit("start_codon", ce - 2, ce, g["strand"], "0", c_pairs)
        else:
            # Incomplete stop: the codon is finished by polyadenylation and does not
            # exist in the genome. Emit the CDS whole and NO stop_codon; writing one
            # would put it outside the exon and kill RiboCode.
            stats["incomplete_stop"] += 1
            emit("CDS", cs, ce, g["strand"], frame, c_pairs)
            if g["strand"] == "+":
                emit("start_codon", cs, cs + 2, g["strand"], "0", c_pairs)
            else:
                emit("start_codon", ce - 2, ce, g["strand"], "0", c_pairs)

    return lines, stats


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--ref-dir", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    d = args.ref_dir or (REFS / args.species)
    genome, gtf = d / "genome.fna", d / "genomic.gtf"
    for p in (genome, gtf):
        if not p.exists():
            sys.exit(f"missing {p}")

    acc = args.accession
    print(f"species={args.species} accession={acc}\n  ref dir {d}")

    # Idempotence: never append twice.
    present = subprocess.run(["grep", "-c", f"^>{acc}", str(genome)],
                             capture_output=True, text=True).stdout.strip()
    if present and present != "0":
        print(f"  {acc} is ALREADY in {genome.name}; nothing to do")
        return 0

    print("  fetching FASTA ...", flush=True)
    fa = fetch(EUTILS, {"db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text"})
    hdr = fa.splitlines()[0]
    seq_len = sum(len(l.strip()) for l in fa.splitlines()[1:] if l.strip())
    print(f"    {hdr[:80]}\n    {seq_len} bp")
    if not hdr.startswith(f">{acc}"):
        sys.exit(f"unexpected FASTA header: {hdr[:120]}")
    if "mitochondrion" not in hdr.lower():
        sys.exit(f"header does not look mitochondrial: {hdr[:120]}")

    print("  fetching GFF3 ...", flush=True)
    gff = fetch(SVIEWER, {"id": acc, "db": "nuccore", "report": "gff3", "retmode": "text"})
    lines, stats = convert(gff, acc)
    print("    emitted: " + ", ".join(f"{k}={v}" for k, v in stats.items() if v))

    # A mammalian mitochondrion is 13 CDS + 22 tRNA + 2 rRNA = 37 genes. Anything
    # else means the record or the conversion is not what was expected.
    if stats["gene"] != 37:
        print(f"    WARNING: {stats['gene']} genes, expected 37 for a standard "
              f"metazoan mitochondrion", file=sys.stderr)
    if stats["transcript"] != stats["gene"]:
        sys.exit(f"transcript count {stats['transcript']} != gene count {stats['gene']}")

    # Every codon feature must be inside its transcript's exons, or the normalizer
    # will drop the transcript and RiboCode would have died anyway.
    ex, cod = {}, []
    for l in lines:
        f = l.split("\t")
        t = re.search(r'transcript_id "([^"]*)"', f[8]).group(1)
        if f[2] == "exon":
            ex.setdefault(t, []).append((int(f[3]), int(f[4])))
        elif f[2] in ("start_codon", "stop_codon"):
            cod.append((t, f[2], int(f[3]), int(f[4])))
    bad = []
    for t, feat, s, e in cod:
        need = e - s + 1
        for a, b in ex.get(t, []):
            if b < s or a > e:
                continue
            need -= min(b, e) - max(a, s) + 1
        if need > 0:
            bad.append(f"{t} {feat} {s}-{e}")
    if bad:
        sys.exit(f"codon(s) outside exons, refusing to write: {bad[:5]}")
    print(f"    all {len(cod)} codon features contained within exons")

    if args.dry_run:
        print("  (dry run, nothing written)")
        print("\n".join(lines[:6]))
        return 0

    before = {"genome_md5": md5(genome), "gtf_md5": md5(gtf),
              "genome_bytes": genome.stat().st_size, "gtf_bytes": gtf.stat().st_size}
    with genome.open("a") as fh:
        if not fa.endswith("\n"):
            fh.write("\n")
        fh.write(fa if fa.endswith("\n") else fa + "\n")
    with gtf.open("a") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  appended {seq_len} bp and {len(lines)} GTF lines")

    # Appending to a FASTA INVALIDATES its .fai: the index records byte offsets and
    # a contig list, so a stale one simply does not contain the new sequence. gffread
    # then fails with "couldn't find fasta record for '<acc>'" even though the record
    # is plainly in the file, and any tool that merely looks the contig up would
    # silently miss it instead. Delete the index and let the next tool rebuild it.
    for idx in (genome.with_suffix(genome.suffix + ".fai"),
                genome.with_suffix(genome.suffix + ".gzi")):
        if idx.exists():
            idx.unlink()
            print(f"  removed now-stale {idx.name} (a FASTA append invalidates it)")

    prov = d / "PROVENANCE.json"
    p = json.loads(prov.read_text()) if prov.exists() else {}
    p.setdefault("supplements", []).append({
        "what": "mitochondrial genome appended from a separate RefSeq record",
        "accession": acc, "definition": hdr.lstrip(">").strip(), "length_bp": seq_len,
        "why": ("the primary assembly ships no mitochondrial sequence, so mito reads "
                "would have nowhere to map and could land on NUMTs instead"),
        "gtf_features": stats,
        "note": ("stop_codon omitted for the CDS whose stop is completed by "
                 "polyadenylation (transl_except / length not divisible by 3); "
                 "writing one would place it outside the exon"),
        "added_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "before": before,
        "after": {"genome_md5": md5(genome), "gtf_md5": md5(gtf),
                  "genome_bytes": genome.stat().st_size, "gtf_bytes": gtf.stat().st_size},
        "reproduce": f"scripts/xspecies/add_mito_contig.py --species {args.species} --accession {acc}",
    })
    prov.write_text(json.dumps(p, indent=2) + "\n")
    print(f"  recorded in {prov}")
    print("\n  NEXT: re-run normalize_annotation.py and REBUILD the STAR/salmon/RiboCode\n"
          "        artifacts for this species; the existing ones predate the mitochondrion.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
