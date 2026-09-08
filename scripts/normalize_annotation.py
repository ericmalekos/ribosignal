#!/usr/bin/env python3
"""Normalize any GTF onto the GENCODE biotype vocabulary and emit the derived
artifacts the Ribo-seq pipeline consumes.

WHY THIS EXISTS
---------------
The pipeline's ncRNA drop step -- the standing rule that rRNA/tRNA/miRNA/Mt loci
never reach an ORF caller -- was written against GENCODE and matches on the
attribute key `gene_type`. Measured on the references fetched for the
cross-species expansion:

    gene_type occurrences in genomes/*/genomic.gtf (RefSeq) ......... 0
    gene_type occurrences in Ensembl 116 GTFs ....................... 0

Both RefSeq and Ensembl use `gene_biotype`. GENCODE is the outlier. So the
existing matcher silently produces an EMPTY drop list on every non-GENCODE
reference, and the postfilter BED it writes has no non-empty guard at all. The
job succeeds, the BED is zero bytes, and the contamination is never removed.

Ensembl was evaluated as an escape route and rejected: it shares RefSeq's
attribute key, its vocabulary is an older GENCODE variant that still needs
mapping (`lincRNA`/`antisense` for `lncRNA`), and for the three primates its
current assemblies are gorGor4 (2014), Pan_tro_3.0 (2016) and Mmul_10 (2019)
against the RefSeq T2T builds already on disk.

FOUR CONVENTIONS ARE SUPPORTED
------------------------------
  gencode   `gene_type` / `transcript_type`. The target vocabulary.
  refseq    `gene_biotype` / `transcript_biotype`. No Mt_rRNA/Mt_tRNA at all --
            mitochondrial rRNA and tRNA are typed plain rRNA/tRNA and are only
            distinguishable by contig, so the mito contig is detected and those
            features are re-typed on the way through. Transcript-level terms
            differ from gene-level ones (`mRNA` not `protein_coding`), the
            generic `transcript` term defers to the gene, and the vocabulary is
            not even self-consistent (`lnc_RNA` in gorilla and zebrafish,
            `lncRNA` in chimp, macaque and worm).
  ensembl   Same keys as RefSeq, older GENCODE vocabulary.
  flybase   No biotype attribute anywhere. The biotype is the feature-type
            column, and `mRNA` appears where a GTF normally says `transcript`.

TWO DROP SETS, AND WHY THEY DIFFER ON PURPOSE
---------------------------------------------
`core` is exactly the set that reproduces the hand-built human and mouse ground
truth on disk (verified: 2,447/2,447 and 2,579/2,579, zero missing, zero extra):

    rRNA, Mt_rRNA, rRNA_pseudogene, Mt_tRNA, miRNA

`extended` adds tRNA, tRNA_pseudogene, SRP_RNA, RNase_P_RNA, RNase_MRP_RNA and
is the default for every non-GENCODE source. That is not an inconsistency, it is
the point: GENCODE's main GTF contains no cytoplasmic tRNA (they ship separately
in gencode.*.tRNAs.gtf.gz and are handled by the bowtie2 contaminant index) and
types 7SL / RNase-P / RNase-MRP as plain misc_RNA. RefSeq annotates all of them
inline -- zebrafish alone has 8,861 tRNA transcripts -- so if they are not
dropped here they are not dropped anywhere. 7SL in particular is the RNA that
crushed the GSE243134 liver universe to 3,321 transcripts when it survived a
biotype filter as lncRNA.

Usage:
    normalize_annotation.py --gtf <in.gtf> --species <label> [--genome <fna>]
        [--source auto|gencode|refseq|ensembl|flybase]
        [--drop-set auto|core|extended] [--out-dir <dir>] [--write-gtf]

Self-test against the on-disk ground truth (run this before trusting it on a new
species):
    normalize_annotation.py --selftest
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths  # noqa: E402

DATA = paths.data_dir()
# GTFs are inputs, not repo content: $RIBO_ANNOT_DIR, else <data>/annotations.
ANNOT = paths._env_path("RIBO_ANNOT_DIR", DATA / "annotations")

# ---------------------------------------------------------------------------
# Drop sets, expressed in NORMALIZED (GENCODE) terms.
# ---------------------------------------------------------------------------
CORE_DROP = frozenset({"rRNA", "Mt_rRNA", "rRNA_pseudogene", "Mt_tRNA", "miRNA"})
EXTENDED_ONLY = frozenset({"tRNA", "tRNA_pseudogene",
                           "SRP_RNA", "RNase_P_RNA", "RNase_MRP_RNA"})
EXTENDED_DROP = CORE_DROP | EXTENDED_ONLY

# Second-tier contaminants: not removed from the transcriptome BAM, but written
# as a genomic BED for the post-alignment `samtools view -L` filter. Matches the
# POSTFILT set in scripts/build_contaminant_index.sbatch.
POSTFILTER = frozenset({"snoRNA", "snRNA", "scaRNA", "misc_RNA",
                        "vault_RNA", "sRNA", "scRNA", "ribozyme"})

# ---------------------------------------------------------------------------
# Vocabulary maps: source term -> GENCODE term.
#
# Built by enumerating the ACTUAL vocabularies of every GTF in scope, not from
# memory. A term absent from the relevant map is a hard error, never a silent
# pass-through: that is what makes this safe to point at an eighth species.
#
# PASSTHROUGH means "no GENCODE equivalent exists, keep the source term". It is
# a deliberate decision recorded per term, not a fallback. Nothing marked
# PASSTHROUGH is in any drop set.
# ---------------------------------------------------------------------------
PASSTHROUGH = "__passthrough__"
# Defer to the gene's biotype: RefSeq uses a generic transcript-level term for
# transcripts whose class is only stated on the gene.
FROM_GENE = "__from_gene__"

REFSEQ_MAP = {
    # coding
    "mRNA": "protein_coding", "protein_coding": "protein_coding",
    # IG/TR segments. RefSeq does not say whether a segment is immunoglobulin or
    # T-cell receptor, and GENCODE's vocabulary (IG_V_gene vs TR_V_gene, IG_J_gene vs
    # TR_J_gene, ...) requires that distinction, so mapping either way would be wrong
    # half the time. PASSTHROUGH keeps the source term; none of these is in a drop set.
    # J_ and D_ segments appear only in the human and mouse RefSeq annotations, which
    # are far richer at these loci than the other seven species.
    "C_gene_segment": PASSTHROUGH, "V_gene_segment": PASSTHROUGH,
    "J_gene_segment": PASSTHROUGH, "D_gene_segment": PASSTHROUGH,
    "C_region": PASSTHROUGH, "V_segment": PASSTHROUGH,
    "vault_RNA": "vault_RNA",
    "Y_RNA": "misc_RNA",       # Ro60 ribonucleoprotein component; misc_RNA in GENCODE
    # long non-coding. RefSeq spells it both ways.
    "lncRNA": "lncRNA", "lnc_RNA": "lncRNA",
    "antisense_RNA": "lncRNA",          # merged into lncRNA in GENCODE v29+
    "telomerase_RNA": "lncRNA",         # TERC is lncRNA in GENCODE
    # small RNA
    "miRNA": "miRNA",
    "primary_transcript": "miRNA",      # pre-miRNA; GENCODE types these miRNA
    "snRNA": "snRNA", "snoRNA": "snoRNA", "scaRNA": "scaRNA",
    "scRNA": "scRNA", "misc_RNA": "misc_RNA",
    "ncRNA": "misc_RNA",                # RefSeq/WormBase catch-all
    "piRNA": PASSTHROUGH,               # 15,363 in C. elegans; no GENCODE class
    # contaminant classes
    "rRNA": "rRNA", "tRNA": "tRNA",
    "SRP_RNA": "SRP_RNA",               # 7SL
    "RNase_P_RNA": "RNase_P_RNA", "RNase_MRP_RNA": "RNase_MRP_RNA",
    "pseudogenic_tRNA": "tRNA_pseudogene", "tRNA_pseudogene": "tRNA_pseudogene",
    "pseudogenic_rRNA": "rRNA_pseudogene", "rRNA_pseudogene": "rRNA_pseudogene",
    # pseudogene / unknown
    "pseudogene": "pseudogene",
    "transcribed_pseudogene": "transcribed_unprocessed_pseudogene",
    "transposable_element": PASSTHROUGH,
    "other": PASSTHROUGH,               # RefSeq's own unknown bucket
    # generic transcript-level term: the gene knows what it is
    "transcript": FROM_GENE,
}

ENSEMBL_MAP = dict(REFSEQ_MAP)
ENSEMBL_MAP.update({
    # miRNA precursors, as Ensembl Metazoa spells them. GENCODE types pre-miRNA as
    # miRNA, and RefSeq's equivalent term `primary_transcript` is mapped the same way.
    "pre_miRNA": "miRNA", "miRNA_primary_transcript": "miRNA",
    "lincRNA": "lncRNA", "antisense": "lncRNA",
    "processed_transcript": "processed_transcript",
    "sense_intronic": "lncRNA", "sense_overlapping": "lncRNA",
    "3prime_overlapping_ncRNA": "lncRNA", "bidirectional_promoter_lncRNA": "lncRNA",
    "macro_lncRNA": "lncRNA", "non_coding": "lncRNA",
    "Mt_rRNA": "Mt_rRNA", "Mt_tRNA": "Mt_tRNA",
    "TEC": "TEC", "ribozyme": "ribozyme", "sRNA": "sRNA",
    "vault_RNA": "vault_RNA", "vaultRNA": "vault_RNA",
    "unprocessed_pseudogene": "unprocessed_pseudogene",
    "processed_pseudogene": "processed_pseudogene",
    "unitary_pseudogene": "unitary_pseudogene",
    "polymorphic_pseudogene": "protein_coding",
    "transcribed_unprocessed_pseudogene": "transcribed_unprocessed_pseudogene",
    "transcribed_processed_pseudogene": "transcribed_processed_pseudogene",
    "transcribed_unitary_pseudogene": "transcribed_unitary_pseudogene",
    "translated_processed_pseudogene": "translated_processed_pseudogene",
    "translated_unprocessed_pseudogene": "translated_unprocessed_pseudogene",
    "nonsense_mediated_decay": "nonsense_mediated_decay",
    "non_stop_decay": "non_stop_decay", "retained_intron": "retained_intron",
    "protein_coding_CDS_not_defined": "protein_coding_CDS_not_defined",
    "protein_coding_LoF": "protein_coding_LoF", "artifact": "artifact",
})
for _p in ("IG_V_gene", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene",
           "IG_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene",
           "IG_D_pseudogene", "IG_J_pseudogene", "TR_V_gene", "TR_C_gene",
           "TR_D_gene", "TR_J_gene", "TR_V_pseudogene", "TR_J_pseudogene"):
    ENSEMBL_MAP[_p] = _p

# FlyBase keeps the biotype in the feature-type column and has no Mt_* classes.
FLYBASE_MAP = {
    "mRNA": "protein_coding", "ncRNA": "misc_RNA", "miRNA": "miRNA",
    "pre_miRNA": "miRNA", "pseudogene": "pseudogene", "tRNA": "tRNA",
    "snoRNA": "snoRNA", "snRNA": "snRNA", "rRNA": "rRNA",
}

# GENCODE is already the target; the map exists so `--source gencode` runs the
# same code path and the self-test exercises it.
GENCODE_MAP: dict[str, str] = {}

MAPS = {"refseq": REFSEQ_MAP, "ensembl": ENSEMBL_MAP,
        "flybase": FLYBASE_MAP, "gencode": GENCODE_MAP}

# Transcript-level feature names to read, per convention. FlyBase is the reason
# this is not just {"transcript"}.
TX_FEATURES = {
    "gencode": {"transcript"}, "refseq": {"transcript"}, "ensembl": {"transcript"},
    "flybase": set(FLYBASE_MAP),
}

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
MITO_HINT = re.compile(r"mitochondri|\bMT\b|chrM", re.I)


def opener(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attrs(field: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(field))


def detect_source(gtf: Path) -> str:
    """Probe the head of the file for which convention is in use."""
    saw_gene_type = saw_gene_biotype = False
    ens = refseq = 0
    with opener(gtf) as fh:
        for i, line in enumerate(fh):
            if i > 200_000:
                break
            if line.startswith("#"):
                continue
            if "gene_type " in line:
                saw_gene_type = True
            if "gene_biotype " in line:
                saw_gene_biotype = True
            m = re.search(r'transcript_id "([^"]*)"', line)
            if m:
                if m.group(1).startswith("ENS"):
                    ens += 1
                elif re.match(r"^[NX][MR]_", m.group(1)):
                    refseq += 1
    if saw_gene_type:
        return "gencode"
    if saw_gene_biotype:
        return "refseq" if refseq >= ens else "ensembl"
    return "flybase"


def detect_mito(gtf: Path, genome: Path | None) -> str | None:
    """Find the mitochondrial contig. RefSeq has no Mt_rRNA/Mt_tRNA biotypes, so
    contig identity is the ONLY way to tell mito rRNA/tRNA from nuclear."""
    if genome and genome.exists():
        with opener(genome) as fh:
            for line in fh:
                if line.startswith(">") and MITO_HINT.search(line):
                    return line[1:].split()[0]
    seen: set[str] = set()
    with opener(gtf) as fh:
        for line in fh:
            if not line.startswith("#"):
                seen.add(line.split("\t", 1)[0])
    for c in ("chrM", "MT", "chrMT", "M", "mitochondrion_genome", "Mito"):
        if c in seen:
            return c
    return None


# ---------------------------------------------------------------------------
# Feature-column repair.
#
# NCBI's RefSeq GTF for C. elegans (GCF_000002985.6, WormBase WS298) writes the
# ISOFORM'S OWN standard_name into the feature-type column where "CDS" belongs:
#
#   NC_003284.9  RefSeq  F31B12.1k  10838615  10838650  .  -  0  gene_id "CELE_F31B12.1"; ...
#                        ^^^^^^^^^ should be CDS
#
# 204,530 of 204,542 CDS rows are hit, across 30,548 distinct names, leaving 12
# genuine "CDS" rows in the whole file. None of the other six RefSeq annotations
# used here (yeast, zebrafish, gorilla, chimp, macaque, human) shows a single odd
# row, so this is one upstream file's defect and not a RefSeq convention.
#
# It is silent and expensive: exon and start_codon rows are untouched, so the GTF
# still builds a correct STAR index, a correct transcriptome FASTA and a correct
# universe. Only the CDS-consuming step notices, and it does not fail -- RiboCode
# simply finds no annotated CDS and types EVERY called ORF "novel". That is what
# it did: 17,745 worm calls, 100% "novel", 0 "annotated", against 5,344/5,358
# annotated for yeast.
#
# The rewrite condition is deliberately narrow, and all three parts held for
# 204,530/204,530 rows: the feature is not a known GTF feature, the frame column
# is a valid 0/1/2 (which only CDS carries), and the feature string equals the
# row's own standard_name. Anything else odd is a HARD ERROR rather than a guess.
STD_FEATURES = frozenset({
    "gene", "transcript", "exon", "CDS", "start_codon", "stop_codon",
    "five_prime_utr", "three_prime_utr", "UTR", "Selenocysteine",
})
_STDNAME_RE = re.compile(r'standard_name "([^"]*)"')


def repair_feature(f: list[str]) -> tuple[str, bool]:
    """Return (feature_type, was_repaired) for one GTF row."""
    feat = f[2]
    if feat in STD_FEATURES:
        return feat, False
    m = _STDNAME_RE.search(f[8])
    if f[7] in ("0", "1", "2") and m and m.group(1) == feat:
        return "CDS", True
    raise SystemExit(
        f"ABORT: unrecognised feature type {feat!r} at {f[0]}:{f[3]}-{f[4]}\n"
        f"  frame={f[7]!r} standard_name={m.group(1) if m else None!r}\n"
        "  Not a known GTF feature, and it does not match the NCBI C. elegans\n"
        "  standard-name-in-feature-column defect (which requires a valid frame and\n"
        "  feature == standard_name). Add it to STD_FEATURES or extend repair_feature()."
    )


class Normalizer:
    def __init__(self, gtf: Path, species: str, source: str, mito: str | None,
                 drop_set: frozenset[str]):
        self.gtf, self.species, self.source = gtf, species, source
        self.mito, self.drop_set = mito, drop_set
        self.vmap = MAPS[source]
        self.txfeat = TX_FEATURES[source]
        self.gene_bt: dict[str, str] = {}
        self.tx: dict[str, dict] = {}
        self.unmapped: Counter = Counter()
        self.bad_strand: dict[str, str] = {}
        self.dup_tx: Counter = Counter()
        self.bad_codon: dict[str, str] = {}
        self.mapped: Counter = Counter()
        self.retyped_mito = 0
        self.repaired_feat = 0

    def _raw_gene_bt(self, a: dict, feat: str) -> str | None:
        if self.source == "gencode":
            return a.get("gene_type")
        if self.source == "flybase":
            return None          # genes carry no type; derived from transcripts
        return a.get("gene_biotype")

    def _raw_tx_bt(self, a: dict, feat: str) -> str | None:
        if self.source == "gencode":
            return a.get("transcript_type")
        if self.source == "flybase":
            return feat
        return a.get("transcript_biotype")

    def _map(self, term: str, gene_id: str) -> str | None:
        if self.source == "gencode":
            return term
        tgt = self.vmap.get(term)
        if tgt is None:
            self.unmapped[term] += 1
            return None
        if tgt == FROM_GENE:
            g = self.gene_bt.get(gene_id)
            tgt = (self.vmap.get(g, PASSTHROUGH) if g else PASSTHROUGH)
            if tgt in (FROM_GENE, PASSTHROUGH):
                tgt = g or "misc_RNA"
        if tgt == PASSTHROUGH:
            tgt = term
        self.mapped[f"{term} -> {tgt}"] += 1
        return tgt

    def _resolve_gene_type(self, gene_id: str, fallback: str) -> str:
        """Gene biotype in GENCODE terms, with PASSTHROUGH resolved to the source term."""
        raw = self.gene_bt.get(gene_id)
        if raw is None:
            return fallback
        tgt = self.vmap.get(raw)
        if tgt is None or tgt == FROM_GENE:
            return raw
        if tgt == PASSTHROUGH:
            return raw
        return tgt

    def pass1_genes(self) -> None:
        """Collect gene-level biotypes so FROM_GENE transcripts can resolve."""
        if self.source in ("gencode", "flybase"):
            return
        with opener(self.gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                a = attrs(f[8])
                gid, bt = a.get("gene_id"), self._raw_gene_bt(a, f[2])
                if gid and bt:
                    self.gene_bt[gid] = bt

    def pass2_transcripts(self) -> None:
        exon_len: dict[str, int] = defaultdict(int)
        exon_iv: dict[str, list[tuple[int, int]]] = defaultdict(list)
        codon_iv: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        with opener(self.gtf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                feat, fixed = repair_feature(f)
                if fixed:
                    self.repaired_feat += 1
                if feat == "exon":
                    a = attrs(f[8])
                    tid = a.get("transcript_id")
                    if tid:
                        exon_len[tid] += int(f[4]) - int(f[3]) + 1
                        exon_iv[tid].append((int(f[3]), int(f[4])))
                    continue
                if feat in ("start_codon", "stop_codon"):
                    a = attrs(f[8])
                    tid = a.get("transcript_id")
                    if tid:
                        codon_iv[tid].append((feat, int(f[3]), int(f[4])))
                    continue
                if feat not in self.txfeat:
                    continue
                a = attrs(f[8])
                tid, gid = a.get("transcript_id"), a.get("gene_id")
                if not tid:
                    continue
                # Strand must be + or -. FlyBase marks the trans-spliced mod(mdg4)
                # transcripts ".", because they are assembled from BOTH strands; there
                # are 9 such records in r6.67. RiboCode's prepare_transcripts raises
                # ValueError('strand is neither "+" nor "-"') and dies on the first one,
                # so they are dropped here with an explicit count rather than left to
                # kill the annotation build.
                if f[6] not in ("+", "-"):
                    self.bad_strand[tid] = f[6]
                    continue
                raw = self._raw_tx_bt(a, feat)
                if raw is None:
                    self.unmapped["<no biotype attribute>"] += 1
                    continue
                bt = self._map(raw, gid or "")
                if bt is None:
                    continue
                # RefSeq has no Mt_* classes: re-type by contig.
                if self.mito and f[0] == self.mito:
                    if bt == "rRNA":
                        bt, self.retyped_mito = "Mt_rRNA", self.retyped_mito + 1
                    elif bt == "tRNA":
                        bt, self.retyped_mito = "Mt_tRNA", self.retyped_mito + 1
                if tid in self.tx:
                    self.dup_tx[tid] += 1
                    continue
                self.tx[tid] = {
                    "tx_id": tid, "gene_id": gid or "",
                    "gene_name": a.get("gene_name") or a.get("gene") or "",
                    "chrom": f[0], "strand": f[6], "start": int(f[3]), "end": int(f[4]),
                    "transcript_type": bt, "raw_type": raw,
                    # Route the gene biotype through the SAME resolver as the
                    # transcript one. A bare self.vmap.get() returns the internal
                    # PASSTHROUGH sentinel for terms with no GENCODE equivalent, and
                    # that string was leaking into the gene_type column and into the
                    # normalized GTF: 15,365 rows of literal "__passthrough__" in
                    # C. elegans (its piRNA genes) and ~200 in each primate
                    # (V_gene_segment / C_gene_segment).
                    "gene_type": (a.get("gene_type", bt) if self.source == "gencode"
                                  else self._resolve_gene_type(gid or "", bt)),
                }
        for tid, rec in self.tx.items():
            rec["length"] = exon_len.get(tid, rec["end"] - rec["start"] + 1)

        # A start/stop codon must lie entirely within the transcript's own exons.
        # Metazoan mitochondrial mRNAs routinely violate this: the transcript ends in
        # T or TA and the stop codon is completed by polyadenylation, so the annotated
        # codon runs 1-2 bp past the last exon. FlyBase annotates the full codon
        # anyway (mt:ND2, mt:CoII, mt:ND4, mt:ND5 in r6.67). RiboCode's
        # prepare_transcripts cannot project such a codon into transcript coordinates
        # and dies with "Can't transform the genomic interval, please check!" -- one
        # bad transcript kills the whole annotation build, so they are dropped here.
        for tid, codons in codon_iv.items():
            if tid not in self.tx:
                continue
            spans = exon_iv.get(tid, [])
            for feat, cs, ce in codons:
                need = set(range(cs, ce + 1))
                for a, b in spans:
                    if b < cs or a > ce:
                        continue
                    need -= set(range(max(a, cs), min(b, ce) + 1))
                    if not need:
                        break
                if need:
                    self.bad_codon[tid] = f"{feat} {cs}-{ce} not contained in exons"
                    break
        for tid in self.bad_codon:
            self.tx.pop(tid, None)

    def run(self) -> None:
        self.pass1_genes()
        self.pass2_transcripts()
        if self.unmapped:
            top = ", ".join(f"{t} (x{n})" for t, n in self.unmapped.most_common(20))
            raise SystemExit(
                f"FATAL: {len(self.unmapped)} biotype term(s) in {self.gtf} have no "
                f"mapping for source={self.source}: {top}\n"
                f"Add each to the {self.source} map in scripts/normalize_annotation.py "
                f"with a deliberate target (or PASSTHROUGH). Refusing to guess."
            )


def merge_intervals(rows: list[tuple[str, int, int, str, str]]) -> list[tuple]:
    """Strand-aware interval merge, so no bedtools dependency."""
    out = []
    by = defaultdict(list)
    for chrom, s, e, name, strand in rows:
        by[(chrom, strand)].append((s, e, name))
    for (chrom, strand), iv in sorted(by.items()):
        iv.sort()
        cs, ce, names = iv[0][0], iv[0][1], {iv[0][2]}
        for s, e, n in iv[1:]:
            if s <= ce:
                ce = max(ce, e)
                names.add(n)
            else:
                out.append((chrom, cs, ce, ",".join(sorted(names)), 0, strand))
                cs, ce, names = s, e, {n}
        out.append((chrom, cs, ce, ",".join(sorted(names)), 0, strand))
    return out


def write_outputs(nz: Normalizer, out_dir: Path, min_ncrna: int) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    tx = {t: r for t, r in nz.tx.items() if t not in nz.dup_tx}
    ncrna = sorted(t for t, r in tx.items() if r["transcript_type"] in nz.drop_set)
    post = [r for r in tx.values() if r["transcript_type"] in POSTFILTER]
    depl = [r for r in tx.values() if r["transcript_type"] in nz.drop_set]

    (out_dir / "ncrna_tx.txt").write_text("".join(f"{t}\n" for t in ncrna))
    with (out_dir / "tx_to_gene.tsv").open("w") as fh:
        for t in sorted(tx):
            fh.write(f"{t}\t{tx[t]['gene_id']}\n")
    with (out_dir / "tx2biotype.tsv").open("w") as fh:
        fh.write("tx_id\tgene_id\tgene_name\tchrom\tstrand\ttranscript_type\tgene_type\tlength\n")
        for t in sorted(tx):
            r = tx[t]
            fh.write(f"{t}\t{r['gene_id']}\t{r['gene_name']}\t{r['chrom']}\t{r['strand']}"
                     f"\t{r['transcript_type']}\t{r['gene_type']}\t{r['length']}\n")
    for name, recs in (("ncrna_deplete.bed", depl), ("ncrna_postfilter.bed", post)):
        merged = merge_intervals([(r["chrom"], r["start"] - 1, r["end"],
                                   r["gene_id"], r["strand"]) for r in recs])
        with (out_dir / name).open("w") as fh:
            for row in merged:
                fh.write("\t".join(str(x) for x in row) + "\n")

    comp = Counter(r["transcript_type"] for r in tx.values())
    report = {
        "species": nz.species, "source_convention": nz.source,
        "gtf": str(nz.gtf), "mito_contig": nz.mito,
        "mito_features_retyped": nz.retyped_mito,
        "feature_column_rows_repaired_to_CDS": nz.repaired_feat,
        "dropped_codon_outside_exons": len(nz.bad_codon),
        "dropped_codon_outside_exons_detail": dict(list(nz.bad_codon.items())[:20]),
        "dropped_duplicate_tx_id": len(nz.dup_tx),
        "dropped_duplicate_tx_ids": sorted(nz.dup_tx)[:50],
        "dropped_invalid_strand": len(nz.bad_strand),
        "dropped_invalid_strand_ids": sorted(nz.bad_strand)[:50],
        "drop_set": sorted(nz.drop_set),
        "n_transcripts": len(tx), "n_genes": len({r["gene_id"] for r in tx.values()}),
        "n_ncrna_tx": len(ncrna),
        "n_deplete_bed_intervals": sum(1 for _ in (out_dir / "ncrna_deplete.bed").open()),
        "n_postfilter_bed_intervals": sum(1 for _ in (out_dir / "ncrna_postfilter.bed").open()),
        "biotype_composition": dict(comp.most_common()),
        "vocabulary_mapping": dict(nz.mapped.most_common()),
        "ncrna_composition": dict(Counter(
            tx[t]["transcript_type"] for t in ncrna).most_common()),
    }
    (out_dir / "normalization_report.json").write_text(json.dumps(report, indent=2) + "\n")

    # Hard gates. An empty or implausible drop list must kill the job: that is
    # the whole failure mode this script exists to prevent.
    problems = []
    if not tx:
        problems.append("no transcripts parsed at all")
    if len(ncrna) < min_ncrna:
        problems.append(f"ncrna_tx.txt has {len(ncrna)} entries, below the floor of {min_ncrna}")
    if report["n_deplete_bed_intervals"] == 0:
        problems.append("ncrna_deplete.bed is empty")
    if report["n_postfilter_bed_intervals"] == 0:
        problems.append("ncrna_postfilter.bed is empty")
    if problems:
        raise SystemExit("FATAL: " + "; ".join(problems) +
                         f"\nSee {out_dir / 'normalization_report.json'}")
    return report


EMPTY_TID_RE = re.compile(r'\s*transcript_id\s+"";')


def write_gtf(nz: Normalizer, dest: Path) -> int:
    """Re-emit the GTF with GENCODE-style attributes added, originals preserved.

    Three normalizations happen here, each fixing a real downstream breakage:

    1. `transcript_id ""` is stripped. RefSeq puts a literal empty transcript_id on
       every `gene` line (6,477 of them in yeast alone); GENCODE omits the attribute
       entirely. gffread reads the empty string as the record's ID and dies with
       "Error: no valid ID found for GFF record" -- on the RAW RefSeq GTF, before
       this script touches it. Without this, no transcriptome FASTA can be built for
       any RefSeq species, and therefore no salmon index.
    2. GENCODE-style `gene_type` / `transcript_type` are appended.
    3. FlyBase transcript-level features are emitted a second time as `transcript`,
       because RiboCode's prepare_transcripts and gffread both look for that feature
       name and FlyBase writes `mRNA` instead.
    """
    n = 0
    with opener(nz.gtf) as src, dest.open("w") as out:
        for line in src:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                out.write(line)
                continue
            f[2] = repair_feature(f)[0]
            a = attrs(f[8])
            if a.get("transcript_id") == "":
                f[8] = EMPTY_TID_RE.sub("", f[8])
                a.pop("transcript_id", None)
            tid = a.get("transcript_id")
            # A transcript dropped for an invalid strand must take its exon/CDS/UTR
            # children with it, or the GTF is left with orphan features.
            if tid and (tid in nz.bad_strand or tid in nz.dup_tx or tid in nz.bad_codon):
                continue
            rec = nz.tx.get(tid) if tid else None
            if rec is None:
                out.write("\t".join(f) + "\n")
                continue
            extra = (f' gene_type "{rec["gene_type"]}";'
                     f' transcript_type "{rec["transcript_type"]}";')
            f[8] = f[8].rstrip() + extra
            # FlyBase writes `mRNA` (and `ncRNA`, `tRNA`, ...) where a GTF normally
            # writes `transcript`, and RiboCode's prepare_transcripts and gffread both
            # look for `transcript`. RENAME it rather than emitting a second copy: the
            # biotype it carried is preserved in the transcript_type attribute added
            # above, and duplicating the line made gffread extract every fly transcript
            # twice (71,107 sequences for 35,736 transcripts).
            if nz.source == "flybase" and f[2] in nz.txfeat:
                f[2] = "transcript"
            out.write("\t".join(f) + "\n")
            n += 1
    return n


def selftest() -> int:
    """Reproduce the hand-built human and mouse ground truth, exactly.

    A normalizer that cannot regenerate the two lists already on disk has no
    business being pointed at a species where no ground truth exists.
    """
    import tempfile
    cases = [
        ("human_v49", ANNOT / "gencode.v49.annotation.gtf",
         DATA / "ncrna_tx_human_v49.txt", DATA / "tx_to_gene_human_v49.tsv"),
        ("mouse_vM38", ANNOT / "gencode.vM38.annotation.gtf",
         DATA / "ncrna_tx_mouse_vM38.txt", DATA / "tx_to_gene_mouse_vM38.tsv"),
    ]
    rc = 0
    for label, gtf, want_nc, want_t2g in cases:
        if not gtf.exists():
            print(f"[{label}] SKIP, missing {gtf}")
            continue
        print(f"[{label}] normalizing {gtf.name} ...", flush=True)
        src = detect_source(gtf)
        mito = detect_mito(gtf, None)
        nz = Normalizer(gtf, label, src, mito, CORE_DROP)
        nz.run()
        with tempfile.TemporaryDirectory() as td:
            rep = write_outputs(nz, Path(td), min_ncrna=1)
            got_nc = {l.strip() for l in (Path(td) / "ncrna_tx.txt").read_text().splitlines() if l.strip()}
            got_t2g = {tuple(l.split("\t")) for l in
                       (Path(td) / "tx_to_gene.tsv").read_text().splitlines() if l}
        exp_nc = {l.strip() for l in want_nc.read_text().splitlines() if l.strip()}
        exp_t2g = {tuple(l.split("\t")) for l in want_t2g.read_text().splitlines() if l}
        ok_src = src == "gencode"
        ok_nc = got_nc == exp_nc
        ok_t2g = got_t2g == exp_t2g
        print(f"  source detected      : {src} {'OK' if ok_src else 'FAIL (expected gencode)'}")
        print(f"  mito contig          : {mito}")
        print(f"  transcripts          : {rep['n_transcripts']}")
        print(f"  ncrna_tx             : got {len(got_nc)}, expected {len(exp_nc)}, "
              f"{'EXACT MATCH' if ok_nc else 'MISMATCH'}")
        if not ok_nc:
            print(f"    missing {len(exp_nc - got_nc)}, extra {len(got_nc - exp_nc)}")
            print(f"    e.g. missing {sorted(exp_nc - got_nc)[:5]} extra {sorted(got_nc - exp_nc)[:5]}")
        print(f"  tx_to_gene           : got {len(got_t2g)}, expected {len(exp_t2g)}, "
              f"{'EXACT MATCH' if ok_t2g else 'MISMATCH'}")
        if not ok_t2g:
            print(f"    missing {len(exp_t2g - got_t2g)}, extra {len(got_t2g - exp_t2g)}")
            print(f"    e.g. missing {sorted(exp_t2g - got_t2g)[:3]} extra {sorted(got_t2g - exp_t2g)[:3]}")
        if not (ok_src and ok_nc and ok_t2g):
            rc = 1
    print("\nSELFTEST", "PASS" if rc == 0 else "FAIL")
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true",
                    help="reproduce the on-disk human and mouse ground truth and exit")
    ap.add_argument("--gtf", type=Path)
    ap.add_argument("--genome", type=Path, help="genome FASTA, used to find the mito contig")
    ap.add_argument("--species")
    ap.add_argument("--source", default="auto",
                    choices=["auto", "gencode", "refseq", "ensembl", "flybase"])
    ap.add_argument("--drop-set", default="auto", choices=["auto", "core", "extended"])
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--write-gtf", action="store_true",
                    help="also emit <species>.normalized.gtf (needed by STAR/RiboCode/gffread)")
    ap.add_argument("--min-ncrna", type=int, default=10,
                    help="hard floor on ncrna_tx.txt; the job fails below it")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if not (args.gtf and args.species):
        ap.error("--gtf and --species are required unless --selftest")

    src = detect_source(args.gtf) if args.source == "auto" else args.source
    # A --genome path that does not exist must not degrade quietly to "no mito contig":
    # an unmatched shell glob is passed through as the literal pattern, detect_mito then
    # finds no header to grep, and the only symptom is 24 worm tRNAs silently keeping
    # `tRNA` instead of `Mt_tRNA`.
    if args.genome is not None and not args.genome.exists():
        raise SystemExit(f"ABORT: --genome {args.genome} does not exist")
    mito = detect_mito(args.gtf, args.genome)
    if args.drop_set == "auto":
        drop = CORE_DROP if src == "gencode" else EXTENDED_DROP
    else:
        drop = CORE_DROP if args.drop_set == "core" else EXTENDED_DROP
    out_dir = args.out_dir or (DATA / "annot" / args.species)

    print(f"species={args.species} source={src} mito={mito} "
          f"drop_set={'core' if drop is CORE_DROP else 'extended'}", flush=True)
    if src != "gencode" and mito is None:
        print("  WARNING: no mitochondrial contig identified. RefSeq has no Mt_rRNA/Mt_tRNA "
              "biotypes, so mito rRNA/tRNA will NOT be re-typed and may escape the drop set.",
              file=sys.stderr)

    nz = Normalizer(args.gtf, args.species, src, mito, drop)
    nz.run()
    rep = write_outputs(nz, out_dir, args.min_ncrna)

    print(f"  transcripts {rep['n_transcripts']}  genes {rep['n_genes']}")
    print(f"  ncrna_tx {rep['n_ncrna_tx']}  "
          f"({', '.join(f'{k}={v}' for k, v in rep['ncrna_composition'].items())})")
    print(f"  mito features re-typed: {rep['mito_features_retyped']}")
    if rep["dropped_codon_outside_exons"]:
        print(f"  dropped {rep['dropped_codon_outside_exons']} transcript(s) whose start/stop codon "
              f"is not contained in their exons (incomplete mito stop codons): "
              f"{', '.join(list(rep['dropped_codon_outside_exons_detail'])[:5])}")
    if rep["dropped_duplicate_tx_id"]:
        print(f"  dropped {rep['dropped_duplicate_tx_id']} transcript id(s) reused at more than "
              f"one locus: {', '.join(rep['dropped_duplicate_tx_ids'][:5])}"
              f"{' ...' if rep['dropped_duplicate_tx_id'] > 5 else ''}")
    if rep["dropped_invalid_strand"]:
        print(f"  dropped {rep['dropped_invalid_strand']} transcript(s) with strand not in +/-: "
              f"{', '.join(rep['dropped_invalid_strand_ids'][:5])}"
              f"{' ...' if rep['dropped_invalid_strand'] > 5 else ''}")
    print(f"  deplete.bed {rep['n_deplete_bed_intervals']}  "
          f"postfilter.bed {rep['n_postfilter_bed_intervals']}")
    if args.write_gtf:
        dest = out_dir / f"{args.species}.normalized.gtf"
        n = write_gtf(nz, dest)
        print(f"  wrote {dest} ({n} annotated feature lines)")
    print(f"  -> {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
