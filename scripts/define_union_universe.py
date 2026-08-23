#!/usr/bin/env python3
"""Define the UNION expressed pc+lncRNA universe across the 8 Chothani tissues (noBrain) + write its
transcript FASTA -- the broader counterpart to the Fibroblast-only define_universe_and_fasta.py.

Universe = protein_coding + lncRNA, length <= 10000, chrom != chrM, salmon meanTPM >= 1 in ANY
of the 8 tissues (Fibroblast, VSMC, ES, Fat, HA_EC, HCAEC, HUVEC, Hepatocytes). Per-tissue mean
TPM is aggregated from the archived decoy-aware Chothani salmon quant (no realignment). Fixes the
Fibroblast-only bias: the Hepatocytes holdout now includes liver-specific tx silent in Fibroblast.

Outputs (data/):
  union_universe.tsv  (tx_id, type, gene_id, gene_name, chrom, length, max_tpm, n_tissues)
  union_universe.fa   (bare versioned-ENST headers, one sequence per line)
"""

# ============================================================================
# CANNOT RUN AS-IS -- input lost 2026-08-15.
#
# SAL/SALMON below points at expression_context_human/data/salmon_quant_chothani_decoy, which went
# with the deleted biotype_probe tree. There is NO surviving copy: no salmon quant anywhere under
# data/ holds the Chothani SRR15513* accessions (data/tpm/salmon_human_rna is a different cohort).
# It was deliberately NOT repointed, because there is nothing to point at.
#
# This script's OUTPUT survives, so nothing downstream is blocked today. What is lost is the
# ability to REBUILD that output from source. Re-deriving it means re-running salmon on the 69
# Chothani RNA libraries; the BAMs are still in data/rnaseq_bam_mm1/.
# ============================================================================
import collections
import glob
import os
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
           "riboseq_signal_model")
ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
           "biotype_probe/expression_context_human")
ANN = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations")
SAL = ECH / "data" / "salmon_quant_chothani_decoy"
PC_FA = ANN / "gencode.v49.pc_transcripts.fa"
LNC_FA = ANN / "gencode.v49.lncRNA_transcripts.fa"
OUT_TSV = NEW / "data" / "union_universe.tsv"
OUT_FA = NEW / "data" / "union_universe.fa"
TISSUES = ["Fibroblast", "VSMC", "ES", "Fat", "HA_EC", "HCAEC", "HUVEC", "Hepatocytes"]
PCL = {"protein_coding", "lncRNA"}
MIN_TPM = 1.0
MAX_LEN = 10000
EXCLUDE_CHROMS = {"chrM"}


def sample_to_tissue():
    s2t = {}
    for ln in open(NEW / "data" / "loto_rnaseq_srr_tissue.tsv").readlines()[1:]:
        c = ln.rstrip("\n").split("\t")
        s2t[c[0]] = c[1]
    for ln in open(NEW / "data" / "fibroblast_rnaseq_srr.txt"):
        srr = ln.strip().split()[0] if ln.strip() else ""
        if srr:
            s2t.setdefault(srr, "Fibroblast")
    return s2t


def per_tissue_tpm():
    s2t = sample_to_tissue()
    tsamp = collections.defaultdict(list)
    for q in glob.glob(str(SAL / "*" / "quant.sf")):
        t = s2t.get(os.path.basename(os.path.dirname(q)))
        if t in TISSUES:
            tsamp[t].append(q)
    tissue_tpm = {}
    for t in TISSUES:
        acc = collections.defaultdict(float)
        n = len(tsamp[t])
        for q in tsamp[t]:
            with open(q) as f:
                f.readline()
                for ln in f:
                    p = ln.split("\t")
                    acc[p[0].split("|")[0]] += float(p[3])   # tx_id = Name before '|', TPM col 3
        tissue_tpm[t] = {tx: v / n for tx, v in acc.items()} if n else {}
    return tissue_tpm


def load_meta():
    meta = {}
    with open(NEW / "data" / "target" / "Fibroblast_psites_summary.tsv") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = {x: i for i, x in enumerate(hdr)}
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            meta[c[0]] = (c[ci["transcript_type"]], c[ci["gene_id"]], c[ci["gene_name"]],
                          c[ci["chrom"]], int(c[ci["length"]]))
    return meta


def stream_fa(path, want, out):
    written = 0
    keep = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                tx = line[1:].split("|")[0]
                keep = tx if tx in want else None
                if keep:
                    out.write(f">{keep}\n")
                    written += 1
            elif keep:
                out.write(line if line.endswith("\n") else line + "\n")
    return written


def main():
    tissue_tpm = per_tissue_tpm()
    meta = load_meta()
    allt = set().union(*[set(tissue_tpm[t]) for t in TISSUES])
    rows = []
    want = set()
    for tx in allt:
        m = meta.get(tx)
        if not m or m[0] not in PCL or m[4] > MAX_LEN or m[3] in EXCLUDE_CHROMS:
            continue
        tpms = [tissue_tpm[t].get(tx, 0.0) for t in TISSUES]
        nt = sum(1 for v in tpms if v >= MIN_TPM)
        if nt == 0:
            continue
        want.add(tx)
        rows.append((tx, m[0], m[1], m[2], m[3], m[4], max(tpms), nt))
    rows.sort()
    with open(OUT_TSV, "w") as o:
        o.write("tx_id\ttranscript_type\tgene_id\tgene_name\tchrom\tlength\tmax_tpm\tn_tissues\n")
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")
    print(f"union universe transcripts: {len(want):,}", file=sys.stderr)
    with open(OUT_FA, "w") as o:
        w1 = stream_fa(PC_FA, want, o)
        w2 = stream_fa(LNC_FA, want, o)
    print(f"FASTA sequences written: {w1 + w2:,} (pc {w1:,} + lncRNA {w2:,})", file=sys.stderr)
    if w1 + w2 != len(want):
        print(f"WARN {len(want)-(w1+w2)} universe tx not in GENCODE FASTAs", file=sys.stderr)


if __name__ == "__main__":
    main()
