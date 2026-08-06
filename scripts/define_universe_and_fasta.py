#!/usr/bin/env python3
"""Define the Fibroblast expressed pc+lncRNA universe and write its transcript FASTA.

Universe = protein_coding + lncRNA transcripts with Fibroblast salmon meanTPM >= 1 and
mature length <= 10000 nt (the RiNALMo extraction cap). This is the transcript set for
per-token FM embeddings and the per-nt signal model (a training subset can be taken
later; extracting the generous set once avoids re-extraction).

Outputs:
  data/fibroblast_universe.tsv  (tx_id, transcript_type, gene_id, gene_name, chrom,
                                 length, mean_tpm, ribo_psites)
  data/fibroblast_universe.fa   (bare versioned-ENST headers, one sequence per line;
                                 the RiNALMo extractor input format)
"""
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ANN = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations")
PC_FA = ANN / "gencode.v49.pc_transcripts.fa"
LNC_FA = ANN / "gencode.v49.lncRNA_transcripts.fa"
OUT_TSV = NEW / "data" / "fibroblast_universe.tsv"
OUT_FA = NEW / "data" / "fibroblast_universe.fa"

MIN_TPM = 1.0
MAX_LEN = 10000
PCL = {"protein_coding", "lncRNA"}
# Mitochondrial mRNAs are mitoribosome-translated (different genetic code, leaderless, no
# cytoplasmic periodicity) -- not this model's target. Drop chrM so the universe / FASTA /
# embeddings never include them (matches dataset.excluded_tx, which also gates every split).
EXCLUDE_CHROMS = {"chrM"}


def load_kv(path, key_col, val_cols, cast):
    d = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[idx[key_col]]] = tuple(cast[j](f[idx[c]]) for j, c in enumerate(val_cols))
    return d


def stream_fa(path, want, out):
    written = 0
    keep = False
    cur = None
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur:
                    out.write(f">{cur}\n{''.join(seq)}\n")
                    written += 1
                cur = line[1:].split("|")[0].strip()
                keep = cur in want
                seq = []
            elif keep:
                seq.append(line.strip())
        if keep and cur:
            out.write(f">{cur}\n{''.join(seq)}\n")
            written += 1
    return written


def main():
    bt = load_kv(NEW / "data" / "tx2biotype.tsv", "tx_id",
                 ["transcript_type", "gene_id", "gene_name", "chrom", "length"],
                 [str, str, str, str, int])
    tpm = load_kv(NEW / "data" / "fibroblast_salmon_mean_tpm.tsv", "tx_id",
                  ["mean_tpm"], [float])
    ribo = load_kv(NEW / "data" / "target" / "Fibroblast_psites_summary.tsv", "tx_id",
                   ["total_psites"], [int])

    rows = []
    want = set()
    for tx, (ttype, gid, gname, chrom, L) in bt.items():
        if ttype not in PCL or L > MAX_LEN or chrom in EXCLUDE_CHROMS:
            continue
        if tpm.get(tx, (0.0,))[0] < MIN_TPM:
            continue
        want.add(tx)
        rows.append((tx, ttype, gid, gname, chrom, L,
                     tpm.get(tx, (0.0,))[0], ribo.get(tx, (0,))[0]))
    rows.sort()

    with OUT_TSV.open("w") as o:
        o.write("tx_id\ttranscript_type\tgene_id\tgene_name\tchrom\tlength\t"
                "mean_tpm\tribo_psites\n")
        for r in rows:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.4f}\t{r[7]}\n")
    print(f"universe transcripts: {len(want):,}", file=sys.stderr)

    with OUT_FA.open("w") as o:
        w1 = stream_fa(PC_FA, want, o)
        w2 = stream_fa(LNC_FA, want, o)
    print(f"FASTA sequences written: {w1 + w2:,} (pc {w1:,} + lncRNA {w2:,})", file=sys.stderr)
    missing = len(want) - (w1 + w2)
    if missing:
        print(f"WARN {missing} universe tx not found in the GENCODE transcript FASTAs",
              file=sys.stderr)


if __name__ == "__main__":
    main()
