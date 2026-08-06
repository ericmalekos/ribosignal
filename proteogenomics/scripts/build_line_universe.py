#!/usr/bin/env python3
"""Build a cell-line's expressed pc+lncRNA universe (tx list + FASTA) from its salmon quant, mirroring the
Fibroblast define_universe_and_fasta rule but parameterized for any line: TPM >= MIN_TPM, transcript_type in
{protein_coding, lncRNA}, mature length <= MAX_LEN (10000, the model's cap), chrom != chrM. Writes
<out_prefix>_universe_tx.txt (versioned ENST, one per line) + <out_prefix>_universe.fa (headers = versioned
ENST, matching the tx list). Reuses tx2biotype.tsv + the GENCODE transcript FASTAs. cas12a env (stdlib)."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ANN = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations")
REF = {
    "human": {
        "tx2bt": NEW / "data" / "tx2biotype.tsv",
        "fastas": [ANN / "gencode.v49.pc_transcripts.fa", ANN / "gencode.v49.lncRNA_transcripts.fa"],
    },
    "mouse": {
        "tx2bt": NEW / "data" / "tx2biotype_mouse.tsv",
        "fastas": [ANN / "gencode.vM38.pc_transcripts.fa", ANN / "gencode.vM38.lncRNA_transcripts.fa"],
    },
}
KEEP_BT = {"protein_coding", "lncRNA"}


def iter_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = ln[1:].split("|")[0].split()[0].strip(); buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salmon", required=True, help="path to quant.sf")
    ap.add_argument("--out_prefix", required=True, help="e.g. .../DoHH2_pilot/DoHH2 -> _universe_tx.txt/.fa")
    ap.add_argument("--species", default="human", choices=["human", "mouse"])
    ap.add_argument("--min_tpm", type=float, default=1.0)
    ap.add_argument("--max_len", type=int, default=10000)
    args = ap.parse_args()
    ref = REF[args.species]

    # tx2biotype: tx_id -> (chrom, transcript_type, length)
    meta = {}
    with open(ref["tx2bt"]) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            meta[p[ci["tx_id"]]] = (p[ci["chrom"]], p[ci["transcript_type"]], int(p[ci["length"]]))

    # expressed tx (TPM >= min): salmon Name = "ENST..|ENSG..|.." -> versioned ENST
    expressed = set()
    with open(args.salmon) as f:
        f.readline()
        for ln in f:
            p = ln.split("\t")
            if float(p[3]) >= args.min_tpm:
                expressed.add(p[0].split("|")[0])
    # apply universe filters
    universe = set()
    for tx in expressed:
        m = meta.get(tx)
        if m and m[1] in KEEP_BT and m[0] != "chrM" and m[2] <= args.max_len:
            universe.add(tx)

    out_tx = Path(f"{args.out_prefix}_universe_tx.txt")
    out_tx.write_text("\n".join(sorted(universe)) + "\n")
    # FASTA: extract universe tx from the transcriptome FASTAs, header = versioned ENST
    out_fa = Path(f"{args.out_prefix}_universe.fa")
    n_fa = 0
    with open(out_fa, "w") as w:
        for fa in ref["fastas"]:
            if not Path(fa).exists():
                continue
            for name, seq in iter_fasta(fa):
                if name in universe:
                    w.write(f">{name}\n{seq}\n"); n_fa += 1
    print(f"expressed(TPM>={args.min_tpm}): {len(expressed):,}  universe(pc/lncRNA, <= {args.max_len} nt, "
          f"no chrM): {len(universe):,}  FASTA written: {n_fa:,}")
    print(f"wrote {out_tx}\nwrote {out_fa}")
    if n_fa != len(universe):
        print(f"WARNING: {len(universe)-n_fa} universe tx missing from FASTA (id mismatch?)")


if __name__ == "__main__":
    main()
