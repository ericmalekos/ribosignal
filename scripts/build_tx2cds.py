#!/usr/bin/env python3
"""Transcript-coordinate 5'UTR / CDS / 3'UTR lengths from the GENCODE v49 GTF.

Needed for the uORF-aware evaluation: uORFs are the most common alternative ORF, and the
whole-transcript profile Pearson is dominated by the main CDS, so it barely reflects 5'UTR
(uORF) prediction. This table lets eval restrict the profile Pearson to the 5'UTR window
[0, cds_start) (and the 3'UTR window for dORFs). Training remains annotation-free; these CDS
coordinates are used ONLY to define eval windows, never fed to the model.

GENCODE v49 annotates undifferentiated `UTR` (not five_prime/three_prime), so 5' vs 3' is
resolved by strand relative to the CDS extent:
  + strand: 5'UTR = exonic bases with genomic pos < min(CDS start); 3'UTR = pos > max(CDS end)
  - strand: 5'UTR = exonic bases with genomic pos > max(CDS end); 3'UTR = pos < min(CDS start)
The first CDS base (transcript orientation) coincides with the start codon for 5'-complete CDS,
so cds_start_tx == utr5_len. tx_len is the summed exon length (matches the mature-transcript L
the embeddings/target are built on). has_start_codon flags 5'-complete models (reliable 5'UTR).

Output: data/tx2cds.tsv with tx_id, strand, tx_len, utr5_len, cds_len, utr3_len, has_start_codon.
Only transcripts with a CDS (coding) are emitted.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

GTF = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "annotations/gencode.v49.annotation.gtf")
OUT = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model/data/tx2cds.tsv")
# --gtf / --out override the defaults so a cross-species held-out (mouse vM38) gets its own
# CDS table with the identical parse logic.
_ap = argparse.ArgumentParser()
_ap.add_argument("--gtf", default=None)
_ap.add_argument("--out", default=None)
_ARGS, _ = _ap.parse_known_args()
if _ARGS.gtf:
    GTF = Path(_ARGS.gtf)
if _ARGS.out:
    OUT = Path(_ARGS.out)


def txid(attr: str) -> str:
    i = attr.find('transcript_id "')
    if i < 0:
        return ""
    i += len('transcript_id "')
    return attr[i:attr.find('"', i)]


def overlap(s, e, a, b):
    """length of [s,e] cap [a,b], all 1-based inclusive."""
    return max(0, min(e, b) - max(s, a) + 1)


def main():
    exons: dict[str, list[tuple[int, int]]] = {}
    strand: dict[str, str] = {}
    cds_lo: dict[str, int] = {}
    cds_hi: dict[str, int] = {}
    has_sc: dict[str, bool] = {}

    t0 = time.time()
    n = 0
    with GTF.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            feat = f[2]
            if feat not in ("exon", "CDS", "start_codon"):
                continue
            s, e, st, attr = int(f[3]), int(f[4]), f[6], f[8]
            t = txid(attr)
            if not t:
                continue
            if feat == "exon":
                exons.setdefault(t, []).append((s, e))
                strand[t] = st
            elif feat == "CDS":
                cds_lo[t] = min(cds_lo.get(t, s), s)
                cds_hi[t] = max(cds_hi.get(t, e), e)
            else:  # start_codon
                has_sc[t] = True
            n += 1
            if n % 5_000_000 == 0:
                print(f"  {n:,} feature lines  {time.time()-t0:.0f}s", file=sys.stderr)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with OUT.open("w") as out:
        out.write("tx_id\tstrand\ttx_len\tutr5_len\tcds_len\tutr3_len\thas_start_codon\n")
        for t, ex in exons.items():
            if t not in cds_lo:                       # non-coding: no CDS, skip
                continue
            lo, hi, stv = cds_lo[t], cds_hi[t], strand[t]
            tx_len = sum(e - s + 1 for s, e in ex)
            five = three = 0
            for s, e in ex:
                below = overlap(s, e, -1, lo - 1)     # genomic pos < cds_lo
                above = overlap(s, e, hi + 1, 10 ** 12)  # genomic pos > cds_hi
                if stv == "+":
                    five += below
                    three += above
                else:
                    five += above
                    three += below
            cds_len = tx_len - five - three
            if cds_len <= 0:
                continue
            out.write(f"{t}\t{stv}\t{tx_len}\t{five}\t{cds_len}\t{three}\t"
                      f"{1 if has_sc.get(t) else 0}\n")
            rows += 1

    print(f"wrote {OUT} ({rows:,} coding transcripts) in {time.time()-t0:.0f}s",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
