#!/usr/bin/env python3
"""Write the transcript FASTA for a universe, subset from GENCODE transcript FASTAs.

`define_universe_and_fasta.py` does this too, but it is Fibroblast-hardcoded with no CLI, so any
other universe needed a source edit. This is the same operation exposed as a tool.

GENCODE transcript FASTA headers are pipe-delimited (`>ENSMUST00000000001.5|ENSMUSG...|...`); the
universe list carries bare versioned tx-ids. The id is taken as the first pipe-field, and the output
header is that bare id so the FASTA is directly usable as `RIBO_ONEHOT_FASTA`.

Missing transcripts are a hard error by default. A universe with no sequence for some of its rows
silently yields a biased subset downstream rather than a crash
(feedback_dump_silent_skip_onehot_fasta), so it must be caught here.
"""
import argparse
import pathlib
import sys


def read_fasta(path):
    """Yield (id, [lines]) from a GENCODE transcript FASTA, keyed by the first pipe-field."""
    tid, buf = None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if tid is not None:
                    yield tid, buf
                tid, buf = ln[1:].split("|")[0].split()[0], []
            else:
                buf.append(ln.rstrip("\n"))
    if tid is not None:
        yield tid, buf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe-tx", required=True, help="one versioned tx-id per line; sets order")
    ap.add_argument("--source-fasta", required=True, nargs="+",
                    help="GENCODE transcript FASTA(s), e.g. pc + lncRNA")
    ap.add_argument("--out", required=True)
    ap.add_argument("--allow-missing", action="store_true",
                    help="warn instead of failing when a universe tx has no sequence")
    a = ap.parse_args()

    order = [l.strip() for l in open(a.universe_tx) if l.strip()]
    want = set(order)
    seqs = {}
    for src in a.source_fasta:
        for tid, buf in read_fasta(src):
            if tid in want and tid not in seqs:
                seqs[tid] = buf

    missing = [t for t in order if t not in seqs]
    if missing:
        msg = (f"{len(missing):,} of {len(order):,} universe tx have no sequence in "
               f"{a.source_fasta} (first: {missing[:3]})")
        if not a.allow_missing:
            sys.exit(f"ERROR: {msg}\nPass --allow-missing only if a biased subset is acceptable.")
        print(f"WARN: {msg}", file=sys.stderr)

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out, "w") as fh:
        for t in order:  # universe order, so the FASTA is row-aligned with the pack
            if t not in seqs:
                continue
            fh.write(f">{t}\n")
            fh.write("\n".join(seqs[t]) + "\n")
            n += 1
    print(f"wrote {out}  ({n:,} sequences of {len(order):,} universe tx)")


if __name__ == "__main__":
    main()
