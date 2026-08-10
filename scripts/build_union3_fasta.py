#!/usr/bin/env python3
"""Write the mouse-liver UNION3 universe FASTA from the GENCODE vM38 pc + lncRNA transcript FASTAs.

The 3x3 factorial scores all three datasets in one shared universe
(data/heldout_refs/mouse_liver_union3_tx.txt, 22,974 transcripts). That universe had a tx list but no
FASTA, and both the ORF track and the model's one-hot input are built from the FASTA.

Header style matches every other universe FASTA in data/heldout_refs/: a bare versioned transcript id
(`>ENSMUST00000130201.8`), sequence unwrapped onto one line. GENCODE ships pipe-delimited headers, so
field 0 is taken and the rest discarded.

HARD-FAILS if any requested transcript is absent. Emitting a FASTA that silently covers a subset is
the documented failure mode here: dump_pred_profiles drops transcripts missing from the one-hot FASTA
WITHOUT erroring, so a short FASTA turns into a biased evaluation subset rather than a crash
(feedback_dump_silent_skip_onehot_fasta).
"""
import argparse
import pathlib
import sys

ANN = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations")


def read_fasta_subset(path, want, out, seen):
    """Stream one GENCODE FASTA, writing any wanted record. Returns count written."""
    n = 0
    keep = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    out.write(f">{keep}\n{''.join(buf)}\n")
                    seen.add(keep)
                    n += 1
                tx = line[1:].split("|")[0].strip()
                keep = tx if (tx in want and tx not in seen) else None
                buf = []
            elif keep:
                buf.append(line.strip())
        if keep:
            out.write(f">{keep}\n{''.join(buf)}\n")
            seen.add(keep)
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tx-list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fastas", nargs="+",
                    default=[str(ANN / "gencode.vM38.pc_transcripts.fa"),
                             str(ANN / "gencode.vM38.lncRNA_transcripts.fa")])
    a = ap.parse_args()

    want = [l.strip() for l in open(a.tx_list) if l.strip()]
    wset = set(want)
    if len(wset) != len(want):
        print(f"note: tx list has {len(want)} lines, {len(wset)} unique", file=sys.stderr)

    for f in a.fastas:
        if not pathlib.Path(f).exists():
            print(f"missing source FASTA: {f}", file=sys.stderr)
            return 1

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    with open(out, "w") as fh:
        for f in a.fastas:
            n = read_fasta_subset(f, wset, fh, seen)
            print(f"  {pathlib.Path(f).name}: {n:,} sequences")

    missing = wset - seen
    print(f"  wrote {len(seen):,} of {len(wset):,} requested -> {out}")
    if missing:
        ex = ", ".join(sorted(missing)[:8])
        print(f"FAIL: {len(missing):,} transcripts absent from the source FASTAs: {ex}",
              file=sys.stderr)
        return 1
    print("  OK: every requested transcript present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
