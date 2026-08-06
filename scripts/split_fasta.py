#!/usr/bin/env python3
"""Split a FASTA into N chunks balanced by total sequence length (greedy longest-first).

Balancing by nt (not sequence count) evens the GPU wall-time across chunks, since RiNALMo
per-token extraction time scales with sequence length (L^2 attention).

Usage: split_fasta.py <input.fa> <out_dir> <n_chunks>  ->  <out_dir>/chunk_{1..N}.fa
"""
import heapq
import sys
from pathlib import Path


def main():
    inp, outdir, n = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
    outdir.mkdir(parents=True, exist_ok=True)

    recs = []
    hdr, seq = None, []
    with open(inp) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr:
                    recs.append((hdr, "".join(seq)))
                hdr, seq = line.rstrip("\n"), []
            else:
                seq.append(line.rstrip("\n"))
    if hdr:
        recs.append((hdr, "".join(seq)))

    recs.sort(key=lambda r: -len(r[1]))
    heap = [(0, i) for i in range(n)]        # (running nt, chunk index)
    heapq.heapify(heap)
    chunks = [[] for _ in range(n)]
    for hdr, s in recs:
        tot, idx = heapq.heappop(heap)
        chunks[idx].append((hdr, s))
        heapq.heappush(heap, (tot + len(s), idx))

    for i in range(n):
        with open(outdir / f"chunk_{i + 1}.fa", "w") as o:
            for hdr, s in chunks[i]:
                o.write(f"{hdr}\n{s}\n")
        nt = sum(len(s) for _, s in chunks[i])
        print(f"chunk_{i + 1}: {len(chunks[i]):,} seqs  {nt:,} nt", file=sys.stderr)


if __name__ == "__main__":
    main()
