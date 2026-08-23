#!/usr/bin/env python3
"""Positional P-site metagene around the annotated start and stop codons, per training tissue.

WHY THIS EXISTS. `S_riboseq_qc` exports per-library QC SCALARS (read-length distribution, frame-0
fraction, P-site offset per read length) but no positional profile, so there was no way to draw the
one plot every Ribo-seq audience expects: density rising at the start codon, 3-nt periodicity through
the CDS, dropping at the stop. RiboCode's `metaplots` computes those quantities internally to pick
offsets and does not persist them.

WHAT THE POSITION MEANS -- THIS IS THE FIRST THING TO STATE IN A CAPTION.
Values are **P-SITES, not read 5' ends**. The pack's `target_counts.npy` is what RiboCode's
`process_bam.py` emits after applying the per-read-length P-site offset chosen by `metaplots`, so
the offset correction is already applied and offset 0 is the **first nucleotide of the start codon**
(the A of the AUG) for the start anchor, and the **first nucleotide of the stop codon** for the stop
anchor. A 5'-end metagene would sit ~12 nt to the left of this one.

DEFINITIONS THAT CHANGE THE PICTURE, SO THEY ARE EXPLICIT:
  * TRANSCRIPT coordinates, not genomic. The packs are per-transcript, so no splicing correction is
    needed or applied.
  * Only transcripts with an annotated CDS AND an annotated start codon (`has_start_codon == 1` in
    tx2cds) contribute. A CDS inferred without a start codon would blur offset 0, which is the one
    position the figure exists to show.
  * A transcript contributes to a given offset only if that offset is inside the transcript, so
    `n_tx` VARIES WITH OFFSET and is reported per row. Dividing the pooled `psites` by a single
    global transcript count would fabricate a decay at the window edges purely from missing UTR.
  * Each transcript is counted ONCE per offset regardless of depth; `psites` is a raw sum. Highly
    translated genes therefore dominate, which is what a conventional metagene shows. `psites_norm`
    is the per-transcript-normalised alternative (each transcript scaled to sum 1 inside its own
    window before pooling), which removes that dominance -- report which one a figure uses.

Output (one row per tissue x anchor x offset):

    tissue  anchor  offset_nt  psites  n_tx  psites_norm

  python3 scripts/build_metagene.py --out results/metagene
  python3 scripts/build_metagene.py --packs data/packed_canon_Hepatocytes --tissue Hepatocytes
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib

import numpy as np

NEW = pathlib.Path(__file__).resolve().parents[1]
# The 8 canonical Chothani training packs. Brain is excluded: it was dropped from training for low
# periodicity, and including it in a "training data" metagene would misrepresent what was trained on.
DEFAULT_PACKS = sorted(p for p in (NEW / "data").glob("packed_canon*") if p.is_dir())
START_WIN = (-50, 150)      # nt relative to the first base of the start codon
STOP_WIN = (-150, 50)       # nt relative to the first base of the stop codon


def tissue_of(pack: pathlib.Path) -> str:
    """packed_canon -> Fibroblast (the base pack), packed_canon_HUVEC -> HUVEC."""
    n = pack.name
    return "Fibroblast" if n == "packed_canon" else n.replace("packed_canon_", "")


def load_cds(path):
    """-> {tx_id: (utr5_len, cds_len)} for transcripts with an annotated start codon."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("has_start_codon") != "1":
                continue
            try:
                u5, cds = int(r["utr5_len"]), int(r["cds_len"])
            except (KeyError, ValueError):
                continue
            if cds >= 3:
                out[r["tx_id"]] = (u5, cds)
    return out


def metagene(pack: pathlib.Path, cds, min_psites=0):
    """-> {anchor: (psites[W], n_tx[W], norm[W])} pooled over the pack's transcripts."""
    order = (pack / "tx_order.txt").read_text().split()
    offs = np.load(pack / "offsets.npy")
    tgt = np.load(pack / "target_counts.npy", mmap_mode="r")
    res = {}
    for anchor, (lo, hi) in (("start", START_WIN), ("stop", STOP_WIN)):
        W = hi - lo + 1
        res[anchor] = [np.zeros(W, np.float64), np.zeros(W, np.int64), np.zeros(W, np.float64)]
    n_used = 0
    for i, tx in enumerate(order):
        c = cds.get(tx)
        if c is None:
            continue
        u5, cl = c
        a, b = int(offs[i]), int(offs[i + 1])
        L = b - a
        # A CDS annotation that overruns the packed transcript means the pack and the annotation
        # disagree about this transcript; skip rather than silently clip into the wrong frame.
        if u5 + cl > L:
            continue
        v = np.asarray(tgt[a:b], dtype=np.float64)
        if min_psites and v.sum() < min_psites:
            continue
        n_used += 1
        for anchor, (lo, hi) in (("start", START_WIN), ("stop", STOP_WIN)):
            zero = u5 if anchor == "start" else u5 + cl - 3
            s, e = zero + lo, zero + hi + 1           # desired slice in transcript coords
            cs, ce = max(0, s), min(L, e)             # clipped to what exists
            if ce <= cs:
                continue
            seg = v[cs:ce]
            d0, d1 = cs - s, (cs - s) + (ce - cs)     # where it lands in the window
            acc, cnt, nrm = res[anchor]
            acc[d0:d1] += seg
            cnt[d0:d1] += 1
            t = seg.sum()
            if t > 0:
                nrm[d0:d1] += seg / t
    return res, n_used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packs", nargs="*", default=[str(p) for p in DEFAULT_PACKS])
    ap.add_argument("--tissue", default=None, help="override the tissue label (single pack only)")
    ap.add_argument("--tx2cds", default=str(NEW / "data/tx2cds.tsv"))
    ap.add_argument("--min-psites", type=int, default=0,
                    help="drop transcripts with fewer than this many P-sites (default 0 = keep all)")
    ap.add_argument("--out", default=str(NEW / "results/metagene"))
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cds = load_cds(a.tx2cds)
    print(f"  {len(cds):,} transcripts with an annotated CDS + start codon")

    rows, pooled, meta = [], {}, []
    for ps in a.packs:
        pack = pathlib.Path(ps)
        if not (pack / "target_counts.npy").exists():
            print(f"  SKIP {pack.name}: no target_counts.npy")
            continue
        t = a.tissue or tissue_of(pack)
        res, n_used = metagene(pack, cds, a.min_psites)
        tot = int(sum(r[0].sum() for r in res.values()))
        print(f"  {t:<14} {n_used:>6,} tx contribute   {tot:>14,} windowed P-sites")
        meta.append(dict(tissue=t, pack=str(pack), n_tx=n_used))
        for anchor, (lo, hi) in (("start", START_WIN), ("stop", STOP_WIN)):
            acc, cnt, nrm = res[anchor]
            for j, o in enumerate(range(lo, hi + 1)):
                rows.append(dict(tissue=t, anchor=anchor, offset_nt=o,
                                 psites=int(acc[j]), n_tx=int(cnt[j]),
                                 psites_norm=float(nrm[j])))
                k = (anchor, o)
                p = pooled.setdefault(k, [0, 0, 0.0])
                p[0] += int(acc[j]); p[1] += int(cnt[j]); p[2] += float(nrm[j])

    cols = ["tissue", "anchor", "offset_nt", "psites", "n_tx", "psites_norm"]
    with open(out / "metagene_by_tissue.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.6f}" if c == "psites_norm" else str(r[c])
                               for c in cols) + "\n")
    with open(out / "metagene_pooled.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for (anchor, o), (p, n, nr) in sorted(pooled.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            fh.write(f"pooled\t{anchor}\t{o}\t{p}\t{n}\t{nr:.6f}\n")

    # Frame bias inside the CDS is the headline QC number a metagene is usually read for.
    frames = {}
    for anchor in ("start", "stop"):
        sel = [r for r in rows if r["anchor"] == anchor and 0 <= r["offset_nt"] < 150] \
            if anchor == "start" else []
        if sel:
            f = [0.0, 0.0, 0.0]
            for r in sel:
                f[r["offset_nt"] % 3] += r["psites"]
            s = sum(f) or 1.0
            frames = {f"frame{i}": f[i] / s for i in range(3)}
    (out / "metagene_meta.json").write_text(json.dumps(dict(
        position_semantics=("P-SITE, not read 5' end. Pack targets come from RiboCode process_bam "
                            "after applying the metaplots-derived per-read-length offset."),
        anchor_zero=("start: first nt of the start codon (A of AUG). "
                     "stop: first nt of the stop codon."),
        coordinates="transcript (not genomic); no splicing correction needed",
        start_window=list(START_WIN), stop_window=list(STOP_WIN),
        transcripts="annotated CDS AND has_start_codon == 1 in data/tx2cds.tsv",
        n_tx_varies_with_offset=("yes -- a transcript contributes to an offset only if that offset "
                                 "exists in it. Divide psites by the per-row n_tx, never by a global "
                                 "count, or the window edges decay for the wrong reason."),
        psites_norm=("each transcript scaled to sum 1 within its own window before pooling; "
                     "removes the dominance of highly translated genes that `psites` has"),
        cds_frame_fraction_first150nt=frames,
        packs=meta), indent=2))
    print(f"\n  wrote {out}/metagene_by_tissue.tsv, metagene_pooled.tsv, metagene_meta.json")
    if frames:
        print(f"  CDS frame fractions (start+0..149 nt, pooled): "
              + "  ".join(f"{k}={v:.3f}" for k, v in frames.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
