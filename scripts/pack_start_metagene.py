#!/usr/bin/env python3
"""Rebuild the start-codon metagene from a pack's stored P-site counts.

The original RiboCode `metaplots` configs for the training packs do not survive (no log before
2026-08-07 carries them), so the P-site offsets those packs were built with cannot be read off
directly. They can be recovered from the packs themselves: `target_counts.npy` holds P-site counts
with the offsets ALREADY APPLIED, so the metagene profile around annotated start codons says where
those offsets actually put the P-sites.

If a pack's offsets were right, the initiating ribosome's P-site sits on the AUG and the metagene
peaks at position 0. A peak off zero means that pack was built with a P-site offset that mislocates
the ribosome by exactly that many nucleotides -- which is a property of the stored training data, not
of any rerun.

Run it on a rebuilt pack too, and the two profiles are the old-vs-new metaplot comparison.
"""
import argparse
import json
import pathlib

import numpy as np


def load_cds(path):
    """tx_id -> (cds_start, cds_end), 1-based inclusive as RiboCode writes them."""
    cds = {}
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                cds[f[0]] = (int(f[1]), int(f[2]))
            except ValueError:
                continue
    return cds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--counts", default="target_counts.npy",
                    help="array inside --pack to profile (default the stored P-sites)")
    ap.add_argument("--window", type=int, default=60)
    ap.add_argument("--min-psites", type=int, default=50,
                    help="skip transcripts below this, so the profile is not dominated by noise")
    ap.add_argument("--label", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    pack = pathlib.Path(a.pack)
    label = a.label or pack.name
    order = [l.strip() for l in open(pack / "tx_order.txt")]
    lengths = np.load(pack / "lengths.npy").astype(np.int64)
    offsets = np.concatenate([[0], np.cumsum(lengths)])
    counts = np.load(pack / a.counts, mmap_mode="r")
    cds = load_cds(a.cds)

    W = a.window
    prof = np.zeros(2 * W + 1, dtype=np.float64)
    n_used = 0
    for i, tx in enumerate(order):
        c = cds.get(tx)
        if c is None:
            continue
        s = c[0] - 1  # to 0-based index of the A in the start codon
        lo, hi = offsets[i], offsets[i + 1]
        if s - W < 0 or s + W + 1 > (hi - lo):
            continue
        seg = np.asarray(counts[lo + s - W: lo + s + W + 1], dtype=np.float64)
        tot = seg.sum()
        if tot < a.min_psites:
            continue
        prof += seg / tot  # normalize per transcript so deep genes do not dominate the peak
        n_used += 1

    if n_used == 0:
        raise SystemExit(f"{label}: no transcripts passed the filters")
    prof /= n_used
    pos = np.arange(-W, W + 1)
    peak = int(pos[int(np.argmax(prof))])
    # Frame phasing: which of the three offsets-mod-3 carries the most signal downstream of the start.
    downstream = prof[W:]
    frames = [float(downstream[f::3].sum()) for f in range(3)]
    best_frame = int(np.argmax(frames))

    res = dict(label=label, pack=str(pack), counts=a.counts, n_transcripts=n_used,
               peak_offset=peak, peak_height=float(prof.max()),
               frame_fractions=[f / sum(frames) for f in frames], dominant_frame=best_frame)
    print(f"  {label:<28} n={n_used:>6}  peak at {peak:+d} nt  "
          f"frame0={res['frame_fractions'][0]:.3f} f1={res['frame_fractions'][1]:.3f} "
          f"f2={res['frame_fractions'][2]:.3f}")
    if a.out:
        p = pathlib.Path(a.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        np.save(p.with_suffix(".npy"), prof)
        p.with_suffix(".json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
