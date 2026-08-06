#!/usr/bin/env python3
"""Empirically fit a Kozak start-context PWM from annotated CDS start codons, to replace the hand-picked
heuristic (0.5*[purine at -3] + 0.5*[G at +4]) in the ORF-track start channel (build_orf_track.py ext mode).

Source: every transcript in the Fibroblast universe with an annotated start codon (tx2cds.tsv
has_start_codon=1; the start A sits at index utr5_len in the mature-mRNA sequence). Context positions
scored = {-6,-5,-4,-3,-2,-1,+4} relative to the A (the ATG triplet at +1/+2/+3 is constant and captured by
the per-codon START_W, so it is excluded). Position-frequency matrix -> log2 odds vs the universe
background nucleotide frequency (pseudocount 1). A candidate start's context score is the summed log-odds;
the applied multiplier (used by build_orf_track --kozak pwm) is score min-max'd onto [0.5, 1.0] using the
p5/p95 of the annotated-start score distribution (so V2 matches the heuristic's [0.5,1] dynamic range).

No test-label leakage: this uses annotation SEQUENCE only, never the Hepatocytes Ribo-seq signal -- the
same posture as the literature-derived heuristic it replaces. cas12a env (numpy). Writes data/kozak_pwm.json
and prints the PWM + a direct check of the heuristic's assumptions (-3 purine share, +4 G share)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
FASTA = NEW / "data/fibroblast_universe.fa"
TX2CDS = NEW / "data/tx2cds.tsv"
OUT = NEW / "data/kozak_pwm.json"
POSITIONS = [-6, -5, -4, -3, -2, -1, 4]     # relative to the start A (=+1); ATG (+1,+2,+3) excluded
NT = "ACGT"
PC = 1.0                                     # pseudocount


def read_fasta(path):
    seqs, h, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(parts)
                h = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if h is not None:
        seqs[h] = "".join(parts)
    return seqs


def rel_index(a, pos):
    """Transcript index of Kozak position `pos` (A = +1, no 0) for a start A at index a."""
    return a + (pos - 1) if pos > 0 else a + pos      # +4 -> a+3 ; -3 -> a-3


def main():
    seqs = read_fasta(FASTA)
    print(f"universe transcripts: {len(seqs):,}")

    # annotated starts: utr5_len = 0-based index of the start A
    starts = {}
    with open(TX2CDS) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[ci["has_start_codon"]] != "1":
                continue
            starts[f[ci["tx_id"]]] = int(f[ci["utr5_len"]])

    # background nt frequency over all universe sequence
    bgc = {n: 0 for n in NT}
    for s in seqs.values():
        su = s.upper().replace("U", "T")
        for n in NT:
            bgc[n] += su.count(n)
    bgtot = sum(bgc.values())
    bg = {n: bgc[n] / bgtot for n in NT}

    # position-frequency counts over annotated starts
    counts = {p: {n: 0 for n in NT} for p in POSITIONS}
    n_used = n_skip_bounds = n_skip_notatg = 0
    for tx, a in starts.items():
        s = seqs.get(tx)
        if s is None:
            continue
        su = s.upper().replace("U", "T")
        if a < 6 or a + 3 >= len(su):        # need -6 and +4 in range
            n_skip_bounds += 1
            continue
        if su[a:a + 3] != "ATG":             # verify the annotated start really is ATG
            n_skip_notatg += 1
            continue
        for p in POSITIONS:
            nt = su[rel_index(a, p)]
            if nt in counts[p]:
                counts[p][nt] += 1
        n_used += 1
    print(f"annotated ATG starts used: {n_used:,} (skipped {n_skip_bounds:,} out-of-bounds, "
          f"{n_skip_notatg:,} start != ATG)")

    # PWM = log2( freq / background )
    pwm = {}
    for p in POSITIONS:
        tot = sum(counts[p].values()) + 4 * PC
        pwm[p] = {n: float(np.log2(((counts[p][n] + PC) / tot) / bg[n])) for n in NT}

    # score every annotated start -> calibration percentiles
    scores = []
    for tx, a in starts.items():
        s = seqs.get(tx)
        if s is None:
            continue
        su = s.upper().replace("U", "T")
        if a < 6 or a + 3 >= len(su) or su[a:a + 3] != "ATG":
            continue
        scores.append(sum(pwm[p][su[rel_index(a, p)]] for p in POSITIONS
                          if su[rel_index(a, p)] in NT))
    scores = np.array(scores)
    p5, p95 = float(np.percentile(scores, 5)), float(np.percentile(scores, 95))

    OUT.write_text(json.dumps({
        "positions": POSITIONS, "nt_order": list(NT), "pseudocount": PC,
        "background": bg, "pwm": {str(p): pwm[p] for p in POSITIONS},
        "calibration": {"p5": p5, "p95": p95, "floor": 0.5,
                        "note": "m = clip(0.5 + 0.5*(score-p5)/(p95-p5), 0.5, 1.0)"},
        "n_annotated_starts": n_used, "source": "annotated ATG starts, Fibroblast universe",
    }, indent=2))

    # ---- report ----
    print(f"\nbackground nt freq: " + "  ".join(f"{n}={bg[n]:.3f}" for n in NT))
    print(f"annotated-start context score: p5={p5:.2f}  p50={np.percentile(scores,50):.2f}  p95={p95:.2f}")
    print("\nPWM (log2 odds vs background); most-favored nt per position in [] :")
    print(f"{'pos':>4} " + "  ".join(f"{n:>6}" for n in NT) + "   consensus")
    for p in POSITIONS:
        best = max(NT, key=lambda n: pwm[p][n])
        print(f"{p:>+4} " + "  ".join(f"{pwm[p][n]:>6.2f}" for n in NT) + f"     [{best}]")

    # direct check of the heuristic's two assumptions
    tot3 = sum(counts[-3].values())
    pur3 = (counts[-3]["A"] + counts[-3]["G"]) / tot3
    tot4 = sum(counts[4].values())
    g4 = counts[4]["G"] / tot4
    print(f"\nHEURISTIC CHECK (what the hand-picked rule assumed):")
    print(f"  -3 purine (A|G) share: {pur3:.3f}   (heuristic rewards purine here)")
    print(f"  +4 G share:            {g4:.3f}   (heuristic rewards G here)")
    print(f"  strongest single position by |log-odds| range: "
          f"{max(POSITIONS, key=lambda p: max(pwm[p].values())-min(pwm[p].values())):+d}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
