#!/usr/bin/env python3
"""Sequence-derived, annotation-free ORF-candidate track for the Fibroblast universe.

Per nucleotide, channels computed from the mature mRNA sequence ALONE (no GTF/CDS), so they
are defined identically for protein_coding CDS, 5'UTR uORFs, 3'UTR dORFs, and lncRNA ORFs.
The model gets the candidate ORF landscape everywhere and must still learn which candidates
are actually translated from the FM embedding + Ribo-seq signal. This is the goal-aligned
alternative to an annotated-CDS track, which would mark only the canonical ORF and blind the
model to the non-canonical translation this project targets.

Two modes (--mode), differing only in the start channel so a model A/B isolates the value of
non-AUG initiation information:

  atg  (v1)  channel 3 = is_atg   : binary, first nt of any AUG codon
  ext  (v2)  channel 3 = start_ext: graded start propensity over AUG + the 9 near-cognate
             starts 1 nt from AUG (CUG/GUG/UUG/ACG/AUA/AUU/AUC/AAG/AGG), weighted by relative
             initiation efficiency (AUG=1, CUG strongest non-AUG). uORFs / dORFs / lncRNA ORFs
             are enriched for these non-AUG starts, so 'ext' is the non-canonical-aware variant.
             The optional Kozak context factor (--kozak) DEFAULTS TO none as of results.md Task 20:
             the ablation showed the hand-picked heuristic is redundant with what the sequence
             backbone learns on its own (a learnable gate rediscovers the Kozak PWM at r=0.807) and
             mildly hurt held-out non-AUG detection, so the default start channel is now codon
             weight only. Pass --kozak heuristic to reproduce the old v2 track.

Channels (both modes):
  0 orf_f0   nt is inside an ATG..in-frame-stop stretch, ATG start pos %3 == 0
  1 orf_f1   ... start pos %3 == 1
  2 orf_f2   ... start pos %3 == 2
  3 is_atg | start_ext   (see modes above)
  4 is_stop  first nt of any stop codon TAA/TAG/TGA

The ATG-based occupancy (channels 0-2) is kept identical across modes so it stays an
informative canonical prior (~37% per frame) rather than washing out; non-AUG initiation
enters only through the graded start channel, which avoids the density explosion of opening
occupancy at 10 start codons. Kozak = 0.5*[purine at -3] + 0.5*[G at +4] of the start codon.
Packed in the SAME tx order + offsets as data/packed/coverage.npy so the loader slices it
exactly like coverage.

Output (default): data/packed/orf_track.npy (atg, int8) or orf_track_v2_<kozak>.npy
(ext, float16, where <kozak> is nokozak / kozak / kozakpwm),
plus a sibling _meta.json.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths  # noqa: E402

# Relative initiation efficiency of AUG + near-cognate starts (T-alphabet), ordinal, from the
# non-AUG initiation literature (CUG strongest non-AUG, then GUG/ACG, then the rest).
START_W = {"ATG": 1.0, "CTG": 0.5, "GTG": 0.35, "ACG": 0.35, "TTG": 0.3,
           "ATA": 0.25, "ATT": 0.2, "ATC": 0.2, "AAG": 0.15, "AGG": 0.15}
A, C, G, T = 65, 67, 71, 84


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


def _codon(c0, c1, c2, cod):
    return (c0 == ord(cod[0])) & (c1 == ord(cod[1])) & (c2 == ord(cod[2]))


def _pwm_multiplier(b, pwm_data, L):
    """Per-candidate-start Kozak multiplier in [floor, 1.0] from an empirical PWM (build_kozak_pwm.py).
    b = uint8 ASCII sequence (length L). Returns a length-(L-2) array aligned to codon-start positions:
    the value at index i is the context multiplier for a start whose A sits at i. Out-of-range context
    positions contribute 0 to the score (partial context near transcript ends), mirroring the heuristic."""
    nt_idx = np.full(L, -1, dtype=np.int64)
    nt_idx[b == A] = 0
    nt_idx[b == C] = 1
    nt_idx[b == G] = 2
    nt_idx[b == T] = 3
    positions = pwm_data["positions"]
    tab = np.array([[pwm_data["pwm"][str(p)][n] for n in "ACGT"] for p in positions], dtype=np.float32)
    score = np.zeros(L, dtype=np.float32)
    for k, p in enumerate(positions):
        off = (p - 1) if p > 0 else p                 # transcript-index offset of Kozak pos p from the A
        src = np.full(L, -1, dtype=np.int64)
        i_lo, i_hi = max(0, -off), min(L, L - off)
        ii = np.arange(i_lo, i_hi)
        src[ii] = nt_idx[ii + off]
        valid = src >= 0
        score[valid] += tab[k][src[valid]]
    cal = pwm_data["calibration"]
    p5, p95, floor = cal["p5"], cal["p95"], cal.get("floor", 0.5)
    m = np.clip(floor + (1.0 - floor) * (score - p5) / (p95 - p5), floor, 1.0)
    return m[:L - 2]


def orf_track(seq, mode, kozak="none", pwm_data=None):
    """(L, 5) ORF-candidate channels from a mature-mRNA sequence."""
    s = seq.upper().replace("U", "T")
    L = len(s)
    dt = np.int8 if mode == "atg" else np.float16
    out = np.zeros((L, 5), dtype=dt)
    if L < 3:
        return out
    b = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    c0, c1, c2 = b[:-2], b[1:-1], b[2:]                       # codon-start-aligned triples
    atg = _codon(c0, c1, c2, "ATG")
    ctg = _codon(c0, c1, c2, "CTG")
    stop = (_codon(c0, c1, c2, "TAA") | _codon(c0, c1, c2, "TAG")
            | _codon(c0, c1, c2, "TGA"))
    out[:L - 2, 4] = stop.view(np.int8) if mode == "atg" else stop.astype(np.float16)

    # start channel (channel 3) and the codons that open an ORF for the occupancy channels.
    start_mask = atg                                          # atg / ext: ATG-only occupancy
    if mode == "atg":
        out[:L - 2, 3] = atg.view(np.int8)
    elif mode == "atgctg":
        # 2-level start: ATG=1.0, CTG=0.5, and occupancy EXTENDED to CTG-initiated ORFs (v3).
        # Tests just the single dominant alternative start (CTG) rather than ext's full 10-codon
        # near-cognate set, whose dilution dropped the ATG lncRNA lift (+0.056 -> +0.033).
        st = np.where(atg, 1.0, 0.0) + np.where(ctg, 0.5, 0.0)
        out[:L - 2, 3] = st.astype(np.float16)
        start_mask = atg | ctg
    else:  # ext (v2): graded start over AUG + 9 near-cognates, modulated by the Kozak context.
        w = np.zeros(L - 2, dtype=np.float32)          # ATG-only occupancy (channels 0-2) unchanged
        for cod, wt in START_W.items():
            w[_codon(c0, c1, c2, cod)] = wt
        if kozak == "none":
            mult = 1.0                                 # codon weight only, no context (Kozak removed)
        elif kozak == "pwm":
            mult = _pwm_multiplier(b, pwm_data, L)      # empirical PWM multiplier (length L-2)
        else:                                          # heuristic: 0.5*[purine@-3] + 0.5*[G@+4]
            purine = (b == A) | (b == G)
            isg = (b == G)
            kz = np.zeros(L - 2, dtype=np.float32)
            idx = np.arange(L - 2)
            m3, p3 = idx - 3, idx + 3
            vm, vp = m3 >= 0, p3 < L
            kz[vm] += 0.5 * purine[m3[vm]]
            kz[vp] += 0.5 * isg[p3[vp]]
            mult = 0.5 + 0.5 * kz
        out[:L - 2, 3] = (w * mult).astype(np.float16)

    # ORF occupancy per frame: open at a start codon, close at the next in-frame stop.
    # start_mask = {ATG} for atg/ext, {ATG,CTG} for atgctg.
    for f in range(3):
        open_start = None
        for i in range(f, L - 2, 3):
            if start_mask[i] and open_start is None:
                open_start = i
            elif stop[i] and open_start is not None:
                out[open_start:i + 3, f] = 1
                open_start = None
        if open_start is not None:
            out[open_start:L, f] = 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["atg", "ext", "atgctg"], default="atg")
    ap.add_argument("--out", default=None, help="output .npy (default per mode)")
    ap.add_argument("--pack", default=None,
                    help="pack dir (offsets/lengths/tx_order) to build the track for "
                         "(default $RIBO_PACK_DIR, else <repo>/data/packed)")
    ap.add_argument("--fasta", default=None,
                    help="universe FASTA matching --pack's tx (default $RIBO_ONEHOT_FASTA)")
    ap.add_argument("--kozak", choices=["none", "heuristic", "pwm"], default="none",
                    help="ext-mode start-channel context model: none (codon weight only, Kozak removed; "
                         "the DEFAULT since results.md Task 20 -- the Kozak factor is redundant with what "
                         "the sequence backbone learns, r=0.807, and mildly hurt held-out non-AUG "
                         "detection), heuristic (0.5*[purine@-3]+0.5*[G@+4], the old v2 default, kept for "
                         "reproducibility), or pwm (empirical Kozak PWM from --pwm_file). Ignored for "
                         "atg/atgctg modes.")
    ap.add_argument("--pwm_file", default=None, help="kozak_pwm.json (required for --kozak pwm)")
    args = ap.parse_args()
    pack = Path(args.pack) if args.pack else paths.pack_dir()
    fasta = Path(args.fasta) if args.fasta else paths.onehot_fasta()
    if not (pack / "tx_order.txt").is_file():
        print(paths.missing(pack, "pack dir", "RIBO_PACK_DIR", "pack"), file=sys.stderr)
        return 1
    if not fasta.is_file():
        print(paths.missing(fasta, "universe FASTA", "RIBO_ONEHOT_FASTA", "fasta"), file=sys.stderr)
        return 1
    pwm_data = None
    if args.kozak == "pwm":
        if not args.pwm_file:
            print("ERROR: --kozak pwm requires --pwm_file", file=sys.stderr)
            return 1
        pwm_data = json.loads(Path(args.pwm_file).read_text())
    # The default name states the START-CONTEXT MODEL, not just the mode. It used to be a bare
    # "orf_track_v2.npy" for every --kozak setting, so a `--kozak none` build and a
    # `--kozak heuristic` build wrote the SAME filename with different content, and telling them
    # apart meant opening the sidecar meta. Feeding a nokozak-trained model the heuristic track
    # is a silent accuracy loss, not an error.
    _kz = {"none": "_nokozak", "heuristic": "_kozak", "pwm": "_kozakpwm"}[args.kozak]
    default_name = {"atg": "orf_track.npy",
                    "ext": f"orf_track_v2{_kz}.npy",
                    "atgctg": "orf_track_v3.npy"}[args.mode]
    out_path = Path(args.out) if args.out else pack / default_name

    offsets = np.load(pack / "offsets.npy")
    lengths = np.load(pack / "lengths.npy")
    order = (pack / "tx_order.txt").read_text().split()
    sum_L = int(offsets[-1])
    print(f"mode={args.mode}  tx: {len(order):,}  sum_L: {sum_L:,}  -> {out_path}",
          file=sys.stderr)

    seqs = read_fasta(fasta)
    missing = [t for t in order if t not in seqs]
    if missing:
        print(f"ERROR: {len(missing)} packed tx missing from FASTA, e.g. {missing[:3]}",
              file=sys.stderr)
        return 1

    dt = np.int8 if args.mode == "atg" else np.float16
    packed = np.zeros((sum_L, 5), dtype=dt)
    t0 = time.time()
    for r, tx in enumerate(order):
        a, b = int(offsets[r]), int(offsets[r + 1])
        L = b - a
        if L != int(lengths[r]):
            print(f"ERROR: {tx} offset len {L} != lengths {int(lengths[r])}", file=sys.stderr)
            return 1
        trk = orf_track(seqs[tx], args.mode, args.kozak, pwm_data)
        if trk.shape[0] != L:
            print(f"ERROR: {tx} seq len {trk.shape[0]} != packed len {L}", file=sys.stderr)
            return 1
        packed[a:b] = trk
        if (r + 1) % 10000 == 0:
            print(f"  [{r+1:,}/{len(order):,}] {time.time()-t0:.0f}s", file=sys.stderr)

    np.save(out_path, packed)
    # float64 accumulator: a float32 sum of ~1e8 values saturates at 2^24 and corrupts the means
    frac = packed.mean(axis=0, dtype=np.float64)
    ch3_name = {"atg": "is_atg", "ext": "start_ext", "atgctg": "start_atgctg"}[args.mode]
    start_w = {"atg": "binary is_atg", "ext": START_W,
               "atgctg": {"ATG": 1.0, "CTG": 0.5}}[args.mode]
    occ_starts = {"atg": "ATG", "ext": "ATG", "atgctg": "ATG|CTG"}[args.mode]
    chan_names = ["orf_f0", "orf_f1", "orf_f2", ch3_name, "is_stop"]
    meta = {
        "mode": args.mode,
        "kozak": args.kozak if args.mode == "ext" else "n/a",
        "pwm_file": args.pwm_file if args.kozak == "pwm" else None,
        "channels": chan_names,
        "occupancy_start_codons": occ_starts,
        "shape": list(packed.shape),
        "dtype": str(packed.dtype),
        "start_weights": start_w,
        # key channel 3 by its mode-specific name (is_atg / start_ext / start_atgctg)
        "mean_per_channel": {c: float(frac[i]) for i, c in enumerate(chan_names)},
    }
    meta_path = out_path.with_name(out_path.stem + "_meta.json")
    meta_path.write_text(json.dumps(meta, indent=2))
    print(f"wrote {out_path} {packed.shape} {packed.dtype} ({packed.nbytes/1e6:.0f} MB); "
          f"means {meta['mean_per_channel']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
