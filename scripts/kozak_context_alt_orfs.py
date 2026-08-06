#!/usr/bin/env python3
"""Quick test: is the Kozak start context different for CTG (non-AUG) ORFs than for ATG ORFs?

The empirical PWM (data/kozak_pwm.json) was fit on canonical annotated CDS ATG starts. V2 applies that
same PWM to CTG candidate starts. This checks whether that is even the right prior: it takes RiboCode's
real CTG-aware Hepatocytes ORF calls (the dropin_ctg real_collapsed.txt), splits them by start codon, and
compares the observed start context (position frequencies + -3 purine / +4 G shares + the annotated-PWM
score distribution) between ATG and CTG ORFs, against the annotated-CDS reference. cas12a env (numpy)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/"
           "expression_context_human")
DEFAULT_CALLS = NEW / "results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/dropin_ctg/real_collapsed.txt"
FASTA = ECH / "data/ribocode_annot/transcripts_sequence.fa"
PWMJSON = NEW / "data/kozak_pwm.json"
POSITIONS = [-6, -5, -4, -3, -2, -1, 4]
NT = "ACGT"


def rel_index(a, pos):
    return a + (pos - 1) if pos > 0 else a + pos


def load_seqs(fasta, need):
    seqs, cur, keep, buf = {}, None, False, []
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None and keep:
                    seqs[cur] = "".join(buf)
                cur = ln[1:].split()[0]; keep = cur in need; buf = []
            elif keep:
                buf.append(ln.strip())
        if cur is not None and keep:
            seqs[cur] = "".join(buf)
    return seqs


def pwm_from_counts(counts, bg):
    pwm = {}
    for p in POSITIONS:
        tot = sum(counts[p].values()) + 4.0
        pwm[p] = {n: float(np.log2(((counts[p][n] + 1.0) / tot) / bg[n])) for n in NT}
    return pwm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default=str(DEFAULT_CALLS))
    ap.add_argument("--min_len", type=int, default=90, help="drop ORFs shorter than this (nt)")
    ap.add_argument("--orf_types", default="",
                    help="comma list of ORF_type to keep (e.g. uORF,Overlap_uORF to control for 5'UTR "
                         "location); empty = all types")
    args = ap.parse_args()
    keep_types = {t for t in args.orf_types.split(",") if t} or None

    pj = json.loads(Path(PWMJSON).read_text())
    bg = pj["background"]
    ann_pwm = {int(p): pj["pwm"][p] for p in pj["pwm"]}
    # annotated reference position frequencies, backed out of PWM + background (freq = 2^logodds * bg)
    ann_freq = {}
    for p in POSITIONS:
        f = {n: (2 ** ann_pwm[p][n]) * bg[n] for n in NT}
        s = sum(f.values()); ann_freq[p] = {n: f[n] / s for n in NT}

    rows = []
    with open(args.calls) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            try:
                rows.append((f[ci["transcript_id"]], int(f[ci["ORF_tstart"]]),
                             int(f[ci["ORF_length"]]), f[ci["ORF_type"]]))
            except (ValueError, KeyError):
                continue
    need = {tx for tx, _, _, _ in rows}
    seqs = load_seqs(FASTA, need)

    groups = {c: {"counts": {p: {n: 0 for n in NT} for p in POSITIONS}, "scores": [], "n": 0}
              for c in ("ATG", "CTG")}
    for tx, tstart, olen, otype in rows:
        if olen < args.min_len:
            continue
        if keep_types is not None and otype not in keep_types:
            continue
        s = seqs.get(tx)
        if s is None:
            continue
        su = s.upper().replace("U", "T")
        a = tstart - 1
        if a < 6 or a + 3 >= len(su):
            continue
        cod = su[a:a + 3]
        if cod not in groups:
            continue
        g = groups[cod]; g["n"] += 1
        sc = 0.0
        for p in POSITIONS:
            nt = su[rel_index(a, p)]
            if nt in NT:
                g["counts"][p][nt] += 1
                sc += ann_pwm[p][nt]
        g["scores"].append(sc)

    def pur3(cnt):
        t = sum(cnt[-3].values()); return (cnt[-3]["A"] + cnt[-3]["G"]) / t if t else float("nan")

    def g4(cnt):
        t = sum(cnt[4].values()); return cnt[4]["G"] / t if t else float("nan")

    print(f"\ncalls file: {args.calls}")
    print(f"min ORF length: {args.min_len} nt\n")
    ann_pur3 = ann_freq[-3]["A"] + ann_freq[-3]["G"]
    ann_g4 = ann_freq[4]["G"]
    ann_score_med = pj["calibration"]["p5"], pj["calibration"]["p95"]
    print(f"{'group':<22} {'n':>7} {'-3 purine':>10} {'+4 G':>7} {'median annPWM score':>20}")
    print(f"{'annotated CDS (ATG)':<22} {pj['n_annotated_starts']:>7} {ann_pur3:>10.3f} "
          f"{ann_g4:>7.3f} {'(p50~1.55 by def)':>20}")
    for c in ("ATG", "CTG"):
        g = groups[c]
        med = float(np.median(g["scores"])) if g["scores"] else float("nan")
        print(f"{'RiboCode ' + c + ' ORFs':<22} {g['n']:>7} {pur3(g['counts']):>10.3f} "
              f"{g4(g['counts']):>7.3f} {med:>20.2f}")

    # per-position consensus log-odds: annotated vs CTG group
    ctg_pwm = pwm_from_counts(groups["CTG"]["counts"], bg)
    atg_pwm = pwm_from_counts(groups["ATG"]["counts"], bg)
    print("\nper-position log2-odds (best nt), annotated CDS | RiboCode-ATG | RiboCode-CTG:")
    print(f"{'pos':>4}   {'annotated':>16}   {'ribo-ATG':>16}   {'ribo-CTG':>16}")
    for p in POSITIONS:
        ab = max(NT, key=lambda n: ann_pwm[p][n])
        tb = max(NT, key=lambda n: atg_pwm[p][n])
        cb = max(NT, key=lambda n: ctg_pwm[p][n])
        print(f"{p:>+4}   {ab}:{ann_pwm[p][ab]:>+5.2f} ({ann_freq[p][ab]:.2f})   "
              f"{tb}:{atg_pwm[p][tb]:>+5.2f}          {cb}:{ctg_pwm[p][cb]:>+5.2f}")

    # how many CTG ORFs would the annotated PWM push toward the floor vs the top?
    if groups["CTG"]["scores"]:
        cs = np.array(groups["CTG"]["scores"])
        p5, p95 = pj["calibration"]["p5"], pj["calibration"]["p95"]
        below = float((cs < p5).mean()); above = float((cs > p95).mean())
        print(f"\nCTG ORFs under the ANNOTATED-PWM calibration: {below:.1%} below p5 (-> floor 0.5x), "
              f"{above:.1%} above p95 (-> 1.0x); ATG ORFs below p5: "
              f"{float((np.array(groups['ATG']['scores'])<p5).mean()):.1%}")


if __name__ == "__main__":
    main()
