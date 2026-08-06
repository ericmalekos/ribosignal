#!/usr/bin/env python3
"""Q3 payoff figure: what start-context weights did the V3 model LEARN, and do they match the
empirical Kozak PWM (V2) and the hand-picked heuristic (V0)?

V3 (`--learn_start_context`) adds a `Conv1d(4 -> 1, kernel_size=10)` over the one-hot slice; its
weight is a (1, 4, 10) tensor = a learned position x nucleotide matrix. With pad (6, 3) the 10 taps
read transcript positions [i-6 .. i+3] relative to a start codon whose A sits at index i, i.e. Kozak
positions [-6,-5,-4,-3,-2,-1,+1,+2,+3,+4] (taps 6,7,8 are the codon A,T,G themselves; the empirical
PWM omits them since START_W already captures the codon). One-hot channel order is A,C,G,T.

This extracts that kernel from best.pt, lines it up against:
  - the empirical PWM  (data/kozak_pwm.json, log2-odds at {-6..-1,+4})
  - the heuristic      (Kozak = 0.5*[purine at -3] + 0.5*[G at +4]; nonzero only at -3 A/G and +4 G)
and (a) writes a 3-panel position x nucleotide heatmap figure, (b) reports the Pearson correlation
between the learned kernel and the empirical PWM over the 7 shared context positions (per-position
mean-centered, since a per-position constant is absorbed by the sigmoid bias), and (c) prints the
learned -3 and +4 nucleotide preferences so the reader can see directly whether the model rediscovered
"-3 purine, +4 G" without being told. cas12a env (numpy + matplotlib + torch-cpu).

Usage: plot_kozak_weights.py [--run <v3_run_dir>] [--pwm <kozak_pwm.json>] [--out <png>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
NT = "ACGT"
# The 10 conv taps, in order, correspond to these Kozak positions (A of the start codon = +1):
KERNEL_POS = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4]
CONTEXT_POS = [-6, -5, -4, -3, -2, -1, 4]   # positions the empirical PWM defines (codon excluded)


def load_learned_kernel(run: Path):
    """Return (W[4,10] numpy, bias float) from the V3 run's best.pt, or None if unavailable."""
    bp = run / "best.pt"
    if not bp.exists():
        return None
    sd = torch.load(bp, map_location="cpu")
    wkey = next((k for k in sd if k.endswith("start_ctx_conv.weight")), None)
    if wkey is None:
        return None
    w = sd[wkey].numpy().reshape(4, 10)          # (out=1, in=4, k=10) -> (4 nt, 10 pos)
    bkey = next((k for k in sd if k.endswith("start_ctx_conv.bias")), None)
    b = float(sd[bkey].numpy().reshape(-1)[0]) if bkey else 0.0
    return w, b


def pwm_matrix(pwm_data):
    """Empirical PWM as a (4, 10) matrix aligned to KERNEL_POS; codon columns (+1,+2,+3) = NaN."""
    m = np.full((4, 10), np.nan, dtype=np.float64)
    pwm = pwm_data["pwm"]
    for j, p in enumerate(KERNEL_POS):
        key = str(p)
        if key in pwm:
            for r, n in enumerate(NT):
                m[r, j] = pwm[key][n]
    return m


def heuristic_matrix():
    """The V0 heuristic as a (4, 10) matrix: kz = 0.5*[purine@-3] + 0.5*[G@+4]. Nonzero only at
    -3 (A,G = 0.5) and +4 (G = 0.5); codon columns NaN to match the PWM panel's masking."""
    m = np.zeros((4, 10), dtype=np.float64)
    for j, p in enumerate(KERNEL_POS):
        if p in (1, 2, 3):
            m[:, j] = np.nan
    j3 = KERNEL_POS.index(-3)
    m[NT.index("A"), j3] = 0.5
    m[NT.index("G"), j3] = 0.5
    j4 = KERNEL_POS.index(4)
    m[NT.index("G"), j4] = 0.5
    return m


def center_cols(m):
    """Mean-center each column (position) over the 4 nucleotides, ignoring NaN columns."""
    out = m.copy()
    for j in range(out.shape[1]):
        col = out[:, j]
        if np.all(np.isfinite(col)):
            out[:, j] = col - col.mean()
    return out


def panel(ax, m, title, vlim):
    disp = np.ma.masked_invalid(m)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="0.85")
    im = ax.imshow(disp, aspect="auto", cmap=cmap, vmin=-vlim, vmax=vlim)
    ax.set_yticks(range(4))
    ax.set_yticklabels(list(NT))
    ax.set_xticks(range(10))
    ax.set_xticklabels([f"{p:+d}" for p in KERNEL_POS])
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Kozak position (start codon A = +1)")
    for j in range(10):
        for r in range(4):
            v = m[r, j]
            if np.isfinite(v):
                ax.text(j, r, f"{v:+.2f}", ha="center", va="center", fontsize=6,
                        color="black" if abs(v) < 0.6 * vlim else "white")
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(NEW / "results/kozak/v3_learned_onehot_fib2hep"))
    ap.add_argument("--pwm", default=str(NEW / "data/kozak_pwm.json"))
    ap.add_argument("--out", default=str(NEW / "figures/kozak/learned_vs_empirical_vs_heuristic.png"))
    args = ap.parse_args()

    run = Path(args.run)
    pwm_data = json.loads(Path(args.pwm).read_text())
    emp = pwm_matrix(pwm_data)
    heur = heuristic_matrix()

    lk = load_learned_kernel(run)
    if lk is None:
        print(f"[warn] no learned kernel at {run}/best.pt yet -- run after V3 training lands; "
              f"plotting empirical + heuristic only.")
        learned = np.full((4, 10), np.nan)
        bias = float("nan")
    else:
        learned, bias = lk

    # Correlation of learned vs empirical over the 7 shared context positions (per-position centered).
    ctx_idx = [KERNEL_POS.index(p) for p in CONTEXT_POS]
    lc = center_cols(learned)[:, ctx_idx].ravel()
    ec = center_cols(emp)[:, ctx_idx].ravel()
    ok = np.isfinite(lc) & np.isfinite(ec)
    if ok.sum() >= 2 and np.std(lc[ok]) > 0 and np.std(ec[ok]) > 0:
        r = float(np.corrcoef(lc[ok], ec[ok])[0, 1])
    else:
        r = float("nan")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 8.0), constrained_layout=True)
    # Learned panel uses its own scale (raw conv weights); PWM/heuristic share the log-odds scale.
    lv = np.nanmax(np.abs(learned)) if np.isfinite(learned).any() else 1.0
    ev = max(np.nanmax(np.abs(emp)), np.nanmax(np.abs(heur)), 1e-6)
    panel(axes[0], learned, f"V3 LEARNED start-context kernel (Conv1d 4x10)  bias={bias:+.2f}", lv)
    panel(axes[1], emp, "V2 EMPIRICAL Kozak PWM (log2-odds, annotated CDS ATG starts)", ev)
    panel(axes[2], heur, "V0 HEURISTIC  kz = 0.5*[purine@-3] + 0.5*[G@+4]", ev)
    fig.suptitle(f"Start-context weights: learned vs empirical vs heuristic\n"
                 f"learned-vs-empirical Pearson r = {r:.3f} over 7 context positions "
                 f"(per-position centered)", fontsize=12)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

    # Direct readout: did the model rediscover "-3 purine, +4 G"?
    if np.isfinite(learned).any():
        def col(pos):
            return learned[:, KERNEL_POS.index(pos)]
        print("\nLearned nucleotide weights (raw conv), key positions:")
        for pos in (-3, 4):
            c = col(pos)
            order = ", ".join(f"{NT[i]}:{c[i]:+.3f}" for i in np.argsort(-c))
            print(f"  Kozak {pos:+d}: {order}   (top = most start-promoting)")
        pur3 = col(-3)[NT.index("A")] + col(-3)[NT.index("G")]
        pyr3 = col(-3)[NT.index("C")] + col(-3)[NT.index("T")]
        g4 = col(4)[NT.index("G")]
        print(f"  -> -3 purine(A+G) sum {pur3:+.3f} vs pyrimidine(C+T) sum {pyr3:+.3f} "
              f"(heuristic predicts purine > pyrimidine)")
        print(f"  -> +4 G weight {g4:+.3f} vs mean of others "
              f"{np.mean([col(4)[i] for i in range(4) if i != NT.index('G')]):+.3f} "
              f"(heuristic predicts G highest)")
    print(f"\nlearned-vs-empirical Pearson r = {r:.3f} (7 context positions, per-position centered)")

    # Persist a compact summary for aggregate_kozak.py / results.md.
    summ = {
        "run": run.name,
        "learned_available": bool(np.isfinite(learned).any()),
        "learned_vs_empirical_pearson_r": r,
        "bias": bias,
        "kernel_pos": KERNEL_POS,
        "nt_order": list(NT),
        "learned_kernel_4x10": (learned.tolist() if np.isfinite(learned).any() else None),
    }
    sp = Path(args.out).with_suffix(".json")
    sp.write_text(json.dumps(summ, indent=2))
    print(f"wrote {sp}")


if __name__ == "__main__":
    main()
