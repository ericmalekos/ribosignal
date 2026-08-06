#!/usr/bin/env python3
"""Example good predictions: observed vs predicted per-nt Ribo-seq P-site profile for one uORF, one CDS,
and one lncRNA ORF from the held-out Hepatocytes test set (the model never saw Hepatocytes). Three rows,
each with a context view (whole/near-transcript, regions annotated) and a feature zoom (5'UTR uORF, or a
3-nt-periodicity window). Predictions come from the already-dumped pred_profiles.npz -- CPU only.

Defaults are the auto-selected examples (select on whole/5'UTR profile Pearson + real signal); override
with --uorf_tx / --cds_tx / --lnc_tx. cas12a env (numpy + matplotlib).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/"
           "expression_context_human")
NPZ = NEW / "results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/dropin/pred_profiles.npz"
CALLS = ECH / "data/ribocode_per_tissue/Hepatocytes/Hepatocytes_collapsed.txt"
OBS_C = "#1f77b4"
PRED_C = "#d62728"
FRAME_C = ["#1b9e77", "#d95f02", "#7570b3"]  # frame 0 (in-frame) / 1 / 2, colorblind-safe (Dark2)


def pearson(a, b):
    a = a.astype(np.float64); b = b.astype(np.float64)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def draw_periodicity(ax, lo, hi, p, q, ref_start, title, xlab=True):
    """Per-nucleotide bars of the PREDICTED P-site fraction, colored by codon frame relative to the ORF
    start ((pos - ref_start) % 3); observed overlaid as a dark step line. A well-translated region has
    frame-0 (in-frame) bars dominating -- the visual signature of 3-nt periodicity."""
    x = np.arange(lo, hi)
    frames = (x - ref_start) % 3
    ax.bar(x, p[lo:hi], width=0.9, color=[FRAME_C[f] for f in frames], linewidth=0)
    ax.step(x, q[lo:hi], where="mid", color="#111", lw=0.8)
    ax.margins(x=0); ax.set_ylabel("P-site fraction", fontsize=8)
    ax.set_title(title, fontsize=9); ax.tick_params(labelsize=7)
    if xlab:
        ax.set_xlabel("transcript position (nt)", fontsize=8)
    handles = [Patch(color=FRAME_C[f], label=f"frame {f}" + (" (in-frame)" if f == 0 else ""))
               for f in range(3)] + [Line2D([0], [0], color="#111", lw=0.8, label="observed")]
    ax.legend(handles=handles, fontsize=6.5, loc="upper right", framealpha=0.9, ncol=2)


ZOOM_W = 60          # nt shown in every right-hand zoom, so the three rows are directly comparable


def best_zoom(c, lo, hi, w=ZOOM_W):
    """Frame-aligned (start at lo+3k) w-nt window inside [lo, hi) with the most SUSTAINED signal:
    maximize the count of nonzero observed positions (breadth), tie-broken by summed signal. This
    lands on a periodic body region rather than a single dominant spike."""
    hi = min(hi, len(c))
    best, bkey = (lo, min(lo + w, hi)), (-1, -1.0)
    for s in range(lo, max(lo + 1, hi - w + 1), 3):
        e = min(s + w, hi)
        seg = c[s:e]
        key = (int((seg > 0).sum()), float(seg.sum()))
        if key > bkey:
            bkey, best = key, (s, e)
    return best


def load_tsv(path, key, cols):
    out = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            out[f[ix[key]]] = tuple(f[ix[c]] for c in cols)
    return out


def orf_windows(path):
    calls = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            calls.setdefault(f[ix["transcript_id"]], []).append(
                (f[ix["ORF_type"]], int(f[ix["ORF_tstart"]]), int(f[ix["ORF_tstop"]])))
    return calls


def region_quality(c, lo, hi, ref_start):
    """(n_nonzero, max_share, frame0_share) of the OBSERVED counts in [lo, hi).

    max_share is the largest single position's share of the region total, and frame0_share is the
    fraction of region signal in frame 0 relative to `ref_start`.
    """
    seg = c[max(lo, 0):hi]
    tot = seg.sum()
    if tot <= 0 or seg.size < 9:
        return 0, 1.0, 0.0
    x = np.arange(max(lo, 0), max(lo, 0) + seg.size)
    f0 = seg[(x - ref_start) % 3 == 0].sum()
    return int((seg > 0).sum()), float(seg.max() / tot), float(f0 / tot)


def autoselect(tx_ids, idx, off, pred, obs, bio, cds, calls, min_obs=200, min_utr_obs=50,
               min_nonzero=40, max_spike=0.15, min_frame0=0.50, min_uorf_share=0.15,
               min_nonzero_uorf=15, max_spike_uorf=0.30, max_uorf_nt=150, min_uorf_r=0.90):
    """Re-pick the best uORF / CDS / lncRNA exemplar for THIS prediction dump.

    The three transcripts a figure showcases are a property of the model that produced it, not of
    the annotation: a different checkpoint fits different transcripts best, so carrying the previous
    model's picks into a new figure would silently show examples chosen for someone else.

    Ranking on Pearson ALONE does not work, and the failure is systematic rather than unlucky. A
    transcript whose signal is one dominant spike scores r ~ 1.0 trivially, because observed and
    predicted agree on the single position that carries all the variance. Selecting on raw Pearson
    against the mamba4 dump returned exactly that: a CDS at r = 1.000 that is a lone spike over a
    flat body, and a lncRNA whose spike sits at the edge of its ORF -- both useless as
    illustrations, and both *better* by the metric than the distributed examples they replaced.

    So the gates below encode what makes a good ILLUSTRATION, and Pearson only ranks what survives:

      min_nonzero  the region must have breadth, not one position
      max_spike    no single position may carry more than this share of the region
      min_frame0   the observed signal must actually be periodic, which is the thing being shown

    uORFs get the SAME gates, only with thresholds scaled to their length (`min_nonzero_uorf`,
    `max_spike_uorf`), plus `min_uorf_share`. An earlier version exempted them from the breadth and
    spike gates on the reasoning that a uORF is short and legitimately peaky. That was wrong for a
    figure: it selected 12-codon uORFs whose entire signal is one initiation peak, so the
    frame-coloured zoom rendered as a single bar and could not show periodicity at all -- the very
    thing the panel exists to show. Requiring breadth instead selects longer uORFs (20-127 codons)
    that are genuinely periodic, of which there are 359 in this dataset, so nothing is lost.

    `min_uorf_share` is still needed on top, because 5'UTR Pearson is even more saturated than the
    whole-transcript version: ~1,961 candidates yield only 8-13 above r = 0.99 and the top few
    differ in the FOURTH decimal place. Those near-ties are not interchangeable as illustrations --
    the uORF's share of transcript signal among them runs from 46% (TBPL1) to 0.5% (EIF5A) -- and a
    0.0002 difference in r once picked PTK2 (5.5% share, visually dwarfed by an unrelated CDS
    spike) over C1orf43 (36.6%).

    -> {"uorf": (tx, r), "cds": (tx, r), "lnc": (tx, r)}
    """
    best = {"uorf": (None, -2.0), "cds": (None, -2.0), "lnc": (None, -2.0)}
    best_uorf_f0 = -1.0          # uORFs rank on periodicity, not on saturated Pearson (see below)
    for tx in tx_ids:
        j = idx[tx]
        c = obs[off[j]:off[j + 1]].astype(np.float64)
        tot = c.sum()
        if tot < min_obs:
            continue
        p = pred[off[j]:off[j + 1]].astype(np.float64)
        q = c / tot
        meta = bio.get(tx)
        if meta is None:
            continue
        _gene, btype = meta
        wins = calls.get(tx, [])
        types = {t for t, _s, _e in wins}

        if btype == "lncRNA":
            w = next((x for x in wins if x[0] in ("novel", "novel_NonPCGs")), None)
            if w is None:
                continue
            lo, hi = w[1] - 1, w[2]
            nz, spike, f0 = region_quality(c, lo, hi, lo)
            if nz < min_nonzero or spike > max_spike or f0 < min_frame0:
                continue
            r = pearson(p, q)
            if r > best["lnc"][1]:
                best["lnc"] = (tx, r)
            continue
        if btype != "protein_coding":
            continue

        cb = cds.get(tx)
        if cb is None:
            continue
        u5, cl = int(cb[0]), int(cb[1])
        cds0, cds1 = u5, u5 + cl

        if "uORF" in types and u5 > 30 and c[:u5].sum() >= min_utr_obs:
            uo = next((x for x in wins if x[0] in ("uORF", "Overlap_uORF")), None)
            if uo is not None:
                ua, ue = uo[1] - 1, uo[2]
                nz, spike, f0 = region_quality(c, ua, ue, ua)
                share = float(c[ua:ue].sum() / tot)
                # Prefer a COMPACT uORF: one long enough to be periodic but short enough that a
                # 60 nt zoom covers a real fraction of it. A 250 nt uORF whose signal sits at the
                # 3' end leaves most of the zoom panel empty.
                if (nz >= min_nonzero_uorf and spike <= max_spike_uorf and f0 >= min_frame0
                        and share >= min_uorf_share and (ue - ua) <= max_uorf_nt):
                    r = pearson(p[:u5], q[:u5])
                    # The panel draws PREDICTED bars coloured by frame, so the prediction's own
                    # periodicity is what the reader sees -- ranking on the observed side alone
                    # picked a uORF at observed f0 = 0.89 whose predicted f0 was only 0.73, and
                    # whose displayed window was 0.62 (a frame-1 bar taller than any in-frame one).
                    # Score the weaker of the two so both sides must be periodic.
                    pf0 = region_quality(p, ua, ue, ua)[2]
                    score = min(f0, pf0)
                    if np.isfinite(r) and r >= min_uorf_r and score > best_uorf_f0:
                        best_uorf_f0 = score
                        best["uorf"] = (tx, r)

        nz, spike, f0 = region_quality(c, cds0, cds1, cds0)
        if nz < min_nonzero or spike > max_spike or f0 < min_frame0:
            continue
        r = pearson(p, q)
        if r > best["cds"][1]:
            best["cds"] = (tx, r)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uorf_tx", default="ENST00000399697.7")   # POLR1D
    ap.add_argument("--cds_tx", default="ENST00000429711.7")    # RPL32
    ap.add_argument("--lnc_tx", default="ENST00000647872.1")    # LINC02693 (signal concentrated on ORF, clean frame)
    ap.add_argument("--npz", default=str(NPZ),
                    help="pred_profiles.npz to plot. The defaults above were auto-selected for the "
                         "ORIGINAL dump; point this at another model and pass --auto so the "
                         "exemplars are re-picked for it.")
    ap.add_argument("--auto", action="store_true",
                    help="re-select the best uORF / CDS / lncRNA exemplar for THIS dump")
    ap.add_argument("--title", default=None, help="override the figure suptitle")
    ap.add_argument("--out", default=str(NEW / "figures/prediction_examples/prediction_examples.png"))
    args = ap.parse_args()

    z = np.load(args.npz, allow_pickle=False)
    tx_ids = [str(t) for t in z["tx_ids"]]
    idx = {t: i for i, t in enumerate(tx_ids)}
    lengths = z["lengths"]; off = np.concatenate([[0], np.cumsum(lengths)])
    pred, obs = z["pred_flat"], z["obs_flat"]
    bio = load_tsv(NEW / "data/tx2biotype.tsv", "tx_id", ("gene_name", "transcript_type"))
    cds = load_tsv(NEW / "data/tx2cds.tsv", "tx_id", ("utr5_len", "cds_len", "utr3_len"))
    calls = orf_windows(CALLS)

    if args.auto:
        sel = autoselect(tx_ids, idx, off, pred, obs, bio, cds, calls)
        missing = [k for k, (t, _r) in sel.items() if t is None]
        if missing:
            raise SystemExit(f"auto-selection found no candidate for: {missing}")
        args.uorf_tx, args.cds_tx, args.lnc_tx = (sel["uorf"][0], sel["cds"][0], sel["lnc"][0])
        for k in ("uorf", "cds", "lnc"):
            t, r = sel[k]
            print(f"auto-selected {k:5s}: {t} {bio[t][0]:<14s} r = {r:.3f}")

    def get(tx):
        j = idx[tx]
        p = pred[off[j]:off[j + 1]].astype(np.float64)
        c = obs[off[j]:off[j + 1]].astype(np.float64)
        return p, c, c / c.sum()

    def draw(ax, lo, hi, p, q, title, xlab=True):
        x = np.arange(lo, hi)
        ax.fill_between(x, q[lo:hi], color=OBS_C, alpha=0.55, lw=0, label="observed")
        ax.plot(x, p[lo:hi], color=PRED_C, lw=1.0, label="predicted")
        ax.margins(x=0); ax.set_ylabel("P-site fraction", fontsize=8)
        ax.set_title(title, fontsize=9)
        if xlab:
            ax.set_xlabel("transcript position (nt)", fontsize=8)
        ax.tick_params(labelsize=7)

    fig, axes = plt.subplots(3, 2, figsize=(13, 9), gridspec_kw={"width_ratios": [1.7, 1]})

    # ---- Row 1: uORF (POLR1D) ----
    tx = args.uorf_tx; p, c, q = get(tx); gn = bio[tx][0]
    u5, cl, _ = (int(v) for v in cds[tx]); cds0, cds1 = u5, u5 + cl
    uo = next((w for w in calls[tx] if w[0] in ("uORF", "Overlap_uORF")), None)
    ua, ue = uo[1] - 1, uo[2]
    r_u5 = pearson(p[:u5], q[:u5])
    ax = axes[0, 0]
    hi = min(len(c), cds1 + 120)
    draw(ax, 0, hi, p, q, f"uORF  |  {gn} {tx}   5'UTR (uORF) profile r = {r_u5:.3f}", xlab=False)
    ax.axvspan(ua, ue, color="#ff9900", alpha=0.30, lw=0)
    ax.axvspan(cds0, hi, color="#cccccc", alpha=0.25, lw=0)
    ax.axvline(cds0, color="#333", ls="--", lw=0.9)
    top = ax.get_ylim()[1]
    ax.text((ua + ue) / 2, top * 0.9, "uORF", ha="center", fontsize=7.5, color="#b36b00")
    ax.text((cds0 + hi) / 2, top * 0.9, "CDS", ha="center", fontsize=7.5, color="#555")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    # Zoom on the uORF ITSELF, in the SAME frame-coloured style and at the same ~60 nt scale as the
    # CDS and lncRNA zooms below, so all three right-hand panels read identically. A uORF is a few
    # dozen nt inside a UTR that can be several hundred, so a full-UTR view renders it as a sliver,
    # and the fill+line style used for the context views cannot show reading frame at all -- which
    # is the whole point of a zoom. Frame is taken relative to the uORF's own start.
    uw = ue - ua
    if uw <= ZOOM_W:                      # short uORF: show all of it, centred, with a little flank
        upad = max(6, (ZOOM_W - uw) // 2)
        uz0, uz1 = max(0, ua - upad), min(len(c), ue + upad)
    else:                                 # long uORF: same rule as the CDS/lncRNA rows
        uz0, uz1 = best_zoom(c, ua, ue)
    ax = axes[0, 1]
    draw_periodicity(ax, uz0, uz1, p, q, ua,
                     f"uORF zoom [{uz0}, {uz1}) -- predicted P-sites colored by frame "
                     f"(uORF {uw} nt, {uw // 3} codons)", xlab=False)
    # Shade only the part of the uORF actually inside the view, and only when it does not cover the
    # whole panel -- a span wider than the axis renders as a solid block and reads as a background.
    if not (ua <= uz0 and ue >= uz1):
        ax.axvspan(max(ua, uz0), min(ue, uz1), color="#ff9900", alpha=0.18, lw=0, zorder=0)
    if uz0 <= ua < uz1:
        ax.axvline(ua, color="#b36b00", ls="--", lw=0.9)
    if uz0 < ue <= uz1:
        ax.axvline(ue, color="#b36b00", ls=":", lw=0.9)

    # ---- Row 2: CDS (RPL32) ----
    tx = args.cds_tx; p, c, q = get(tx); gn = bio[tx][0]
    u5, cl, _ = (int(v) for v in cds[tx]); cds0, cds1 = u5, u5 + cl
    r_whole = pearson(p, q)
    ax = axes[1, 0]
    draw(ax, 0, len(c), p, q, f"CDS  |  {gn} {tx}   whole-transcript profile r = {r_whole:.3f}", xlab=False)
    ax.axvspan(cds0, cds1, color="#cccccc", alpha=0.30, lw=0)
    ax.axvline(cds0, color="#333", ls="--", lw=0.9); ax.axvline(cds1, color="#333", ls=":", lw=0.9)
    top = ax.get_ylim()[1]
    ax.text((cds0 + cds1) / 2, top * 0.9, "CDS", ha="center", fontsize=7.5, color="#555")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    zlo, zhi = best_zoom(c, cds0, cds1)
    draw_periodicity(axes[1, 1], zlo, zhi, p, q, cds0,
                     f"CDS zoom [{zlo}, {zhi}) -- predicted P-sites colored by frame", xlab=False)

    # ---- Row 3: lncRNA (LINC02693) ----
    tx = args.lnc_tx; p, c, q = get(tx); gn = bio[tx][0]
    orf = calls[tx][0]; oa, oe = orf[1] - 1, orf[2]
    r_whole = pearson(p, q)
    ax = axes[2, 0]
    draw(ax, 0, len(c), p, q, f"lncRNA  |  {gn} {tx}   whole-transcript profile r = {r_whole:.3f}")
    ax.axvspan(oa, oe, color="#2ca02c", alpha=0.22, lw=0)
    ax.axvline(oa, color="#2ca02c", ls="--", lw=0.9); ax.axvline(oe, color="#2ca02c", ls=":", lw=0.9)
    top = ax.get_ylim()[1]
    ax.text((oa + oe) / 2, top * 0.9, "translated ORF", ha="center", fontsize=7.5, color="#1a7a1a")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    zlo, zhi = best_zoom(c, oa, oe)
    draw_periodicity(axes[2, 1], zlo, zhi, p, q, oa,
                     f"ORF zoom [{zlo}, {zhi}) -- predicted P-sites colored by frame")

    fig.suptitle(args.title or ("Observed vs predicted Ribo-seq P-site profile -- held-out "
                                "Hepatocytes (one-hot orf_v2_attn, tissue never seen in training)"),
                 fontsize=11, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    # every figure in this project ships raster + vector; slides need the vector copy
    pdf = out.with_suffix(".pdf")
    if pdf != out:
        fig.savefig(pdf, bbox_inches="tight")
        print(f"wrote {pdf}")


if __name__ == "__main__":
    main()
