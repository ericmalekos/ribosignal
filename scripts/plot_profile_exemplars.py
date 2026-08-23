#!/usr/bin/env python3
"""Plot observed vs predicted per-nt P-site profiles for a chosen set of transcripts.

Driven entirely by a transcript list, so the figure is re-specifiable without touching code. Pick rows
from `results/profile_exemplars/*_exemplars.tsv` (produced by `select_profile_exemplars.py`, which
ranks for agreement that SURVIVES deleting the top observed peaks) and pass them here.

Each transcript gets two panels:
  left   the whole transcript, observed (grey fill) vs predicted (line), the called ORF shaded
  right  a zoom on a window inside the ORF, observed as outline bars and predicted coloured BY FRAME,
         so the 3-nt periodicity is visible rather than asserted

Every panel's data is also written to a TSV next to the figure (one row per nucleotide plotted), so
the figure can be restyled in any other tool without re-running the model.

Selection inputs:
  --tsv FILE       exemplar TSV; rows are used in file order unless --tx is given
  --tx  ID [ID..]  explicit transcript IDs (must appear in --tsv, which supplies their coordinates)
  --dataset KEY    restrict --tsv rows to one dataset
  --per-class N    take the best N per (dataset, class) instead of the raw top rows

Usage:
  plot_profile_exemplars.py --tsv results/profile_exemplars/human_hepatocytes_exemplars.tsv --per-class 1
  plot_profile_exemplars.py --tsv results/profile_exemplars/all_exemplars.tsv \\
      --tx ENST00000908311.1 ENSMUST00000022416.15 --out figures/profile_exemplars_mixed
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
# Must match select_profile_exemplars.DATASETS.
def dumps_for(model):
    """Dump dir per dataset, as a function of architecture (see select_profile_exemplars)."""
    return {
        "human_hepatocytes":
            f"results/loto/orf_v2_{model}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin",
        "human_ruizorera": f"results/heldout/human_ruizorera/released_{model}/dropin",
        "mouse_wang_liver": f"results/liver_released/{model}_mouse_wang_liver",
        "mouse_janich_liver": f"results/liver_released/{model}_mouse_janich_liver_decon",
        "mouse_gse243134_liver": f"results/liver_released/{model}_mouse_gse243134_liver",
    }


DUMPS = dumps_for("mamba4")
# One-hot universe FASTA per dataset -- EXACTLY the file each dump was produced with, so the sequence
# under a profile is the sequence the model actually saw. Taking the "obvious" species-level FASTA
# instead would silently mis-pair transcripts whose universe differs between arms.
FASTA = {
    "human_hepatocytes": "data/union_universe.fa",
    "human_ruizorera": "data/union_universe.fa",
    "mouse_wang_liver": "data/heldout_refs/mouse_wang_expressed_universe.fa",
    "mouse_janich_liver": "data/heldout_refs/mouse_janich_decon_universe.fa",
    "mouse_gse243134_liver": "data/heldout_refs/mouse_gse243134_universe.fa",
}
CODON2AA = {}
for _i, _b1 in enumerate("TCAG"):
    for _b2 in "TCAG":
        for _b3 in "TCAG":
            CODON2AA[_b1 + _b2 + _b3] = (
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
                [("TCAG".index(_b1) * 16) + ("TCAG".index(_b2) * 4) + "TCAG".index(_b3)])

# Observed = dark blue, predicted = dark red, mirrored against each other in the whole-transcript
# panel. OBS_C is used ONLY there and in the legend: the ORF zoom draws observed as unfilled outline
# bars so the frame-coloured predicted bars stay readable behind them.
OBS_C, PRED_C = "#1F4E79", "#C0392B"
FRAME_C = ["#1B7F5A", "#E67E22", "#7D8CA3"]     # frame 0 (in-frame), 1, 2
_CACHE: dict[str, dict] = {}
_SEQ: dict[str, dict] = {}


def load_seq(key, want):
    """Sequences for `want` transcripts from that dataset's universe FASTA.

    Streamed and filtered to the requested IDs: union_universe.fa is 208 MB and holding all of it is
    pointless when a figure needs a handful of transcripts.
    """
    if key not in FASTA:
        return {}
    cache = _SEQ.setdefault(key, {})
    need = {t for t in want if t not in cache}
    if not need:
        return cache
    path = NEW / FASTA[key]
    if not path.exists():
        print(f"  warning: {path} missing; no sequence track for {key}", file=sys.stderr)
        return cache
    cur, buf = None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in need:
                    cache[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            elif cur in need:
                buf.append(line.strip())
    if cur in need:
        cache[cur] = "".join(buf).upper()
    return cache


def codon_at(seq, orf_start, pos):
    """(codon, aa, codon_index) for the in-frame codon beginning at `pos`, else (None, None, None)."""
    if seq is None or pos + 3 > len(seq) or (pos - orf_start) % 3 != 0:
        return None, None, None
    cod = seq[pos:pos + 3]
    return cod, CODON2AA.get(cod, "X"), (pos - orf_start) // 3


def load_dump(key):
    if key in _CACHE:
        return _CACHE[key]
    if key not in DUMPS:
        raise SystemExit(f"unknown dataset {key!r}; known: {', '.join(sorted(DUMPS))}")
    z = np.load(NEW / DUMPS[key] / "pred_profiles.npz", allow_pickle=False)
    ids = [str(t) for t in z["tx_ids"]]
    off = np.concatenate([[0], np.cumsum(z["lengths"])])
    _CACHE[key] = {"idx": {t: i for i, t in enumerate(ids)}, "off": off,
                   "obs": z["obs_flat"], "pred": z["pred_flat"]}
    return _CACHE[key]


def get_profiles(key, tx):
    d = load_dump(key)
    i = d["idx"].get(tx)
    if i is None:
        raise SystemExit(f"{tx} not present in dump for {key}")
    a, b = d["off"][i], d["off"][i + 1]
    return d["obs"][a:b].astype(float), d["pred"][a:b].astype(float)


def read_tsv(path, dataset=None):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if dataset:
        rows = [r for r in rows if r["dataset"] == dataset]
    for r in rows:
        for k in ("orf_tstart", "orf_tstop", "tx_len", "psites_total"):
            r[k] = int(r[k])
        for k in ("r_full", "r_drop3", "r_drop_orf3", "top1_frac", "f0_obs_in_orf", "score"):
            r[k] = float(r[k])
        r["cls"] = "lncRNA" if "lncRNA" in r.get("tx_biotype", "") else r["orf_type"]
    return rows


def pick_zoom(obs, ts, te, width, mode="start"):
    """Which `width` nt of the ORF to show. See --zoom-mode.

    `dense` was the original and only behaviour, and it is NOT neutral: because agreement tracks read
    depth, the densest window correlates better than the ORF as a whole in 9 of 12 current exemplars
    (by up to +0.26 r) and carries 2.4-3.1x the signal of an average window on long ORFs. That is a
    best-case view. `start` is the default instead: it involves no selection on the data at all, it
    always shows the initiator codon, and it is the same window definition for every transcript, so
    panels are comparable. `dense` remains available for slides where legibility matters more.
    """
    ts, te = max(0, ts), min(len(obs), te)
    if te - ts <= width:
        return ts, te
    k = min(width, te - ts)
    if mode == "start":
        return ts, ts + k
    if mode == "center":
        s = ts + (te - ts - k) // 2
        return s, s + k
    csum = np.concatenate([[0.0], np.cumsum(obs[ts:te])])
    dens = csum[k:] - csum[:-k]                 # sliding-window sums of OBSERVED counts
    if mode == "median":
        # the window whose density is closest to the median density -- typical, not best-case
        s = ts + int(np.argmin(np.abs(dens - np.median(dens))))
    else:                                        # dense
        s = ts + int(np.argmax(dens))
    return s, s + k


def plot_one(ax_full, ax_zoom, row, zoom_width, rows_out, style="mirror",
             seq=None, seq_track="codon", zoom_xticks=False, zoom_mode="start"):
    obs, pred = get_profiles(row["dataset"], row["tx_id"])
    ts, te = row["orf_tstart"] - 1, row["orf_tstop"]
    # Both tracks are normalised to their own total so shape is comparable; the raw counts and the
    # scale factors go into the TSV so nothing is lost.
    o_n = obs / obs.sum() if obs.sum() else obs
    p_n = pred / pred.sum() if pred.sum() else pred
    x = np.arange(len(obs))

    if style == "mirror":
        # Observed up, predicted down. In overlay the two tracks agree so closely on a good exemplar
        # that the predicted line completely hides the observed fill, and the panel reads as though
        # only one track were plotted. Mirroring keeps both visible and makes disagreement the thing
        # that stands out.
        ax_full.fill_between(x, o_n, step="mid", color=OBS_C, linewidth=0, label="observed")
        ax_full.fill_between(x, -p_n, step="mid", color=PRED_C, linewidth=0, alpha=0.85,
                             label="predicted")
        ax_full.axhline(0, color="#444", lw=0.5)
        lim = max(o_n.max(), p_n.max()) * 1.08 or 1.0
        ax_full.set_ylim(-lim, lim)
        ax_full.set_ylabel("P-site fraction\n(obs up / pred down)", fontsize=7)
    else:
        ax_full.fill_between(x, o_n, step="mid", color=OBS_C, linewidth=0, label="observed")
        ax_full.plot(x, p_n, color=PRED_C, lw=0.65, label="predicted")
        ax_full.set_ylabel("P-site fraction", fontsize=7)
    ax_full.axvspan(ts, te, color="#000000", alpha=0.055, lw=0)
    ax_full.set_xlim(0, len(obs))
    ax_full.tick_params(labelsize=6)
    ax_full.spines[["top", "right"]].set_visible(False)
    ax_full.set_title(
        f"{row['cls']}  |  {row['gene']}  {row['tx_id']}   "
        f"r={row['r_full']:.3f}, r(top-3 removed)={row['r_drop_orf3']:.3f}, "
        f"top1={row['top1_frac']:.3f}", fontsize=7.2, loc="left")

    zs, ze = pick_zoom(obs, ts, te, zoom_width, zoom_mode)
    xs = np.arange(zs, ze)
    frames = (xs - ts) % 3
    for f in (0, 1, 2):
        m = frames == f
        if m.any():
            ax_zoom.bar(xs[m], p_n[zs:ze][m], width=0.9, color=FRAME_C[f], linewidth=0,
                        label=f"frame {f}" + (" (in-frame)" if f == 0 else ""))
    ax_zoom.bar(xs, o_n[zs:ze], width=0.9, facecolor="none", edgecolor="#222", linewidth=0.45,
                label="observed")
    ax_zoom.set_xlim(zs - 0.5, ze - 0.5)
    ax_zoom.tick_params(labelsize=6)
    ax_zoom.spines[["top", "right"]].set_visible(False)
    # The codon track under this axis IS the x axis; numeric positions on top of it are redundant
    # clutter, and the panel title already states the window. Restore them with --zoom-xticks.
    if seq is not None and seq_track != "none" and not zoom_xticks:
        ax_zoom.set_xticks([])
    MODE_LAB = {"start": "from ORF start", "dense": "DENSEST window (best case)",
                "median": "median-density window", "center": "ORF centre"}
    ax_zoom.set_title(f"ORF zoom ({zs}-{ze}, {MODE_LAB[zoom_mode]}) -- predicted coloured by frame",
                      fontsize=7.2, loc="left")

    # ---- sanity check: a called ORF must translate without internal stops ----
    # Free and decisive. If the ORF coordinates, the reading frame, or the transcript-to-FASTA pairing
    # were wrong, the translation would be peppered with stops. Getting M...* with none in between
    # confirms all three at once, so this runs on every exemplar rather than being trusted.
    if seq is not None:
        prot = [codon_at(seq, ts, q)[1] for q in range(ts, te) if (q - ts) % 3 == 0]
        prot = [a for a in prot if a]
        internal = sum(1 for a in prot[:-1] if a == "*")
        if internal:
            print(f"  WARNING {row['tx_id']} ({row['dataset']}): called ORF translates with "
                  f"{internal} INTERNAL stop codon(s) -- coordinates, frame or FASTA pairing is "
                  f"wrong for this transcript", file=sys.stderr)
        row["_prot_len"] = len(prot)
        row["_prot_ok"] = (not internal) and bool(prot) and prot[0] == "M" and prot[-1] == "*"

    # ---- sequence track under the zoom ----
    if seq is not None and seq_track != "none":
        cods = [(q, *codon_at(seq, ts, q)) for q in range(zs, ze)
                if (q - ts) % 3 == 0 and q + 3 <= len(seq)]
        cods = [c for c in cods if c[1] is not None]
        # Codon triplets are 3 characters wide against a 3-nt slot, so they crowd about three times
        # sooner than single-letter amino acids do. Thin on the width of what is actually drawn.
        wide = seq_track in ("codon", "both")
        cap = 22 if wide else 26
        every = 1 if len(cods) <= cap else (2 if len(cods) <= cap * 2 else 5)
        if seq_track == "both":
            y_top, y_bot = -0.19, -0.30            # aa over codon
        else:
            y_top, y_bot = -0.19, -0.19
        for q, cod, aa, ci in cods:
            if ci % every:
                continue
            xc = q + 1                              # centre of the codon
            stop, start = aa == "*", ci == 0
            col = "#B03A2E" if stop else ("#1B7F5A" if start else "#333")
            if seq_track in ("aa", "both"):
                ax_zoom.text(xc, y_top, aa, transform=ax_zoom.get_xaxis_transform(),
                             ha="center", va="top", fontsize=6.6, family="monospace",
                             fontweight="bold" if (stop or start) else "normal", color=col)
            if seq_track in ("codon", "both"):
                ax_zoom.text(xc, y_bot, cod, transform=ax_zoom.get_xaxis_transform(),
                             ha="center", va="top", fontsize=5.6, family="monospace",
                             fontweight="bold" if (stop or start) else "normal",
                             color=col if seq_track == "codon" else "#777")
        lab = {"aa": "amino acid", "codon": "codon", "both": "amino acid / codon"}[seq_track]
        note = f"{lab}, frame of the called ORF" + ("" if every == 1 else f" (every {every})")
        # transAxes, NOT get_xaxis_transform(): the latter takes x in DATA coordinates, so x=0.0
        # lands at transcript position 0 -- thousands of nt left of a zoom window -- and matplotlib
        # expands the whole figure to contain it.
        ax_zoom.text(0.0, y_bot - 0.10, note, transform=ax_zoom.transAxes, ha="left", va="top",
                     fontsize=5.6, color="#888", style="italic")

    for i in range(len(obs)):
        cod_i = codon_at(seq, ts, i) if (seq is not None and ts <= i < te) else None
        rows_out.append({
            "dataset": row["dataset"], "tx_id": row["tx_id"], "gene": row["gene"],
            "class": row["cls"], "zoom_mode": zoom_mode, "position_nt": i,
            "in_called_orf": int(ts <= i < te), "in_zoom_window": int(zs <= i < ze),
            "frame_rel_orf": int((i - ts) % 3) if ts <= i < te else -1,
            "obs_psites": int(obs[i]), "pred_raw": float(pred[i]),
            "obs_frac": float(o_n[i]), "pred_frac": float(p_n[i]),
            # `nt` is given for every position, but codon/aa ONLY inside the called ORF. Extending
            # the ORF's reading frame across the whole transcript yields codons at negative indices
            # that translate to plausible-looking residues while meaning nothing -- a reader
            # filtering on `aa != ""` would silently pick up 5'UTR "protein".
            "nt": (seq[i] if seq is not None and i < len(seq) else ""),
            "codon": (cod_i[0] or "") if cod_i else "",
            "aa": (cod_i[1] or "") if cod_i else "",
            "codon_idx_in_orf": (cod_i[2] if cod_i and cod_i[2] is not None else ""),
        })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--tx", nargs="*", default=None)
    ap.add_argument("--dataset")
    ap.add_argument("--per-class", type=int, default=0)
    ap.add_argument("--zoom-width", type=int, default=60)
    ap.add_argument("--zoom-mode", choices=["start", "dense", "median", "center"], default="start",
                    help="which window of the ORF the right panel shows. start (default) = the first "
                         "--zoom-width nt from the initiator: no selection on the data, always shows "
                         "the start codon, identical definition for every transcript. dense = the "
                         "highest-observed-signal window, which is a BEST-CASE view (2.4-3.1x the "
                         "signal of an average window on long ORFs, and higher r than the whole ORF "
                         "in 9 of 12 exemplars). median = a typical-density window. center = middle.")
    ap.add_argument("--seq-track", choices=["codon", "aa", "both", "none"], default="codon",
                    help="sequence shown under the ORF zoom: codon triplets (default), "
                         "single-letter amino acids, both, or none")
    ap.add_argument("--zoom-xticks", action="store_true",
                    help="keep numeric position ticks on the ORF zoom axis. Off by default: the "
                         "codon track replaces them and the panel title already states the range. "
                         "The whole-transcript panel always keeps its position axis.")
    ap.add_argument("--style", choices=["mirror", "overlay"], default="mirror",
                    help="mirror = observed up / predicted down (default, both visible); "
                         "overlay = observed fill with the predicted line on top")
    ap.add_argument("--out", default=str(NEW / "figures/profile_exemplars/profile_exemplars"))
    ap.add_argument("--model", default="mamba4", choices=["mamba4", "attn"],
                    help="architecture whose dump supplies the predicted profile. MUST match the "
                         "model the --tsv exemplars were selected from.")
    a = ap.parse_args()
    globals()["DUMPS"] = dumps_for(a.model)
    # A TSV selected from one architecture plotted against another's dump would silently pair the
    # wrong predicted profile with the right observed one, so the mismatch is refused rather than
    # drawn. Older TSVs have no `model` column and are let through with a warning.
    _mods = {r.get("model") for r in read_tsv(a.tsv, a.dataset) if r.get("model")}
    if _mods and _mods != {a.model}:
        raise SystemExit(f"--model {a.model} but {a.tsv} was selected from {sorted(_mods)}. "
                         f"Re-run select_profile_exemplars.py --model {a.model}, or pass "
                         f"--model {sorted(_mods)[0]}.")
    if not _mods:
        print(f"  WARNING: {a.tsv} has no `model` column (pre-2026-08-15). Assuming {a.model}.")

    rows = read_tsv(a.tsv, a.dataset)
    if not rows:
        raise SystemExit(f"no rows in {a.tsv}" + (f" for dataset {a.dataset}" if a.dataset else ""))
    if a.tx:
        by = {r["tx_id"]: r for r in rows}
        missing = [t for t in a.tx if t not in by]
        if missing:
            raise SystemExit(f"not in {a.tsv}: {', '.join(missing)}")
        sel = [by[t] for t in a.tx]
    elif a.per_class:
        seen, sel = {}, []
        for r in rows:                       # TSV is already score-sorted
            k = (r["dataset"], r["cls"])
            if seen.get(k, 0) >= a.per_class:
                continue
            seen[k] = seen.get(k, 0) + 1
            sel.append(r)
    else:
        sel = rows[:3]

    n = len(sel)
    row_h = 1.75 + (0.50 if a.seq_track == "both" else 0.32 if a.seq_track != "none" else 0)
    fig, axes = plt.subplots(n, 2, figsize=(11.5, row_h * n), squeeze=False,
                             gridspec_kw=dict(width_ratios=[2.1, 1.0]))
    rows_out = []
    # One FASTA pass per dataset for all of its selected transcripts.
    seqs = {}
    for key in {r["dataset"] for r in sel}:
        seqs[key] = load_seq(key, [r["tx_id"] for r in sel if r["dataset"] == key])
    for i, r in enumerate(sel):
        seq = seqs.get(r["dataset"], {}).get(r["tx_id"])
        if seq is None and a.seq_track != "none":
            print(f"  warning: no sequence for {r['tx_id']} ({r['dataset']}); "
                  "plotting without a sequence track", file=sys.stderr)
        plot_one(axes[i][0], axes[i][1], r, a.zoom_width, rows_out, a.style, seq, a.seq_track,
                 a.zoom_xticks, a.zoom_mode)
    axes[-1][0].set_xlabel("transcript position (nt)", fontsize=7)
    # Only label the zoom axis when it still shows numeric positions; with the codon track in their
    # place a "transcript position" label under the codons would be describing the wrong thing.
    if a.seq_track == "none" or a.zoom_xticks:
        pad = {"both": 34, "aa": 22, "codon": 22, "none": 4}[a.seq_track]
        axes[-1][1].set_xlabel("transcript position (nt)", fontsize=7, labelpad=pad)
    pred_handle = (Patch(facecolor=PRED_C, label="predicted (plotted downward)")
                   if a.style == "mirror" else Line2D([], [], color=PRED_C, lw=1.2, label="predicted"))
    fig.legend(handles=[Patch(facecolor=OBS_C, label="observed P-sites"),
                        pred_handle,
                        Patch(facecolor=FRAME_C[0], label="predicted, in-frame"),
                        Patch(facecolor=FRAME_C[1], label="predicted, frame 1"),
                        Patch(facecolor=FRAME_C[2], label="predicted, frame 2")],
               loc="upper center", ncol=5, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, 1.005))
    fig.tight_layout(rect=(0, 0, 1, 0.975))

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    tsv = f"{out}.tsv"
    cols = ["dataset", "tx_id", "gene", "class", "zoom_mode", "position_nt", "in_called_orf", "in_zoom_window",
            "frame_rel_orf", "nt", "codon", "aa", "codon_idx_in_orf",
            "obs_psites", "pred_raw", "obs_frac", "pred_frac"]
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows_out:
            fh.write("\t".join(f"{r[c]:.8g}" if isinstance(r[c], float) else str(r[c])
                               for c in cols) + "\n")
    print(f"wrote {out}.pdf/.png and {tsv} ({len(rows_out):,} nt rows)")
    for r in sel:
        print(f"  {r['dataset']:<22}{r['cls']:<14}{r['tx_id']:<21}{str(r['gene'])[:12]:<13}"
              f"r={r['r_full']:.3f}  r_drop_orf3={r['r_drop_orf3']:.3f}  "
              f"top1={r['top1_frac']:.3f}  {r['psites_total']:,} P-sites")
    return 0


if __name__ == "__main__":
    sys.exit(main())
