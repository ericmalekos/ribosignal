#!/usr/bin/env python3
"""Figure C10: input saliency -- what SEQUENCE the model reads to place periodic P-sites.

For each exemplar transcript (the uORF / CDS / lncRNA ORF from prediction_examples, held-out Hepatocytes),
we set the target = the predicted profile log-probability summed over the ORF's in-frame (frame-0)
positions, and backprop to the one-hot sequence input. |d target / d one-hot| per nucleotide = how much
each position's identity drives the model's decision to put periodic ribosome signal on this ORF.

Expectation (and the interpretability payoff): saliency concentrates at the START codon + its Kozak
context, and shows 3-nt-periodic structure through the ORF body -- i.e. the model has learned to read the
start context and reading frame, not just copy RNA-seq coverage.

Panel per exemplar: (left) per-nt sequence saliency along the transcript with the ORF region shaded;
(right) zoom on +/- 30 nt around the start codon, saliency bars colored by frame, start codon marked.

CPU-only (~5M params). Reuses the project's PackedStore / RiboDataset / RiboSignalModel. cas12a env
(torch 2.x CPU).

CHECKPOINT (repointed 2026-08-08) = the released union attn model. This previously read
`orf_v2_attn_onehot_holdout_Hepatocytes` with `orf_track_v2.npy`. BOTH had to change, together: the
released models trained on the no-Kozak track (feedback_kozak_never_set), so pairing a released
checkpoint with the Kozak-heuristic track would feed the model an input channel it never saw and the
saliency map would be an artifact of that mismatch rather than of the model.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/C10_saliency"
RUN = NEW / "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
# THREE env vars, set together, exactly as train_loto_union.sbatch set them. They all name the same
# universe, so any one of them left at its default silently swaps the universe out from under the
# other two. dataset.py defaults are the FIBROBLAST universe (`fibroblast_universe.fa`, unsuffixed
# packs), and data/packed/orf_track_v2_nokozak.npy is ALSO kozak=none -- so a Fibroblast-track run
# is internally consistent, produces a plausible saliency map, and passes a kozak check. `setdefault`
# is deliberate: an explicit env var from the caller still wins.
os.environ.setdefault("RIBO_PACK_SUFFIX", "union")
os.environ.setdefault("RIBO_ORF_TRACK", str(NEW / "data/packed_union/orf_track_v2.npy"))
os.environ.setdefault("RIBO_ONEHOT_FASTA", str(NEW / "data/union_universe.fa"))
sys.path.insert(0, str(NEW / "scripts"))
from dataset import BACKENDS, PackedStore, RiboDataset, tissue_pack_dir  # noqa: E402
from model import RiboSignalModel  # noqa: E402

# exemplar (tx, gene, ORF tstart, tstop [python half-open], label) -- from Hepatocytes RiboCode calls
EXEMPLARS = [
    ("ENST00000399697.7", "POLR1D", 147, 183, "uORF"),        # 5'UTR uORF (1-based 148..183)
    ("ENST00000429711.7", "RPL32", 77, 485, "CDS"),           # annotated CDS
    ("ENST00000647872.1", "LINC02693", 258, 924, "lncRNA ORF"),  # novel lncRNA ORF
]
FRAME_C = {0: "#2C6FBB", 1: "#E08A2B", 2: "#7B5EA7"}
SEQ_C = "#2C6FBB"


def load_cfg():
    aj = RUN / "args.json"
    cfg = json.load(open(aj)) if aj.exists() else {}
    # defaults = the deployed onehot orf_v2_attn config
    return dict(channels=cfg.get("channels", 256), n_blocks=cfg.get("n_blocks", 10),
                n_attn_layers=cfg.get("n_attn_layers", 2), n_heads=cfg.get("n_heads", 8),
                input_mode=cfg.get("input_mode", "both"), cov_norm=cfg.get("cov_norm", "global_mean"),
                use_orf=cfg.get("use_orf_track", True),
                learn_start_context=cfg.get("learn_start_context", False))


def main():
    cfg = load_cfg()
    store = PackedStore(pack=tissue_pack_dir("Hepatocytes"), tx_index=BACKENDS["onehot"])
    txs = [e[0] for e in EXEMPLARS if e[0] in store.row]
    ds = RiboDataset(store, txs, input_mode=cfg["input_mode"], use_orf=cfg["use_orf"],
                     cov_norm=cfg["cov_norm"])
    d_emb = store.load_emb(txs[0]).shape[1]
    extra_in = store.orf_channels if cfg["use_orf"] else 0
    model = RiboSignalModel(d_emb=d_emb, channels=cfg["channels"], n_blocks=cfg["n_blocks"],
                            dropout=0.0, extra_in=extra_in, n_attn_layers=cfg["n_attn_layers"],
                            n_heads=cfg["n_heads"], learn_start_context=cfg["learn_start_context"])
    model.load_state_dict(torch.load(RUN / "best.pt", map_location="cpu"))
    model.eval()

    rows = {e[0]: e for e in EXEMPLARS}
    results = {}
    for k, tx in enumerate(txs):
        item = ds[k]
        feats = item["feats"].unsqueeze(0).clone().requires_grad_(True)   # (1, L, C)
        mask = torch.ones(1, feats.shape[1], dtype=torch.bool)
        logits, _ = model(feats, mask)
        logp = torch.log_softmax(logits[0], dim=0)                        # (L,)
        _, _, ts, te, _ = rows[tx]
        fr0 = torch.arange(ts, te, 3)                                     # in-frame positions
        target = logp[fr0].sum()
        model.zero_grad(set_to_none=True)
        if feats.grad is not None:
            feats.grad = None
        target.backward()
        g = feats.grad[0].detach().numpy()                               # (L, C)
        sal_seq = np.abs(g[:, :4]).sum(axis=1)                           # per-nt one-hot saliency
        sal_cov = np.abs(g[:, -1])                                       # coverage-channel saliency
        results[tx] = dict(sal_seq=sal_seq, sal_cov=sal_cov, ts=ts, te=te)
        print(f"  {tx} ({rows[tx][1]} {rows[tx][4]}): L={len(sal_seq)} "
              f"start-codon saliency rank { (sal_seq >= sal_seq[ts]).sum() }/{len(sal_seq)} "
              f"(1 = highest); mean-in-ORF/mean-out = "
              f"{sal_seq[ts:te].mean() / max(sal_seq.mean(), 1e-9):.2f}x")

    # ---- plot ----
    fig, axes = plt.subplots(len(txs), 2, figsize=(9.6, 2.5 * len(txs)),
                             gridspec_kw=dict(width_ratios=[2.1, 1]))
    if len(txs) == 1:
        axes = axes[None, :]
    for i, tx in enumerate(txs):
        r = results[tx]
        gene, lab = rows[tx][1], rows[tx][4]
        sal = r["sal_seq"]
        ts, te = r["ts"], r["te"]
        L = len(sal)
        axL, axR = axes[i]

        # left: full-transcript saliency, ORF shaded
        axL.fill_between(np.arange(L), sal, color=SEQ_C, lw=0, alpha=0.85)
        axL.axvspan(ts, te, color="#f2c94c", alpha=0.25, lw=0)
        axL.axvline(ts, color="#c0392b", lw=1.0, ls="--")
        axL.set_xlim(0, L)
        axL.set_ylabel("seq saliency", fontsize=8)
        axL.set_title(f"{gene}  ({lab}, {tx})", fontsize=8.5, loc="left")
        axL.spines[["top", "right"]].set_visible(False)
        if i == len(txs) - 1:
            axL.set_xlabel("transcript position (nt)", fontsize=8.5)

        # right: start-codon zoom, saliency bars colored by frame
        lo, hi = max(0, ts - 30), min(L, ts + 30)
        pos = np.arange(lo, hi)
        cols = [FRAME_C[(p - ts) % 3] for p in pos]
        axR.bar(pos - ts, sal[lo:hi], color=cols, width=0.9, lw=0)
        axR.axvline(0, color="#c0392b", lw=1.0, ls="--")
        axR.text(0, axR.get_ylim()[1] * 0.98, "start", color="#c0392b", fontsize=6.8,
                 ha="center", va="top")
        axR.set_xlim(-30, 30)
        axR.set_xlabel("nt relative to start", fontsize=8) if i == len(txs) - 1 else None
        axR.set_title("start-codon zoom (frame-colored)", fontsize=7.8, loc="left")
        axR.spines[["top", "right"]].set_visible(False)

    from matplotlib.patches import Patch
    handles = [Patch(color=FRAME_C[f], label=f"frame {f}") for f in (0, 1, 2)]
    fig.legend(handles=handles, frameon=False, fontsize=7.5, ncol=3, loc="upper right",
               bbox_to_anchor=(0.99, 1.005))
    fig.suptitle("What sequence the model reads: input saliency concentrates at the start context "
                 "(held-out Hepatocytes)", fontsize=10, y=1.01)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"C10_saliency.{ext}", dpi=300, bbox_inches="tight")
    (HERE / "C10_values.json").write_text(json.dumps({
        tx: {"gene": rows[tx][1], "label": rows[tx][4], "ts": int(results[tx]["ts"]),
             "te": int(results[tx]["te"]),
             "start_saliency_rank": int((results[tx]["sal_seq"] >= results[tx]["sal_seq"][results[tx]["ts"]]).sum()),
             "L": int(len(results[tx]["sal_seq"])),
             "mean_in_orf_over_mean": float(results[tx]["sal_seq"][results[tx]["ts"]:results[tx]["te"]].mean()
                                            / max(results[tx]["sal_seq"].mean(), 1e-9))}
        for tx in txs}, indent=2) + "\n")
    print(f"wrote C10_saliency.pdf/.png for {len(txs)} exemplars")


if __name__ == "__main__":
    main()
