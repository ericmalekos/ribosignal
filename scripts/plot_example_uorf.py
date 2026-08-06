#!/usr/bin/env python3
"""Observed vs predicted per-nt Ribo-seq profile for one example transcript, illustrating what a
good whole-transcript and 5'UTR (uORF) profile-Pearson score look like.

Three panels: (top) the whole transcript with the 5'UTR window shaded; (bottom-left) the 5'UTR
(uORF) zoom; (bottom-right) a CDS 3-nt-periodicity zoom. Default: the orf_v2_attn fold-0 model on
CKS1B (ENST00000368439.5), which scores whole-tx r=0.99 and uORF r=0.99. Runs on CPU in cas12a
(torch + numpy + matplotlib). The ORF-track env is set before importing the dataset loader.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")

ap = argparse.ArgumentParser()
ap.add_argument("--run", default=str(NEW / "results" / "improve" / "orf_v2_attn_f0"))
ap.add_argument("--tx", default="ENST00000368439.5")
ap.add_argument("--orf_track", default=str(NEW / "data" / "packed" / "orf_track_v2.npy"))
ap.add_argument("--out", default=None)
args = ap.parse_args()
# the dataset loader reads RIBO_ORF_TRACK on import/build, so set it before importing dataset
os.environ.setdefault("RIBO_ORF_TRACK", args.orf_track)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

sys.path.insert(0, str(NEW / "scripts"))
from dataset import BACKENDS, PackedStore, RiboDataset, load_split  # noqa: E402
from model import RiboSignalModel  # noqa: E402
from train import pearson  # noqa: E402


def gene_name(tx):
    with (NEW / "data" / "tx2biotype.tsv").open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ti, gi = h.index("tx_id"), h.index("gene_name")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ti] == tx:
                return f[gi]
    return "?"


def cds_bounds(tx):
    with (NEW / "data" / "tx2cds.tsv").open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len", "utr3_len")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["tx_id"]] == tx:
                return int(f[ix["utr5_len"]]), int(f[ix["cds_len"]]), int(f[ix["utr3_len"]])
    raise SystemExit(f"{tx} not in tx2cds.tsv")


def main():
    run = Path(args.run)
    cfg = json.loads((run / "args.json").read_text())
    store = PackedStore(tx_index=BACKENDS[cfg.get("emb_backend", "rinalmo")])
    _, _, test_ids = load_split(store, test_fold=cfg.get("test_fold", 0))
    if args.tx not in test_ids:
        raise SystemExit(f"{args.tx} not in the held-out test set of {run.name}")
    k = test_ids.index(args.tx)

    use_orf = cfg.get("use_orf_track", False)
    ds = RiboDataset(store, test_ids, input_mode=cfg.get("input_mode", "both"), use_orf=use_orf,
                     cov_norm=cfg.get("cov_norm", "raw"))
    d_emb = store.load_emb(test_ids[0]).shape[1]
    model = RiboSignalModel(
        d_emb=d_emb, channels=cfg["channels"], n_blocks=cfg["n_blocks"], dropout=0.0,
        extra_in=store.orf_channels if use_orf else 0,
        n_attn_layers=cfg.get("n_attn_layers", 0), n_heads=cfg.get("n_heads", 8),
        mixer=cfg.get("mixer", "transformer"), mamba_d_state=cfg.get("mamba_d_state", 16),
        mamba_d_conv=cfg.get("mamba_d_conv", 4), mamba_expand=cfg.get("mamba_expand", 2))
    model.load_state_dict(torch.load(run / "best.pt", map_location="cpu"))
    model.eval()

    it = ds[k]
    L = it["length"]
    with torch.no_grad():
        logits, _ = model(it["feats"].unsqueeze(0), torch.ones(1, L, dtype=torch.bool))
        p = torch.softmax(logits[0, :L], dim=0).numpy()
    c = it["counts"].numpy()[:L].astype(float)
    q = c / c.sum()

    u5, cl, _u3 = cds_bounds(args.tx)
    cds_start, cds_end = u5, u5 + cl
    gname = gene_name(args.tx)
    r_whole = pearson(p, q)
    r_uorf = pearson(p[:cds_start], c[:cds_start])
    u5_ps = int(c[:cds_start].sum())

    def draw(ax, lo, hi):
        x = np.arange(lo, hi)
        ax.fill_between(x, q[lo:hi], color="#1f77b4", alpha=0.55, lw=0, label="observed")
        ax.plot(x, p[lo:hi], color="#d62728", lw=0.9, label="predicted")
        ax.margins(x=0)
        ax.set_ylabel("P-site fraction")

    fig = plt.figure(figsize=(12, 6.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])
    axT, axU, axC = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    draw(axT, 0, L)
    axT.axvspan(0, cds_start, color="#ffd479", alpha=0.35, lw=0)
    axT.axvline(cds_start, color="#333", ls="--", lw=1)
    axT.axvline(cds_end, color="#333", ls=":", lw=1)
    top = axT.get_ylim()[1]
    axT.text(cds_start / 2, top * 0.88, "5'UTR (uORF)", ha="center", fontsize=8, color="#8a6d00")
    axT.text((cds_start + cds_end) / 2, top * 0.88, "CDS", ha="center", fontsize=8, color="#333")
    axT.text((cds_end + L) / 2, top * 0.88, "3'UTR", ha="center", fontsize=8, color="#555")
    maxfrac = c.max() / c.sum()
    axT.set_title(f"{gname}  {args.tx}   L={L} nt   whole r = {r_whole:.3f}   "
                  f"(max single nt = {maxfrac:.1%} of P-sites)", fontsize=10)
    axT.legend(loc="upper right", fontsize=8, framealpha=0.9)

    draw(axU, 0, cds_start)
    axU.set_title(f"5'UTR / uORF window  [0, {cds_start})   uORF r = {r_uorf:.3f}   "
                  f"({u5_ps} P-sites)", fontsize=10)
    axU.set_xlabel("transcript position (nt)")

    zlo, zhi = cds_start, min(cds_start + 66, cds_end, L)
    draw(axC, zlo, zhi)
    axC.set_title(f"CDS 3-nt periodicity zoom  [{zlo}, {zhi})", fontsize=10)
    axC.set_xlabel("transcript position (nt)")

    fig.suptitle(f"Observed vs predicted Ribo-seq P-site profile  ({run.name})",
                 fontsize=11, y=1.0)
    fig.tight_layout()
    out = Path(args.out) if args.out else (run / "figures" / f"example_{args.tx}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"{gname} {args.tx}: whole r={r_whole:.3f}  uORF r={r_uorf:.3f}  "
          f"u5_ps={u5_ps}  L={L}  cds=[{cds_start},{cds_end})")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
