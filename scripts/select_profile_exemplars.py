#!/usr/bin/env python3
"""Rank observed-vs-predicted profile exemplars that are NOT carried by a single spike.

THE PROBLEM THIS SOLVES. A per-transcript Pearson r between observed and predicted P-site profiles can
be driven almost entirely by one enormous peak: if 60% of a transcript's observed signal sits on one
nucleotide and the model puts a peak there too, r is ~0.95 and the rest of the profile is irrelevant to
it. Such a transcript looks like a spectacular prediction and demonstrates almost nothing. The
exemplars in the current figure were picked on r alone and are vulnerable to exactly this.

THE CRITERION. Recompute r with the top observed positions DELETED. If agreement survives that, it is
distributed agreement -- the model got the shape, not one lucky peak.

  r_full        Pearson(obs, pred) over the whole transcript
  r_drop1       Pearson after removing the single highest-observed position
  r_drop3       ... the top 3
  r_drop_orf3   ... the top 3 WITHIN the called ORF (the region a figure actually zooms into)
  top1_frac     obs.max() / obs.sum()          -> want LOW
  spread90      fraction of nonzero positions needed to reach 90% of the signal -> want HIGH
  f0_obs        observed in-frame fraction inside the called ORF -> want HIGH (real periodicity)

A good exemplar has high r_drop3, low top1_frac, high spread90, real periodicity, and enough depth to
be believable. The composite `score` below is deliberately simple and its terms are printed alongside
it, so a human can re-rank on any single column instead of trusting the blend.

Scored per (dataset, transcript, called ORF). ORF class comes from that dump's own
`real_collapsed.txt`, so uORF / CDS / lncRNA / novel exemplars can be chosen per class, and the
observed calls are the same ones every other analysis uses.

Usage:
  select_profile_exemplars.py                       # all configured datasets -> one TSV each + combined
  select_profile_exemplars.py --dataset mouse_wang_liver --top 40
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
OUTDIR = NEW / "results/profile_exemplars"

# (dataset key, species, dump dir holding pred_profiles.npz + real_collapsed.txt), as a function of
# the architecture. The paths were mamba4 literals, so the shipped exemplar TSVs were all mamba4 with
# nothing on the files to say so -- a poster built as all-attn had no attn observed-vs-predicted
# source at all. `--model attn` selects the attn dumps, which already exist for every dataset here.
def datasets_for(model):
    return [
        ("human_hepatocytes", "human",
         f"results/loto/orf_v2_{model}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin"),
        ("human_ruizorera", "human", f"results/heldout/human_ruizorera/released_{model}/dropin"),
        ("mouse_wang_liver", "mouse", f"results/liver_released/{model}_mouse_wang_liver"),
        ("mouse_janich_liver", "mouse", f"results/liver_released/{model}_mouse_janich_liver_decon"),
        ("mouse_gse243134_liver", "mouse", f"results/liver_released/{model}_mouse_gse243134_liver"),
    ]


DATASETS = datasets_for("mamba4")
MIN_PSITES = 200      # below this the observed profile is too sparse to judge a shape
MIN_ORF_NT = 90       # matches the drop-in comparison's min_len


def pearson(a, b):
    if len(a) < 8:
        return float("nan")
    a = a.astype(np.float64); b = b.astype(np.float64)
    sa, sb = a.std(), b.std()
    if sa == 0 or sb == 0:
        return float("nan")
    return float(((a - a.mean()) * (b - b.mean())).mean() / (sa * sb))


def drop_top(obs, pred, k):
    """Pearson with the k highest-OBSERVED positions removed (chosen on obs, applied to both)."""
    if len(obs) <= k + 8:
        return float("nan")
    keep = np.ones(len(obs), bool)
    keep[np.argsort(obs)[-k:]] = False
    return pearson(obs[keep], pred[keep])


def load_calls(path):
    """{tx_id: (orf_type, tstart, tstop, gene_name, transcript_type)}, longest call per tx.

    transcript_type is carried so a lncRNA exemplar can be selected: RiboCode labels an ORF on a
    lncRNA as `novel`, which does not distinguish it from a novel ORF on a protein-coding transcript.
    """
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            try:
                a, b = int(f[c["ORF_tstart"]]), int(f[c["ORF_tstop"]])
            except (ValueError, KeyError):
                continue
            if b - a + 1 < MIN_ORF_NT:
                continue
            tx = f[c["transcript_id"]]
            prev = out.get(tx)
            if prev is None or (b - a) > (prev[2] - prev[1]):
                out[tx] = (f[c["ORF_type"]], a, b,
                           f[c["gene_name"]] if "gene_name" in c else "NA",
                           f[c["transcript_type"]] if "transcript_type" in c else "NA")
    return out


def score_dataset(key, species, d, top):
    dd = NEW / d
    npz = dd / "pred_profiles.npz"
    calls_p = dd / "real_collapsed.txt"
    if not npz.exists() or not calls_p.exists():
        print(f"  {key}: SKIP (missing npz or real_collapsed.txt)", file=sys.stderr)
        return []
    z = np.load(npz, allow_pickle=False)
    tx_ids = [str(t) for t in z["tx_ids"]]
    lengths = z["lengths"]; off = np.concatenate([[0], np.cumsum(lengths)])
    obs_flat, pred_flat = z["obs_flat"], z["pred_flat"]
    calls = load_calls(calls_p)

    rows = []
    for i, tx in enumerate(tx_ids):
        call = calls.get(tx)
        if call is None:
            continue
        a, b = off[i], off[i + 1]
        obs = obs_flat[a:b]; pred = pred_flat[a:b]
        tot = int(obs.sum())
        if tot < MIN_PSITES:
            continue
        otype, ts, te = call[0], call[1], call[2]
        ts = max(0, ts - 1); te = min(len(obs), te)          # collapsed coords are 1-based inclusive
        if te - ts < MIN_ORF_NT:
            continue
        o_orf, p_orf = obs[ts:te], pred[ts:te]
        if o_orf.sum() <= 0:
            continue
        nz = obs[obs > 0]
        srt = np.sort(nz)[::-1]
        cum = np.cumsum(srt) / srt.sum()
        spread90 = float((np.searchsorted(cum, 0.90) + 1) / len(nz))
        f0 = float(o_orf[0::3].sum() / o_orf.sum())
        rows.append({
            "dataset": key, "species": species, "tx_id": tx, "gene": call[3],
            "tx_biotype": call[4], "orf_type": otype, "orf_tstart": ts + 1, "orf_tstop": te, "tx_len": int(len(obs)),
            "psites_total": tot, "psites_in_orf": int(o_orf.sum()),
            "r_full": pearson(obs, pred),
            "r_drop1": drop_top(obs, pred, 1),
            "r_drop3": drop_top(obs, pred, 3),
            "r_orf": pearson(o_orf, p_orf),
            "r_drop_orf3": drop_top(o_orf, p_orf, 3),
            "top1_frac": float(obs.max() / tot),
            "spread90": spread90,
            "f0_obs_in_orf": f0,
        })

    for r in rows:
        # Composite: distributed agreement, penalised for spike dominance, requiring real periodicity.
        # Every term is in the TSV so this can be ignored and re-ranked by hand.
        rd = r["r_drop_orf3"]
        r["score"] = float("nan") if not np.isfinite(rd) else (
            rd * (1.0 - min(r["top1_frac"] / 0.35, 1.0)) * min(r["f0_obs_in_orf"] / 0.5, 1.0))
    rows = [r for r in rows if np.isfinite(r.get("score", float("nan")))]
    rows.sort(key=lambda r: -r["score"])
    print(f"  {key}: {len(rows):,} scoreable transcripts (>= {MIN_PSITES} P-sites, ORF >= {MIN_ORF_NT} nt)",
          file=sys.stderr)
    return rows[:top] if top else rows


# `model` is FIRST-CLASS, not optional: write_tsv emits exactly COLS, so stamping rows without
# listing it here would have silently dropped it (it did, on the first attn run).
COLS = ["model", "dataset", "species", "tx_id", "gene", "tx_biotype", "orf_type", "orf_tstart",
        "orf_tstop", "tx_len",
        "psites_total", "psites_in_orf", "r_full", "r_drop1", "r_drop3", "r_orf", "r_drop_orf3",
        "top1_frac", "spread90", "f0_obs_in_orf", "score"]


def write_tsv(rows, path):
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in COLS) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", help="only this dataset key")
    ap.add_argument("--top", type=int, default=200, help="rows kept per dataset (0 = all)")
    ap.add_argument("--outdir", default=str(OUTDIR))
    ap.add_argument("--model", default="mamba4", choices=["mamba4", "attn"],
                    help="which architecture's dumps to score. Default mamba4 (historical).")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    datasets = datasets_for(a.model)
    print(f"  model = {a.model}")

    combined = []
    for key, species, d in datasets:
        if a.dataset and key != a.dataset:
            continue
        if not (NEW / d / "pred_profiles.npz").exists():
            print(f"  SKIP {key}: no dump at {d}")
            continue
        rows = score_dataset(key, species, d, a.top)
        if not rows:
            continue
        # The model is stamped on every row, so an exemplar TSV can never again be architecture-
        # ambiguous the way the shipped mamba4 ones were.
        for r in rows:
            r["model"] = a.model
        write_tsv(rows, out / f"{key}_exemplars.tsv")
        combined.extend(rows)

    if combined:
        write_tsv(sorted(combined, key=lambda r: -r["score"]), out / "all_exemplars.tsv")
        print(f"\nwrote {out}/  ({len(combined):,} rows across "
              f"{len({r['dataset'] for r in combined})} datasets)")
        # A per-(dataset, class) shortlist is what a figure actually needs.
        print("\nTop candidate per dataset x ORF class "
              "(r_drop_orf3 = agreement with the top 3 observed peaks deleted):")
        hdr = (f"{'dataset':<22}{'class':<14}{'tx_id':<20}{'gene':<12}"
               f"{'r_full':>7}{'r_dropORF3':>11}{'top1':>7}{'f0':>6}{'psites':>8}")
        print(hdr); print("-" * len(hdr))
        seen = set()
        for r in sorted(combined, key=lambda r: -r["score"]):
            cls = ("lncRNA" if "lncRNA" in str(r["tx_biotype"]) else r["orf_type"])
            k = (r["dataset"], cls)
            if k in seen:
                continue
            seen.add(k)
            print(f"{r['dataset']:<22}{cls:<14}{r['tx_id']:<20}{str(r['gene'])[:11]:<12}"
                  f"{r['r_full']:>7.3f}{r['r_drop_orf3']:>11.3f}{r['top1_frac']:>7.3f}"
                  f"{r['f0_obs_in_orf']:>6.2f}{r['psites_total']:>8,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
