#!/usr/bin/env python3
"""Final comparison table for the Hepatocytes-holdout architecture variants: the deployed 2-layer
transformer (standard) vs the expanded 4-layer transformer (attn4) vs the bidirectional Mamba mixer
(mamba). Same recipe otherwise (one-hot + mm1 coverage + no-Kozak track, noBrain 7 tissues).

Per model, on the Hepatocytes holdout:
  - profile accuracy: best val_pearson (+ epoch), model params
  - ORF-call quality at real depth (pred_obsdepth vs real): F1/precision/recall, CDS vs non-canon
  - standalone over-calling (pred_preddepth): total calls + novel count
Reads history.json + dropin/{real,pred_obsdepth,pred_preddepth}_collapsed.txt. Missing dropin ->
'pending'. Writes results/arch_compare_table.md.
"""
import json
import os

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
LOTO = f"{NEW}/results/loto"
MODELS = [
    ("onehot 2-attn (baseline)", "orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes"),
    ("onehot attn4", "orf_v2_attn4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes"),
    ("onehot mamba", "orf_v2_mamba_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes"),
    ("rinalmo (1280d)", "orf_v2_attn_rinalmo_noBrain_nokozak_mm1_holdout_Hepatocytes"),
    ("orthrus (512d)", "orf_v2_attn_orthrus_noBrain_nokozak_mm1_holdout_Hepatocytes"),
    ("hydrarna (1024d)", "orf_v2_attn_hydrarna_noBrain_nokozak_mm1_holdout_Hepatocytes"),
]


def load_calls(path):
    d = {}
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            oid = p[ci["ORF_ID"]]
            parts = oid.rsplit("_", 3)
            key = "_".join(parts[-3:]) if len(parts) >= 4 else oid
            d[key] = p[ci["ORF_type"]]
    return d


def cds(d):
    return {k for k, t in d.items() if t == "annotated"}


def nc(d):
    return {k for k, t in d.items() if t != "annotated"}


def prf(pred, ref):
    i = len(pred & ref)
    p = i / len(pred) if pred else 0.0
    r = i / len(ref) if ref else 0.0
    return p, r, (2 * p * r / (p + r) if (p + r) else 0.0)


def best_val(run):
    h = f"{run}/history.json"
    if not os.path.exists(h):
        return None, None, None
    hist = json.load(open(h))
    vps = [(e.get("val", {}).get("all", {}).get("pearson_median"), e["epoch"]) for e in hist]
    vps = [(v, e) for v, e in vps if v is not None]
    if not vps:
        return None, None, len(hist)
    bv, be = max(vps, key=lambda x: x[0])
    return bv, be, len(hist)


def n_params(run):
    try:
        import torch
        sd = torch.load(f"{run}/best.pt", map_location="cpu")
        return sum(v.numel() for v in sd.values() if hasattr(v, "numel"))
    except Exception:
        return None


def main():
    rows = []
    for label, sub in MODELS:
        run = f"{LOTO}/{sub}"
        bv, be, ne = best_val(run)
        real = load_calls(f"{run}/dropin/real_collapsed.txt")
        obs = load_calls(f"{run}/dropin/pred_obsdepth_collapsed.txt")
        pred = load_calls(f"{run}/dropin/pred_preddepth_collapsed.txt")
        row = {"label": label, "val": bv, "epoch": be, "nep": ne,
               "params": n_params(run) if bv is not None else None}
        if real and obs:
            row["all"] = prf(set(obs), set(real))
            row["cds"] = prf(cds(obs), cds(real))
            row["nc"] = prf(nc(obs), nc(real))
        if pred:
            row["pred_total"] = len(pred)
            row["pred_novel"] = sum(1 for t in pred.values() if t == "novel")
        rows.append(row)

    note = ("Same recipe; only the mixer differs. pred_obsdepth = shape at the real per-tx depth"
            " (isolates SHAPE quality); pred_preddepth = fully standalone.")
    cols = ("| model | params | val_pearson (ep) | obs F1 all | obs F1 CDS (P/R) |"
            " obs F1 non-canon (P/R) | standalone calls | novel |")
    out = ["# Arch-variant comparison -- Hepatocytes holdout (one-hot + mm1 + no-Kozak, noBrain)\n",
           note + "\n", cols, "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["val"] is None:
            out.append(f"| {r['label']} | training | pending | -- | -- | -- | -- | -- |")
            continue
        pm = f"{r['params']/1e6:.1f}M" if r.get("params") else "?"
        vp = f"{r['val']:.4f} (e{r['epoch']})"
        if "all" in r:
            fa = f"{r['all'][2]:.3f}"
            fc = f"{r['cds'][2]:.3f} ({r['cds'][0]:.2f}/{r['cds'][1]:.2f})"
            fn = f"{r['nc'][2]:.3f} ({r['nc'][0]:.2f}/{r['nc'][1]:.2f})"
        else:
            fa = fc = fn = "dropin pending"
        pt = f"{r.get('pred_total', '--')}"
        nv = f"{r.get('pred_novel', '--')}"
        out.append(f"| {r['label']} | {pm} | {vp} | {fa} | {fc} | {fn} | {pt} | {nv} |")
    out.append("\nHigher val_pearson = better per-nt profile. obs F1 = predicted SHAPE vs the real"
               " ORF calls at real depth. Lower standalone calls/novel at similar CDS = less"
               " over-calling. (Prior: attn4 + mamba trade recall for precision.)")
    txt = "\n".join(out)
    open(f"{NEW}/results/arch_compare_table.md", "w").write(txt + "\n")
    print(txt)
    print(f"\nwrote {NEW}/results/arch_compare_table.md")


if __name__ == "__main__":
    main()
