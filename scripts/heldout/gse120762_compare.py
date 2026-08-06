#!/usr/bin/env python3
"""Compare the deployed model's GSE120762 predictions against observed Ribo-seq, per condition
(NT vs LPS), and assess the NT->LPS differential translation the model recovers.

Inputs (self-contained -- the dump already paired predicted vs observed on a shared per-nt axis):
  proteogenomics/data/gse120762/pred_{nt,lps}/pred_profiles.npz
    tx_ids (U), lengths (int), pred_flat (f32), obs_flat (i32 observed P-sites), pred_total (f64)
  data/packed_heldout_mouse_gse120762_{nt,lps}/{coverage.npy,offsets.npy,lengths.npy,tx_order.txt}
    -> per-tx RNA-seq coverage total (the ONLY condition-specific model input)

Outputs (results/gse120762/):
  per_condition.json         median/mean per-tx Pearson(pred,obs) + localization AUROC, per condition
  nt_vs_lps_differential.tsv  per shared tx: obs/pred/rna totals + CPM + log2FC(LPS/NT) each
  nt_vs_lps_summary.json      correlations: pred-dFC vs obs-dFC, rna-dFC vs obs-dFC, pred vs rna

Key interpretation: because sequence+ORF are identical across conditions, pred differential is
driven entirely by the RNA-seq differential; so pred-dFC ~ rna-dFC, and pred-dFC vs obs-dFC
measures how much of the LPS translational program is transcriptionally encoded.
No scipy dependency (numpy-only Pearson/Spearman).
"""
import json, os
import numpy as np

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
G = f"{NEW}/proteogenomics/data/gse120762"
OUT = f"{NEW}/results/gse120762"
os.makedirs(OUT, exist_ok=True)
MIN_PSITES = 50   # a tx is "scorable" for the profile Pearson if observed P-sites >= this


def pearson(a, b):
    a = a.astype(np.float64); b = b.astype(np.float64)
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return pearson(ra, rb)


def auroc(scores, labels):
    """Rank-based AUROC (per-nt frame-0 localization). labels in {0,1}."""
    labels = labels.astype(bool)
    n1 = labels.sum(); n0 = labels.size - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    order = np.argsort(scores)
    ranks = np.empty(scores.size, dtype=np.float64); ranks[order] = np.arange(1, scores.size + 1)
    # average ties (frame scores are near-continuous; cheap tie handling)
    return float((ranks[labels].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def load_pred(cond):
    d = np.load(f"{G}/pred_{cond}/pred_profiles.npz", allow_pickle=True)
    tx = np.array([str(x) for x in d["tx_ids"]])
    L = d["lengths"].astype(np.int64)
    off = np.zeros(len(L) + 1, dtype=np.int64); off[1:] = np.cumsum(L)
    return dict(tx=tx, L=L, off=off, pred=d["pred_flat"], obs=d["obs_flat"])


def load_rna_total(cond, tx_order):
    """Per-tx RNA-seq coverage total, aligned to tx_order (the pred tx order)."""
    pk = f"{NEW}/data/packed_heldout_mouse_gse120762_{cond}"
    cov = np.load(f"{pk}/coverage.npy", mmap_mode="r")
    Lp = np.load(f"{pk}/lengths.npy")
    offp = np.zeros(len(Lp) + 1, dtype=np.int64); offp[1:] = np.cumsum(Lp.astype(np.int64))
    ptx = [l.split("\t")[0].strip() for l in open(f"{pk}/tx_order.txt")]
    idx = {t: i for i, t in enumerate(ptx)}
    tot = np.zeros(len(tx_order), dtype=np.float64)
    for j, t in enumerate(tx_order):
        i = idx.get(t)
        if i is not None:
            tot[j] = float(np.asarray(cov[offp[i]:offp[i + 1]]).sum())
    return tot


per_cond = {}
store = {}
for cond in ["nt", "lps"]:
    P = load_pred(cond)
    store[cond] = P
    tx, L, off, pred, obs = P["tx"], P["L"], P["off"], P["pred"], P["obs"]
    rs, loc = [], []
    obs_tot = np.zeros(len(tx)); pred_tot = np.zeros(len(tx))
    for i in range(len(tx)):
        s, e = off[i], off[i + 1]
        o = obs[s:e]; p = pred[s:e]
        ot = int(o.sum()); obs_tot[i] = ot; pred_tot[i] = float(p.sum())
        if ot >= MIN_PSITES:
            r = pearson(p, o)
            if not np.isnan(r):
                rs.append(r)
            # localization: per-nt frame-0 of the ORF -- proxy = obs nonzero as "footprint present"
            lab = (o > 0).astype(np.int8)
            a = auroc(p, lab)
            if not np.isnan(a):
                loc.append(a)
    rs = np.array(rs); loc = np.array(loc)
    P["obs_tot"] = obs_tot; P["pred_tot"] = pred_tot
    per_cond[cond] = dict(
        n_tx=int(len(tx)), n_scorable=int(len(rs)),
        profile_pearson_median=float(np.median(rs)), profile_pearson_mean=float(rs.mean()),
        frac_r_gt_0p5=float((rs > 0.5).mean()), frac_r_gt_0p7=float((rs > 0.7).mean()),
        footprint_auroc_median=float(np.median(loc)),
        obs_psites_total=int(obs_tot.sum()),
    )
    print(f"[{cond}] scorable={len(rs)}  median profile r={np.median(rs):.3f}  "
          f"frac r>0.5={np.mean(rs > 0.5):.3f}  median footprint AUROC={np.median(loc):.3f}")

json.dump(per_cond, open(f"{OUT}/per_condition.json", "w"), indent=2)

# ---- NT vs LPS differential (shared tx) ----------------------------------------
tx_nt, tx_lps = store["nt"]["tx"], store["lps"]["tx"]
int_tx = np.array(sorted(set(tx_nt) & set(tx_lps)))
i_nt = {t: j for j, t in enumerate(tx_nt)}; i_lps = {t: j for j, t in enumerate(tx_lps)}
rna_nt = load_rna_total("nt", int_tx)
rna_lps = load_rna_total("lps", int_tx)


def cpm(x):
    x = np.asarray(x, dtype=np.float64); s = x.sum()
    return x / s * 1e6 if s > 0 else x


obs_nt = np.array([store["nt"]["obs_tot"][i_nt[t]] for t in int_tx])
obs_lps = np.array([store["lps"]["obs_tot"][i_lps[t]] for t in int_tx])
pred_nt = np.array([store["nt"]["pred_tot"][i_nt[t]] for t in int_tx])
pred_lps = np.array([store["lps"]["pred_tot"][i_lps[t]] for t in int_tx])

ocn, ocl = cpm(obs_nt), cpm(obs_lps)
pcn, pcl = cpm(pred_nt), cpm(pred_lps)
rcn, rcl = cpm(rna_nt), cpm(rna_lps)
l2 = lambda a, b: np.log2((b + 1.0) / (a + 1.0))
d_obs, d_pred, d_rna = l2(ocn, ocl), l2(pcn, pcl), l2(rcn, rcl)

# restrict differential correlations to tx with real Ribo signal in at least one condition
expr = (obs_nt + obs_lps) >= MIN_PSITES
summary = dict(
    n_shared_tx=int(len(int_tx)), n_expressed=int(expr.sum()),
    pred_dFC_vs_obs_dFC_pearson=pearson(d_pred[expr], d_obs[expr]),
    pred_dFC_vs_obs_dFC_spearman=spearman(d_pred[expr], d_obs[expr]),
    rna_dFC_vs_obs_dFC_pearson=pearson(d_rna[expr], d_obs[expr]),
    rna_dFC_vs_obs_dFC_spearman=spearman(d_rna[expr], d_obs[expr]),
    pred_dFC_vs_rna_dFC_pearson=pearson(d_pred[expr], d_rna[expr]),
    obs_total_nt=int(obs_nt.sum()), obs_total_lps=int(obs_lps.sum()),
)
json.dump(summary, open(f"{OUT}/nt_vs_lps_summary.json", "w"), indent=2)

with open(f"{OUT}/nt_vs_lps_differential.tsv", "w") as fh:
    fh.write("tx\tobs_nt\tobs_lps\tpred_nt\tpred_lps\trna_nt\trna_lps\t"
             "d_obs_log2FC\td_pred_log2FC\td_rna_log2FC\texpressed\n")
    for j, t in enumerate(int_tx):
        fh.write(f"{t}\t{obs_nt[j]:.0f}\t{obs_lps[j]:.0f}\t{pred_nt[j]:.1f}\t{pred_lps[j]:.1f}\t"
                 f"{rna_nt[j]:.0f}\t{rna_lps[j]:.0f}\t{d_obs[j]:.3f}\t{d_pred[j]:.3f}\t{d_rna[j]:.3f}\t"
                 f"{int(expr[j])}\n")

print("\n=== NT vs LPS differential (expressed shared tx, n=%d) ===" % summary["n_expressed"])
print("  pred-dFC vs obs-dFC : Pearson %.3f  Spearman %.3f"
      % (summary["pred_dFC_vs_obs_dFC_pearson"], summary["pred_dFC_vs_obs_dFC_spearman"]))
print("  rna-dFC  vs obs-dFC : Pearson %.3f  Spearman %.3f"
      % (summary["rna_dFC_vs_obs_dFC_pearson"], summary["rna_dFC_vs_obs_dFC_spearman"]))
print("  pred-dFC vs rna-dFC : Pearson %.3f  (should be high -- pred is RNA-driven)"
      % summary["pred_dFC_vs_rna_dFC_pearson"])
print(f"\nwrote {OUT}/{{per_condition.json,nt_vs_lps_summary.json,nt_vs_lps_differential.tsv}}")
