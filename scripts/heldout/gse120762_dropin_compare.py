#!/usr/bin/env python3
"""Compare the model's DROP-IN ORF calls to the observed calls, per GSE120762 condition, and
ask whether the model reproduces (a) the ORF-type composition within each condition and (b) the
NT->LPS shift in that composition.

All variants come from ribocode_dropin.py on the SAME universe tx set + SAME caller/params, so
'real' (RiboCode on observed P-sites) is the fair reference for the predicted variants:
  real            observed P-sites          -> the reference calls (universe-restricted)
  pred_obsdepth   predicted SHAPE x real depth   -> isolates profile-shape quality
  pred_preddepth  predicted shape x predicted total -> standalone (no real Ribo data)

Inputs:  proteogenomics/data/gse120762/pred_{nt,lps}/dropin/{real,pred_obsdepth,pred_preddepth}_collapsed.txt
Outputs (results/gse120762/):
  dropin_composition.tsv   ORF_type x condition x variant counts
  dropin_agreement.tsv     per-condition per-type detection precision/recall/F1 (pred vs real, by ORF_ID)
  dropin_compare_summary.txt  readable report
"""
import math
import os
from collections import Counter, defaultdict

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
G = f"{NEW}/proteogenomics/data/gse120762"
OUT = f"{NEW}/results/gse120762"
os.makedirs(OUT, exist_ok=True)
CONDS = ["nt", "lps"]
VARIANTS = ["real", "pred_obsdepth", "pred_preddepth"]


def load_calls(cond, variant):
    """-> {ORF_ID: ORF_type} for one condition/variant, or None if the file is absent."""
    f = f"{G}/pred_{cond}/dropin/{variant}_collapsed.txt"
    if not os.path.exists(f):
        return None
    d = {}
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[p[ci["ORF_ID"]]] = p[ci["ORF_type"]]
    return d


calls = {(c, v): load_calls(c, v) for c in CONDS for v in VARIANTS}
missing = [k for k, v in calls.items() if v is None]
if missing:
    raise SystemExit(f"drop-in outputs not ready yet: {missing}")

all_types = sorted({t for d in calls.values() for t in d.values()},
                   key=lambda t: -sum(sum(1 for x in calls[(c, 'real')].values() if x == t) for c in CONDS))

# ---- composition: counts by ORF_type x condition x variant ---------------------
with open(f"{OUT}/dropin_composition.tsv", "w") as fh:
    fh.write("ORF_type\t" + "\t".join(f"{c}_{v}" for c in CONDS for v in VARIANTS) + "\n")
    for t in all_types:
        row = [str(sum(1 for x in calls[(c, v)].values() if x == t)) for c in CONDS for v in VARIANTS]
        fh.write(t + "\t" + "\t".join(row) + "\n")

# ---- detection agreement (pred vs real, by ORF_ID) per condition, per type -----
def prf(pred_ids, real_ids):
    inter = len(pred_ids & real_ids)
    p = inter / len(pred_ids) if pred_ids else 0.0
    r = inter / len(real_ids) if real_ids else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f


agree = {}
with open(f"{OUT}/dropin_agreement.tsv", "w") as fh:
    fh.write("condition\tvariant\tORF_type\tn_pred\tn_real\tprecision\trecall\tF1\n")
    for c in CONDS:
        real = calls[(c, "real")]
        real_ids_all = set(real)
        for v in ("pred_obsdepth", "pred_preddepth"):
            pred = calls[(c, v)]
            pred_ids_all = set(pred)
            p, r, f = prf(pred_ids_all, real_ids_all)
            agree[(c, v, "ALL")] = (len(pred), len(real), p, r, f)
            fh.write(f"{c}\t{v}\tALL\t{len(pred)}\t{len(real)}\t{p:.3f}\t{r:.3f}\t{f:.3f}\n")
            for t in all_types:
                pid = {i for i, tt in pred.items() if tt == t}
                rid = {i for i, tt in real.items() if tt == t}
                if not pid and not rid:
                    continue
                p2, r2, f2 = prf(pid, rid)
                fh.write(f"{c}\t{v}\t{t}\t{len(pid)}\t{len(rid)}\t{p2:.3f}\t{r2:.3f}\t{f2:.3f}\n")

# ---- readable report -----------------------------------------------------------
def comp(cond, variant):
    return Counter(calls[(cond, variant)].values())


L = []
L.append("GSE120762 drop-in: model ORF calls vs observed, per condition\n")
L.append(f"{'ORF_type':14s}" + "".join(f"{c+'/'+v[:8]:>16s}" for c in CONDS for v in ("real", "pred_preddepth")))
for t in all_types:
    cells = ""
    for c in CONDS:
        cr = sum(1 for x in calls[(c, 'real')].values() if x == t)
        cp = sum(1 for x in calls[(c, 'pred_preddepth')].values() if x == t)
        cells += f"{cr:>8d}{cp:>8d}"
    L.append(f"{t:14s}{cells}")
L.append("")
L.append("detection F1 (pred vs real, all ORF types, by ORF_ID):")
for c in CONDS:
    for v in ("pred_obsdepth", "pred_preddepth"):
        n_pred, n_real, p, r, f = agree[(c, v, "ALL")]
        L.append(f"  {c:4s} {v:14s}  P={p:.3f} R={r:.3f} F1={f:.3f}  (n_pred={n_pred} n_real={n_real})")
L.append("")
# NT->LPS shift fidelity: does pred reproduce the direction of the real per-type change?
L.append("NT->LPS log2 fold-change per ORF_type (does the model track the observed shift?):")
L.append(f"  {'ORF_type':14s}{'real':>10s}{'pred_preddepth':>16s}")
for t in all_types:
    def lg(v):
        n = sum(1 for x in calls[('nt', v)].values() if x == t)
        l = sum(1 for x in calls[('lps', v)].values() if x == t)
        return math.log2((l + 1) / (n + 1))
    L.append(f"  {t:14s}{lg('real'):>10.2f}{lg('pred_preddepth'):>16.2f}")
report = "\n".join(L) + "\n"
open(f"{OUT}/dropin_compare_summary.txt", "w").write(report)
print(report)
print(f"wrote {OUT}/dropin_composition.tsv, dropin_agreement.tsv, dropin_compare_summary.txt")
