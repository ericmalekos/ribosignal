#!/usr/bin/env python3
"""GSE120762 NT vs LPS ORF-call differential: how the RiboCode ORF-call landscape shifts
between the two conditions -- by ORF_type (annotated/uORF/dORF/novel/...) and host
transcript_type (protein_coding vs lncRNA). This is the biologically meaningful "differential"
for an ORF-detection project: which ORFs get called, not quantitative fold-change of depth.

Inputs (already produced by the RiboCode held-out runs):
  data/heldout_psites/mouse_gse120762_{nt,lps}/mouse_gse120762_{nt,lps}_collapsed.txt

Outputs (results/gse120762/):
  orf_differential_by_type.tsv   ORF_type x condition counts, fraction, log2 fold-change
  orf_differential_lncrna.tsv    ORF_type x {protein_coding, lncRNA} per condition
  orf_call_overlap.json          shared / NT-only / LPS-only by ORF_ID, per ORF_type
  orf_differential_summary.txt   readable report incl. the depth control

Depth control: report each condition's total called ORFs AND total in-frame P-sites; a uniform
depth increase would scale every ORF_type equally, so disproportionate non-canonical enrichment
is the real signal.
"""
import json
import os
from collections import Counter, defaultdict

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
OUT = f"{NEW}/results/gse120762"
os.makedirs(OUT, exist_ok=True)
CONDS = ["nt", "lps"]
import math


def load(cond):
    """-> list of dicts (one per called ORF), keyed by header names."""
    f = f"{NEW}/data/heldout_psites/mouse_gse120762_{cond}/mouse_gse120762_{cond}_collapsed.txt"
    rows = []
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            rows.append(p)
    return rows, ci


data = {c: load(c) for c in CONDS}

# ---- 1. ORF_type composition + fold-change ------------------------------------
import numpy as np
by_type = {}
psite_total = {}
for c in CONDS:
    rows, ci = data[c]
    by_type[c] = Counter(r[ci["ORF_type"]] for r in rows)
    # library depth control = TOTAL P-sites on the universe (pack target), the cleanest measure
    # (total-in-called-ORFs would be circular: more calls -> more summed P-sites).
    psite_total[c] = int(np.load(f"{NEW}/data/packed_heldout_mouse_gse120762_{c}/"
                                 "target_counts.npy").sum())
depth_ratio = psite_total["lps"] / psite_total["nt"]

all_types = sorted(set(by_type["nt"]) | set(by_type["lps"]),
                   key=lambda t: -max(by_type["nt"][t], by_type["lps"][t]))
tot = {c: sum(by_type[c].values()) for c in CONDS}
with open(f"{OUT}/orf_differential_by_type.tsv", "w") as fh:
    fh.write("ORF_type\tNT\tLPS\tNT_frac\tLPS_frac\tlog2_LPS_over_NT\n")
    for t in all_types:
        n, l = by_type["nt"][t], by_type["lps"][t]
        lg = math.log2((l + 1) / (n + 1))
        fh.write(f"{t}\t{n}\t{l}\t{n/tot['nt']:.4f}\t{l/tot['lps']:.4f}\t{lg:.3f}\n")

# ---- 2. lncRNA vs protein_coding host, per ORF_type ---------------------------
def bt(tt):
    return "lncRNA" if tt == "lncRNA" else ("protein_coding" if tt == "protein_coding" else "other")


lnc = {c: defaultdict(Counter) for c in CONDS}   # lnc[c][ORF_type][biotype]
for c in CONDS:
    rows, ci = data[c]
    for r in rows:
        lnc[c][r[ci["ORF_type"]]][bt(r[ci["transcript_type"]])] += 1
with open(f"{OUT}/orf_differential_lncrna.tsv", "w") as fh:
    fh.write("ORF_type\tbiotype\tNT\tLPS\n")
    for t in all_types:
        for b in ("protein_coding", "lncRNA", "other"):
            fh.write(f"{t}\t{b}\t{lnc['nt'][t][b]}\t{lnc['lps'][t][b]}\n")

# lncRNA-hosted ORFs specifically (the "lncRNA ORF" differential)
lnc_orfs = {c: sum(v["lncRNA"] for v in lnc[c].values()) for c in CONDS}

# ---- 3. shared / NT-only / LPS-only by ORF_ID, per ORF_type -------------------
ids = {}
type_of = {}
for c in CONDS:
    rows, ci = data[c]
    ids[c] = set(r[ci["ORF_ID"]] for r in rows)
    for r in rows:
        type_of[r[ci["ORF_ID"]]] = r[ci["ORF_type"]]
shared = ids["nt"] & ids["lps"]
nt_only = ids["nt"] - ids["lps"]
lps_only = ids["lps"] - ids["nt"]


def type_breakdown(idset):
    return dict(Counter(type_of[i] for i in idset))


overlap = {
    "n_shared": len(shared), "n_nt_only": len(nt_only), "n_lps_only": len(lps_only),
    "shared_by_type": type_breakdown(shared),
    "nt_only_by_type": type_breakdown(nt_only),
    "lps_only_by_type": type_breakdown(lps_only),
}
json.dump(overlap, open(f"{OUT}/orf_call_overlap.json", "w"), indent=2)

# ---- 4. readable report -------------------------------------------------------
lines = []
lines.append("GSE120762 mouse BMDM +/- LPS -- observed RiboCode ORF-call differential\n")
lines.append(f"total ORFs called   NT={tot['nt']:,}   LPS={tot['lps']:,}   "
             f"(LPS/NT = {tot['lps']/tot['nt']:.2f}x)")
lines.append(f"library P-sites     NT={psite_total['nt']:,}   LPS={psite_total['lps']:,}   "
             f"(LPS/NT = {depth_ratio:.2f}x)  [depth control -- the real confounder]")
lines.append("")
lines.append(f"{'ORF_type':16s}{'NT':>8s}{'LPS':>8s}{'log2(LPS/NT)':>14s}")
for t in all_types:
    n, l = by_type["nt"][t], by_type["lps"][t]
    lines.append(f"{t:16s}{n:>8d}{l:>8d}{math.log2((l+1)/(n+1)):>14.2f}")
lines.append("")
lines.append(f"lncRNA-hosted ORFs   NT={lnc_orfs['nt']}   LPS={lnc_orfs['lps']}   "
             f"(LPS/NT = {lnc_orfs['lps']/max(lnc_orfs['nt'],1):.2f}x)")
lines.append("")
lines.append(f"ORF-call overlap (by ORF_ID):  shared={len(shared):,}  "
             f"NT-only={len(nt_only):,}  LPS-only={len(lps_only):,}")
lines.append(f"  LPS-only by type: {overlap['lps_only_by_type']}")
lines.append(f"  NT-only  by type: {overlap['nt_only_by_type']}")
lines.append("")
# depth-normalized read: does a type grow FASTER than the 2.7x LIBRARY-DEPTH increase?
# >1 = enriched beyond depth (candidate real LPS induction); <=1 = explained by depth / saturated.
# CAVEAT: ORF detection is superlinear in depth for weak ORFs, so even >1 needs a depth-MATCHED
# subsample (downsample LPS to NT depth + re-run RiboCode) to separate induction from sensitivity.
lines.append(f"enrichment vs {depth_ratio:.2f}x library-depth increase "
             f"(>1 = beyond depth; needs depth-matched subsample to confirm):")
for t in all_types:
    n, l = by_type["nt"][t], by_type["lps"][t]
    if n >= 10:
        lines.append(f"  {t:16s} {(l/n)/depth_ratio:>5.2f}x")
lines.append(f"  {'lncRNA-hosted':16s} {(lnc_orfs['lps']/max(lnc_orfs['nt'],1))/depth_ratio:>5.2f}x")
report = "\n".join(lines) + "\n"
open(f"{OUT}/orf_differential_summary.txt", "w").write(report)
print(report)
print(f"wrote {OUT}/orf_differential_*.{{tsv,json,txt}}")
