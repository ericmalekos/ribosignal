#!/usr/bin/env python3
"""Are the model's "extra" NT ORF calls (predicted standalone but absent from NT's own shallow
observed calls) spurious, or real ORFs that shallow NT Ribo-seq under-detects? Test: check them
against the DEEP LPS observed calls (same cell type, ~2.7x deeper) as a higher-quality reference.

An extra call that reappears in LPS-real is a real, translatable ORF the shallow NT data simply
missed -- i.e. the model recovered it, it is not an over-call. The validated fraction is a LOWER
bound on how many extras are real (some unvalidated ones may be NT-specific or below detection
even at LPS depth).

Inputs (ribocode_dropin outputs): proteogenomics/data/gse120762/pred_{nt,lps}/dropin/*_collapsed.txt
Output: results/gse120762/overcall_check.txt
"""
import os
import sys
from collections import Counter

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
G = f"{NEW}/proteogenomics/data/gse120762"
OUT = f"{NEW}/results/gse120762"
SUF = sys.argv[1] if len(sys.argv) > 1 else ""   # "" = deployed, "_nokozak" = mm1 no-kozak model


def load(cond, variant):
    f = f"{G}/pred_{cond}{SUF}/dropin/{variant}_collapsed.txt"
    d = {}
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[p[ci["ORF_ID"]]] = p[ci["ORF_type"]]
    return d


nt_pred = load("nt", "pred_preddepth")   # standalone NT prediction (the "over-caller")
nt_real = load("nt", "real")             # shallow NT observed (14M P-sites)
lps_real = load("lps", "real")           # DEEP LPS observed (38M P-sites) -- the reference
lps_ids = set(lps_real)
nt_real_ids = set(nt_real)

extra = {i: t for i, t in nt_pred.items() if i not in nt_real_ids}   # predicted, not in shallow NT

L = []
L.append("GSE120762: are the model's extra NT calls real? (validate vs deep LPS observed)\n")
L.append(f"NT standalone predicted : {len(nt_pred):,}")
L.append(f"  confirmed by NT-real  : {len(set(nt_pred) & nt_real_ids):,}")
L.append(f"  EXTRA (not in NT-real): {len(extra):,}")
val = sum(1 for i in extra if i in lps_ids)
L.append(f"    of which in DEEP LPS-real: {val:,}  ({val/len(extra)*100:.1f}%)  <- lower bound on 'real'")
L.append("")
# baseline: how often does an NT-real call itself reappear in LPS-real? (cross-condition ceiling)
base = len(nt_real_ids & lps_ids)
L.append(f"baseline (reproducibility ceiling): NT-real calls also in LPS-real = "
         f"{base:,}/{len(nt_real_ids):,} ({base/len(nt_real_ids)*100:.1f}%)")
L.append("  -> if the extras validate near this rate, they are as real as the confirmed calls")
L.append("")
L.append(f"{'ORF_type':14s}{'extra':>8s}{'in_LPS':>8s}{'valid%':>8s}{'   NTreal->LPS%':>16s}")
for t in sorted(set(extra.values()), key=lambda t: -Counter(extra.values())[t]):
    ex_t = [i for i, tt in extra.items() if tt == t]
    v = sum(1 for i in ex_t if i in lps_ids)
    # baseline per type: NT-real calls of this type that reappear in LPS-real
    ntr_t = [i for i, tt in nt_real.items() if tt == t]
    b = sum(1 for i in ntr_t if i in lps_ids)
    bp = f"{b/len(ntr_t)*100:.0f}% (n={len(ntr_t)})" if ntr_t else "n=0"
    L.append(f"{t:14s}{len(ex_t):>8d}{v:>8d}{v/len(ex_t)*100:>7.1f}%{bp:>16s}")
report = "\n".join(L) + "\n"
os.makedirs(OUT, exist_ok=True)
open(f"{OUT}/overcall_check{SUF}.txt", "w").write(report)
print(f"[model: {'mm1 no-kozak' if SUF else 'deployed heuristic'}]")
print(report)
