#!/usr/bin/env python3
"""Three-way mouse-liver ORF-call concordance: Wang vs Janich vs GSE243134, and the model against
each. Precision AND recall, both directions, by ORF class.

WHY THREE. The whole calibration argument rests on a ceiling estimated from ONE pair (Wang vs
Janich): CDS F1 0.735, non-canonical 0.505. With a single pair there is no way to tell whether 0.505
is a property of non-canonical ORF calling or a quirk of those two datasets. GSE243134 supplies a
third independent liver Ribo-seq (21 WT runs, 273 M unique RPF) and therefore two more pairs.

Its RNA-seq was rejected on library chemistry, but observed-vs-observed concordance needs no RNA-seq
at all -- only the Ribo-seq -- so the arm that was shelved with its RNA can still carry a ceiling
point. No model arm is built ON GSE243134 for the same reason (a model arm would need its
expression), so the model is scored AGAINST its calls using the Wang- and Janich-driven predictions.

BOTH DIRECTIONS ARE REPORTED for every pair, and that is the point of the design rather than
decoration: precision and recall swap when the roles swap, so an asymmetry between A->B and B->A is
a DEPTH difference (the deeper set calls more, so it recalls more of the shallower one and is less
precise against it), not a reproducibility limit. Reading only one direction confuses the two.

Keying is genomic -- (gene_id, ORF_gstop) -- because RiboCode's collapse keeps the longest ORF per
genomic stop, so which isoform carries a call depends on the transcript set. All three observed sets
come from full-annotation RiboCode runs on the same prepared vM38 annotation, so no gene restriction
is applied between them; the model arms ARE restricted to their own universe's gene space, and
n_ref is printed for every comparison.

Usage: threeway_liver_concordance.py [--min-len 90] [--pval 0.05] [--out <dir>]
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from compare_dropin_calls import NONCANON, load_calls, prf  # noqa: E402

ECH = ROOT.parent.parent / "experiments" / "biotype_probe" / "expression_context_human"
OBS = {
    "Wang": ROOT / "data/heldout_psites/mouse_wang_liver/mouse_wang_liver_collapsed.txt",
    "Janich": ECH / "data/ribocode_mouse_liver/Liver_5samp/Liver_5samp_collapsed.txt",
    "GSE243134": ROOT / "data/heldout_psites/mouse_gse243134_liver/mouse_gse243134_liver_collapsed.txt",
}
# BOTH arms per the standing two-arm rule. theta=0.02 is the CDS-anchored calibrated operating point
# from the O2 sweep (results/o2_liver/o2_poisson_summary.txt).
# RELEASED models only. The previous entries pointed at results/o2_liver/, which came from the
# pre-union Fibroblast-universe checkpoint and had no mamba4 arm -- see scripts/dump_liver_released.
LR = ROOT / "results/liver_released"
MODEL = {
    "attn/Wang th=1":      LR / "attn_mouse_wang_liver/pred_preddepth_collapsed.txt",
    "mamba4/Wang th=1":    LR / "mamba4_mouse_wang_liver/pred_preddepth_collapsed.txt",
    "attn/Janich th=1":    LR / "attn_mouse_janich_liver_decon/pred_preddepth_collapsed.txt",
    "mamba4/Janich th=1":  LR / "mamba4_mouse_janich_liver_decon/pred_preddepth_collapsed.txt",
    "attn/GSE243134 th=1": LR / "attn_mouse_gse243134_liver/pred_preddepth_collapsed.txt",
    "mamba4/GSE243134 th=1": LR / "mamba4_mouse_gse243134_liver/pred_preddepth_collapsed.txt",
}
CLASSES = {"CDS": {"annotated"}, "non-canonical": NONCANON}

# Expression restriction. The O2 dumps used --tx_list and so bypassed the min_signal gate every
# other eval applies: the model predicts over the whole 221,835-tx annotation while the observed
# sets span ~17k transcripts. Scoring unrestricted therefore charges the model for calling ORFs on
# transcripts no reference could ever contain, and the resulting precision (0.02 on non-canonical) is
# an artifact of the transcript space, not a property of the model. Standing rule:
# feedback_orf_call_transcript_space.
#
# Janich's quant is the DECONTAMINATED one (mean TPM over 7 replicates): its shipped merged quant put
# 52-70% of TPM into two structural transcripts, which suppressed ~4,900 real transcripts below the
# TPM >= 1 cut. GSE243134 contributes no expression set -- its totalRNA was rejected on chemistry --
# so the shared space is Wang AND Janich, and GSE243134 is scored inside it.
SALMON_WANG = ROOT / "data/salmon_quant_wang_liver/wang_merged_quant.sf"
JANICH_DECON = ROOT / "data/salmon_quant_janich_liver_decontam"
TX2GENE = ROOT / "data/tx2biotype_mouse.tsv"


def salmon_tpm(path):
    """{bare versioned tx id: TPM}; quant.sf carries full pipe-delimited GENCODE headers."""
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ni, ti = hdr.index("Name"), hdr.index("TPM")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            d[f[ni].split("|")[0]] = float(f[ti])
    return d


def mean_tpm_dir(d):
    """Mean TPM across every quant.sf under `d` (missing in a replicate counts as 0)."""
    quants = sorted(Path(d).glob("*/quant.sf"))
    acc, n = {}, len(quants)
    for q in quants:
        for tx, v in salmon_tpm(q).items():
            acc[tx] = acc.get(tx, 0.0) + v
    return {tx: v / n for tx, v in acc.items()}, n


def expressed_genes(min_tpm):
    """Genes with >= min_tpm in BOTH Wang and (decontaminated) Janich."""
    w = {t for t, v in salmon_tpm(SALMON_WANG).items() if v >= min_tpm}
    jt, njan = mean_tpm_dir(JANICH_DECON)
    j = {t for t, v in jt.items() if v >= min_tpm}
    shared = w & j
    t2g = {}
    with open(TX2GENE) as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            t2g[f[ci["tx_id"]]] = f[ci["gene_id"]]
    genes = {t2g[t] for t in shared if t in t2g}
    return genes, {"wang_tx": len(w), "janich_tx": len(j), "janich_reps": njan,
                   "shared_tx": len(shared), "shared_genes": len(genes)}


def by_class(pred, real, want):
    pk = {k for k in pred if pred[k]["type"] in want}
    rk = {k for k in real if real[k]["type"] in want}
    return prf(pk, rk)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-len", type=int, default=90)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--min-tpm", type=float, default=1.0,
                    help="expression cut for the shared Wang+Janich restriction; 0 = unrestricted")
    ap.add_argument("--out", default=str(ROOT / "results/o2_liver"))
    a = ap.parse_args()
    keep_genes, xinfo = (None, {})
    if a.min_tpm > 0:
        keep_genes, xinfo = expressed_genes(a.min_tpm)
        print(f"expression restriction (TPM >= {a.min_tpm}, Wang AND decontaminated Janich):")
        for k, v in xinfo.items():
            print(f"  {k:<16}{v:>10,}")
    kw = dict(key_mode="genomic", max_pval=a.pval, min_len=a.min_len, keep_genes=keep_genes)

    calls, missing = {}, []
    for name, f in {**OBS, **MODEL}.items():
        if not Path(f).exists():
            missing.append(f"{name}: {f}")
            continue
        calls[name] = load_calls(f, **kw)
    if missing:
        print("MISSING (absent from the table, not imputed):")
        for m in missing:
            print(f"  {m}")

    rows = []
    # ---- observed vs observed: every ordered pair, so precision/recall asymmetry is visible ----
    for a_, b_ in itertools.permutations([k for k in OBS if k in calls], 2):
        r = {"kind": "observed", "pred": a_, "ref": b_, "overall": prf(set(calls[a_]), set(calls[b_]))}
        for cl, want in CLASSES.items():
            r[cl] = by_class(calls[a_], calls[b_], want)
        rows.append(r)
    # ---- model vs each observed set ----
    for m in [k for k in MODEL if k in calls]:
        for b_ in [k for k in OBS if k in calls]:
            r = {"kind": "model", "pred": m, "ref": b_, "overall": prf(set(calls[m]), set(calls[b_]))}
            for cl, want in CLASSES.items():
                r[cl] = by_class(calls[m], calls[b_], want)
            rows.append(r)
    # ---- model vs the three-way CORE (called independently by all three experiments) ----
    obs_present = [k for k in OBS if k in calls]
    if len(obs_present) == 3:
        core_keys = set(calls[obs_present[0]]) & set(calls[obs_present[1]]) & set(calls[obs_present[2]])
        core = {k: calls[obs_present[0]][k] for k in core_keys}
        calls["CORE(all 3)"] = core
        for m in [k for k in MODEL if k in calls]:
            r = {"kind": "model-vs-core", "pred": m, "ref": "CORE(all 3)",
                 "overall": prf(set(calls[m]), core_keys)}
            for cl, want in CLASSES.items():
                r[cl] = by_class(calls[m], core, want)
            rows.append(r)

    print(f"\ncall-set sizes (pval<={a.pval}, min_len={a.min_len}, genomic keying)")
    for k, v in calls.items():
        print(f"  {k:<22}{len(v):>8,}")

    for kind, title in (("observed", "OBSERVED vs OBSERVED  (the ceiling, both directions)"),
                        ("model", "MODEL vs each observed set"),
                        ("model-vs-core", "MODEL vs the three-way CORE")):
        sel = [r for r in rows if r["kind"] == kind]
        if not sel:
            continue
        print(f"\n=== {title} ===")
        print(f"{'pred':<22}{'ref':<14}{'n_pred':>8}{'n_ref':>8}"
              f"{'P':>7}{'R':>7}{'F1':>7} | {'CDS P':>7}{'CDS R':>7}{'CDS F1':>8}"
              f" | {'nc P':>7}{'nc R':>7}{'nc F1':>7}")
        for r in sel:
            o, c, n = r["overall"], r["CDS"], r["non-canonical"]
            print(f"{r['pred']:<22}{r['ref']:<14}{o['n_pred']:>8,}{o['n_real']:>8,}"
                  f"{o['precision']:>7.3f}{o['recall']:>7.3f}{o['f1']:>7.3f}"
                  f" | {c['precision']:>7.3f}{c['recall']:>7.3f}{c['f1']:>8.3f}"
                  f" | {n['precision']:>7.3f}{n['recall']:>7.3f}{n['f1']:>7.3f}")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "threeway_liver_concordance.json").write_text(json.dumps(
        {"min_len": a.min_len, "pval": a.pval, "min_tpm": a.min_tpm, "restriction": xinfo,
         "sizes": {k: len(v) for k, v in calls.items()}, "rows": rows, "missing": missing}, indent=2))
    print(f"\nwrote {out}/threeway_liver_concordance.json")


if __name__ == "__main__":
    main()
