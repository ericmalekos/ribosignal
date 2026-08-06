#!/usr/bin/env python3
"""Check 7 / Objective O2: does the model predict mouse-liver ORF calls about as well as one real
liver Ribo-seq experiment predicts another's? The between-experiment agreement is the reproducibility
CEILING (biology + technical noise); if the model (sequence + RNA-seq, no Ribo) approaches it, that is
strong de novo evidence.

  ceiling = F1(Wang_obs, Janich_obs)                 two real liver Ribo datasets
  model   = F1(Wang_model_preddepth, Janich_obs)     model on Wang's RNA vs the OTHER real dataset
  (also model vs Wang_obs, and the with-Ribo pred_obsdepth, for context)

ORFs are matched on GENOMIC coordinates (gstart_gstop_len, the tail of RiboCode's ORF_ID), which is
annotation-version independent -- so a gene-version mismatch between the two pipelines does not break it.
Reported overall and split CDS (annotated) vs non-canonical.
"""
import argparse

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
ECH = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human"
WANG_OBS = f"{NEW}/data/heldout_psites/mouse_wang_liver/mouse_wang_liver_collapsed.txt"
JANICH_OBS = f"{ECH}/data/ribocode_mouse_liver/Liver_5samp/Liver_5samp_collapsed.txt"
WANG_MODEL = f"{NEW}/results/o2_liver/wang_nokozak/dropin/pred_preddepth_collapsed.txt"
WANG_MODEL_OBSDEPTH = f"{NEW}/results/o2_liver/wang_nokozak/dropin/pred_obsdepth_collapsed.txt"
OUT = f"{NEW}/results/o2_liver/o2_summary.txt"


def load(path):
    """-> {coord_key: ORF_type}. coord_key = gstart_gstop_len (annotation-version independent)."""
    d = {}
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


def f1(pred, ref):
    inter = len(pred & ref)
    p = inter / len(pred) if pred else 0.0
    r = inter / len(ref) if ref else 0.0
    return p, r, (2 * p * r / (p + r) if (p + r) else 0.0)


def cls(d, which):
    if which == "cds":
        return {k for k, t in d.items() if t == "annotated"}
    return {k for k, t in d.items() if t != "annotated"}


def load_opt(p):
    try:
        return load(p)
    except FileNotFoundError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=WANG_MODEL,
                    help="de-novo (pred_preddepth) call file for the model arm")
    ap.add_argument("--model-obsdepth", default=WANG_MODEL_OBSDEPTH, dest="model_obsdepth")
    ap.add_argument("--model-name", default="Wang_model", dest="model_name")
    ap.add_argument("--other-name", default="Janich_obs", dest="other_name")
    ap.add_argument("--swap", action="store_true",
                    help="SYMMETRIC arm: model built on JANICH RNA, so score it against WANG "
                         "observed (the other experiment). Mirrors the default arm.")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    wang = load(WANG_OBS)
    jan = load(JANICH_OBS)
    model = load_opt(a.model)         # None until the dump lands -> ceiling still computes
    model_obs = load_opt(a.model_obsdepth)
    have_model = model is not None
    # `other` = the experiment the model is scored AGAINST (must not be the one it was built on);
    # `self_` = the matched experiment, reported for context.
    other, self_ = (wang, jan) if a.swap else (jan, wang)

    head = f"Wang_obs={len(wang)}  Janich_obs={len(jan)}"
    head += (f"  {a.model_name}={len(model)} calls" if have_model
             else "  (model pending -- ceiling only)")
    L = ["O2: mouse-liver ORF-call reproducibility ceiling vs model (de novo)\n", head + "\n"]
    if a.swap:
        L.append("SYMMETRIC arm: model built on Janich RNA, scored against Wang observed.\n")
    for grp, sel in [("ALL", None), ("CDS", "cds"), ("non-canonical", "noncan")]:
        def S(d, sel=sel):
            return set(d) if sel is None else cls(d, "cds" if sel == "cds" else "nc")
        _, _, cf = f1(S(wang), S(jan))     # CEILING: two real experiments
        line = f"[{grp}]  ceiling F1(Wang,Janich)={cf:.3f}"
        if have_model:
            m = S(model)
            _, _, mf = f1(m, S(other))     # model vs the OTHER real dataset
            _, _, sf = f1(m, S(self_))     # model vs the matched dataset
            line += f"   model(deNovo)F1={mf:.3f}   model(self)F1={sf:.3f}"
            if cf:
                line += f"   [model/ceiling={mf / cf:.2f}]"
            if model_obs is not None:
                _, _, mof = f1(S(model_obs), S(other))
                line += f"   with-Ribo F1={mof:.3f}"
        L.append(line)
    report = "\n".join(L) + (
        "\n\nCeiling = agreement between two real liver Ribo-seq datasets (the bar). "
        "model/ceiling ~1.0 => model is as informative as a 2nd real experiment.\n")
    open(a.out, "w").write(report)
    print(report)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
