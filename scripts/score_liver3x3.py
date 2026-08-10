#!/usr/bin/env python3
"""Score the mouse-liver 3x3 factorial: precision / recall / F1 across experiments and models.

The grid is 3 Ribo-seq datasets x 3 RNA-seq inputs, all pooled from BAMs through ONE pipeline onto
one shared universe (22,974 transcripts), so a difference in the table is a difference between
experiments rather than between code paths.

Two questions are answered separately:

  A. Ribo-vs-Ribo. How well does each dataset's OBSERVED call set reproduce each other dataset's?
     This is the ceiling: no model is involved, and it bounds what any predictor can be expected to
     score against a given reference.

  B. Model-vs-Ribo. For each (model, RNA input), how well do the predicted calls reproduce each
     Ribo dataset's observed calls? Reported for both prediction arms:
       pred_obsdepth   predicted shape scaled by the reference dataset's observed depth
       pred_preddepth  fully standalone -- no Ribo-seq at inference

Filtering comes from `compare_dropin_calls.build_loader`, the same code path behind every other
drop-in number in this project. Reimplementing it once produced a table that disagreed with the
project convention by 0.14 F1, so it is imported, never copied.

Keying is GENOMIC (gene_id, ORF_gstop). Transcript keying would charge a model for picking a
different representative isoform of the same ORF, which is a collapse artefact and not a biological
disagreement.

Every row states n_ref, the size of the reference call set it was scored against
(feedback_orf_call_transcript_space): a precision computed against an unstated transcript space is
not interpretable.
"""
import argparse
import itertools
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_dropin_calls import build_loader  # noqa: E402

DATASETS = ["janich", "gse243134", "wang"]
MODELS = ["attn", "mamba4"]


def prf(pred, ref, subset=None):
    """Precision / recall / F1 of pred against ref, optionally restricted to one ORF class."""
    if subset:
        pred = {k: v for k, v in pred.items() if v["type"] in subset}
        ref = {k: v for k, v in ref.items() if v["type"] in subset}
    if not ref or not pred:
        return dict(precision=0.0, recall=0.0, f1=0.0, n_pred=len(pred), n_ref=len(ref), tp=0)
    tp = len(set(pred) & set(ref))
    p = tp / len(pred)
    r = tp / len(ref)
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return dict(precision=p, recall=r, f1=f, n_pred=len(pred), n_ref=len(ref), tp=tp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="results/mouse_liver_3x3")
    ap.add_argument("--out", required=True, help="output directory for TSVs")
    ap.add_argument("--classes", default="annotated",
                    help="comma-separated ORF classes for the restricted table (default annotated)")
    # build_loader defaults to the HUMAN tx2biotype; this factorial is mouse. It errors loudly on a
    # zero-overlap map rather than silently scoring nothing, but the right map still has to be passed.
    ap.add_argument("--tx2gene",
                    default=str(pathlib.Path(__file__).resolve().parents[1]
                                / "data" / "tx2biotype_mouse.tsv"),
                    help="tx->gene map matching the assembly (mouse for this factorial)")
    a = ap.parse_args()

    res = pathlib.Path(a.results)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cls = tuple(a.classes.split(","))

    # ---- load every cell's calls under its own dump's filtering ----
    # `real` depends only on the Ribo source, so it is read once per Ribo dataset from the diagonal
    # cell. The off-diagonal copies were verified identical by construction (same P-sites), and that
    # is asserted below rather than assumed.
    real, loaders = {}, {}
    for ds in DATASETS:
        d = res / f"mamba4_ribo-{ds}_rna-{ds}"
        if not (d / "pred_profiles.npz").exists():
            print(f"missing dump for {ds}: {d}", file=sys.stderr)
            return 1
        lc, _, _ = build_loader(d / "pred_profiles.npz", tx2gene=a.tx2gene)
        loaders[ds] = lc
        real[ds] = lc(d / "real_collapsed.txt", is_pred=False)
        print(f"  real[{ds}]: {len(real[ds]):,} calls")

    # Consistency: the same Ribo source must give the same observed calls whatever RNA was in the pack.
    incons = []
    for ds, rna in itertools.product(DATASETS, DATASETS):
        d = res / f"mamba4_ribo-{ds}_rna-{rna}"
        f = d / "real_collapsed.txt"
        if not f.exists():
            continue
        got = loaders[ds](f, is_pred=False)
        if set(got) != set(real[ds]):
            incons.append(f"ribo={ds} rna={rna}: {len(got)} vs {len(real[ds])}")
    print(f"  consistency: {'OK -- real identical across RNA columns' if not incons else 'FAILED'}")
    for m in incons:
        print(f"    {m}")

    # ---- A. Ribo vs Ribo (the ceiling) ----
    rows = []
    for ref, pred in itertools.product(DATASETS, DATASETS):
        if ref == pred:
            continue
        s = prf(real[pred], real[ref])
        sc = prf(real[pred], real[ref], subset=cls)
        rows.append(dict(reference=ref, compared=pred, **{k: s[k] for k in
                    ("precision", "recall", "f1", "n_pred", "n_ref", "tp")},
                    f1_annotated=sc["f1"], n_ref_annotated=sc["n_ref"]))
    write_tsv(out / "ribo_vs_ribo.tsv", rows)

    # ---- B. model vs Ribo, both arms ----
    rows = []
    for model, rna, ref in itertools.product(MODELS, DATASETS, DATASETS):
        d = res / f"{model}_ribo-{ref}_rna-{rna}"
        if not (d / "pred_profiles.npz").exists():
            continue
        lc, _, _ = build_loader(d / "pred_profiles.npz", tx2gene=a.tx2gene)
        for arm in ("pred_obsdepth", "pred_preddepth"):
            f = d / f"{arm}_collapsed.txt"
            if not f.exists():
                continue
            pc = lc(f, is_pred=True)
            s = prf(pc, real[ref])
            sc = prf(pc, real[ref], subset=cls)
            rows.append(dict(model=model, rna_input=rna, ribo_reference=ref, arm=arm,
                             **{k: s[k] for k in
                                ("precision", "recall", "f1", "n_pred", "n_ref", "tp")},
                             precision_annotated=sc["precision"], recall_annotated=sc["recall"],
                             f1_annotated=sc["f1"], n_ref_annotated=sc["n_ref"]))
    write_tsv(out / "model_vs_ribo.tsv", rows)

    # ---- matched vs mismatched RNA, the question the factorial exists to answer ----
    summ = []
    for model, ref in itertools.product(MODELS, DATASETS):
        for arm in ("pred_obsdepth", "pred_preddepth"):
            sel = [r for r in rows if r["model"] == model and r["ribo_reference"] == ref
                   and r["arm"] == arm]
            if not sel:
                continue
            m = [r for r in sel if r["rna_input"] == ref]
            x = [r for r in sel if r["rna_input"] != ref]
            if not m or not x:
                continue
            summ.append(dict(model=model, ribo_reference=ref, arm=arm,
                             f1_matched_rna=m[0]["f1"],
                             f1_mismatched_rna_mean=sum(r["f1"] for r in x) / len(x),
                             delta=m[0]["f1"] - sum(r["f1"] for r in x) / len(x),
                             n_ref=m[0]["n_ref"]))
    write_tsv(out / "matched_vs_mismatched_rna.tsv", summ)

    (out / "scoring_meta.json").write_text(json.dumps(dict(
        universe_tx=22974, key="genomic (gene_id, ORF_gstop)",
        filtering="compare_dropin_calls.build_loader (pval<=0.05, ORF>=90nt, enrichment>=0.5x)",
        classes_restricted=list(cls), consistency_ok=not incons), indent=2))
    print(f"\nwrote 3 TSVs + scoring_meta.json to {out}")
    return 0


def write_tsv(path, rows):
    if not rows:
        print(f"  (no rows for {path.name})")
        return
    hdr = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])
                               for k in hdr) + "\n")
    print(f"  wrote {path.name}  ({len(rows)} rows)")


if __name__ == "__main__":
    sys.exit(main())
