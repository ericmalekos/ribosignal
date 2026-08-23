#!/usr/bin/env python3
"""Over-call validation: are the model's experiment-unsupported ORF calls false positives?

THE QUESTION. The mouse-liver UpSet (figures/P9_upset_mouse3x3) shows ORFs called by the model and by
NO single Ribo-seq experiment. Two explanations are indistinguishable from that plot: they are model
false positives, or they are real ORFs no single experiment was deep enough to detect. Every observed
reference in the 3x3 is one dataset deep, so none of them can arbitrate.

THE ARBITER. `results/merged_liver_ribocode` is a JOINT RiboCode call over all 28 canonically-aligned
mouse-liver libraries (Janich 5 + GSE243134 21 + Wang 2) -- 26,388 raw calls, +77% to +102% over any
single dataset. An ORF the model called that even the merged reference misses is a much stronger
false-positive candidate.

THE CONTROL THAT MAKES IT INTERPRETABLE. A raw recovery rate means nothing on its own: the merged
reference is deeper than any single dataset, so it recovers more of EVERYTHING. The comparison that
carries information is the model's unsupported calls against the OBSERVED datasets' equally
unsupported calls -- ORFs called by exactly one experiment and missed by the other two. Those are
known-real (an experiment saw them) and equally unsupported by the rest of the 3x3, so they measure
how often "unsupported by the other arms" simply means "too deep for them". Reporting the model's
number without this control would attribute the depth effect to the model.

  recovery(model-unique)  ~ recovery(experiment-unique)  -> the extras behave like real, deep ORFs
  recovery(model-unique) << recovery(experiment-unique)  -> the extras are enriched for FPs

Both are reported per ORF class, because lncRNA (`novel`) and uORF calls fail in opposite directions
and a pooled rate hides that.

Filtering and keying come from `compare_dropin_calls.build_loader` -- the merged reference is
restricted to the SAME transcript space as every other number here (feedback_orf_call_transcript_space),
and matched on (gene_id, ORF_gstop) so an isoform-representative difference is not scored as a
disagreement.

  python3 scripts/score_merged_liver_overcall.py --out results/merged_liver_ribocode/overcall
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_dropin_calls import build_loader  # noqa: E402

NEW = pathlib.Path(__file__).resolve().parents[1]
DATASETS = ["janich", "gse243134", "wang"]
MODELS = ["mamba4", "attn"]
MODEL_RNA = "gse243134"       # matches figures/P9_upset_mouse3x3
ARM = "pred_preddepth"        # the standalone arm -- what the application actually deploys
CLASSES = ["annotated", "novel", "uORF", "internal", "Overlap_uORF", "dORF", "Overlap_dORF"]


def cell(model, ribo, rna):
    return NEW / f"results/mouse_liver_3x3_canon/{model}_ribo-{ribo}_rna-{rna}"


def rate(sub, merged):
    """Fraction of `sub` that the merged reference recovers."""
    n = len(sub)
    hit = len(set(sub) & set(merged))
    return dict(n=n, recovered=hit, rate=hit / n if n else float("nan"))


def by_class(sub, merged, calls):
    out = {}
    for c in CLASSES:
        s = {k for k in sub if calls[k]["type"] == c}
        if s:
            out[c] = rate(s, merged)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged",
                    default=str(NEW / "results/merged_liver_ribocode/merged_liver_collapsed.txt"))
    ap.add_argument("--out", default=str(NEW / "results/merged_liver_ribocode/overcall"))
    ap.add_argument("--tx2gene", default=str(NEW / "data/tx2biotype_mouse.tsv"))
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    base = cell("mamba4", MODEL_RNA, MODEL_RNA)
    lc, keep_tx, keep_genes = build_loader(str(base / "pred_profiles.npz"), key="genomic",
                                           tx2gene=a.tx2gene)
    print(f"transcript space: {len(keep_tx):,} tx / {len(keep_genes):,} genes  (from {base.name})")

    merged = lc(a.merged, is_pred=False)
    print(f"merged 28-library reference: {len(merged):,} calls in that space")

    obs = {ds: lc(str(cell("mamba4", ds, ds) / "real_collapsed.txt"), is_pred=False)
           for ds in DATASETS}
    for ds in DATASETS:
        print(f"  observed[{ds}]: {len(obs[ds]):,}")
    pred = {}
    for m in MODELS:
        f = cell(m, MODEL_RNA, MODEL_RNA) / f"{ARM}_collapsed.txt"
        if f.exists():
            pred[m] = lc(str(f), is_pred=True)
            print(f"  model[{m}] {ARM}: {len(pred[m]):,}")

    # Attributes for class lookup: merged is the most complete, then the observed sets, then the
    # model. A key present in several carries the same ORF type by construction (same collapsed
    # schema), so first-wins is safe and keeps every key resolvable.
    calls = {}
    for src in [merged, *obs.values(), *pred.values()]:
        for k, v in src.items():
            calls.setdefault(k, v)

    any_obs = set().union(*(set(o) for o in obs.values()))
    rows, detail = [], {}

    # -- the CONTROL: observed calls unique to one experiment, missed by the other two --
    for ds in DATASETS:
        others = set().union(*(set(obs[o]) for o in DATASETS if o != ds))
        uniq = set(obs[ds]) - others
        r = rate(uniq, merged)
        detail[f"obs_unique_{ds}"] = dict(**r, by_class=by_class(uniq, merged, calls))
        rows.append(dict(kind="experiment-unique (control)", source=ds, **r))

    # -- the QUESTION: model calls no experiment supports --
    for m, pc in pred.items():
        uniq = set(pc) - any_obs
        r = rate(uniq, merged)
        detail[f"model_unique_{m}"] = dict(**r, by_class=by_class(uniq, merged, calls))
        rows.append(dict(kind="model-unique (question)", source=m, **r))
        # And the model calls that ARE observation-supported, as an upper anchor: these should be
        # recovered at close to 100%, and if they are not the merged reference is the problem.
        sup = set(pc) & any_obs
        rows.append(dict(kind="model, experiment-supported (anchor)", source=m, **rate(sup, merged)))

    # -- shared over-calls: both models, no experiment. A prediction shared by two architectures is
    #    either a stronger real signal or a shared systematic error; the merged reference separates
    #    those too.
    if len(pred) == 2:
        both = (set(pred[MODELS[0]]) & set(pred[MODELS[1]])) - any_obs
        r = rate(both, merged)
        detail["model_unique_both"] = dict(**r, by_class=by_class(both, merged, calls))
        rows.append(dict(kind="model-unique, BOTH models", source="+".join(MODELS), **r))

        # -- and the EXCLUSIVE remainder per model: called by this architecture, not the other, and
        #    by no experiment. `model-unique (question)` already covers each model's full
        #    unsupported set (which INCLUDES the shared ones), so this is the only stratum that was
        #    genuinely never scored: 2,056 attn-unsupported = 1,504 shared + 552 attn-exclusive.
        #    It answers a different question from the shared bar -- whether a call the OTHER
        #    architecture declined to make is worse than one both agreed on.
        for i, m in enumerate(MODELS):
            other = MODELS[1 - i]
            excl = (set(pred[m]) - set(pred[other])) - any_obs
            r = rate(excl, merged)
            detail[f"model_unique_{m}_exclusive"] = dict(**r, by_class=by_class(excl, merged, calls))
            rows.append(dict(kind="model-unique, THIS MODEL ONLY", source=m, **r))

    # -- close off the "no data there" explanation --
    # A model-unique call the merged reference misses is only evidence of a false positive if the
    # merged reference HAD signal at that locus. Split on whether merged called any ORF on the same
    # GENE: if it did, the locus is demonstrably deep enough and the specific ORF was rejected on its
    # own merits; if it did not, the gene may simply be too lowly translated for any reference and the
    # miss is uninformative. Without this split the headline is confounded with expression.
    merged_genes = {k[0] for k in merged}
    for m, pc in pred.items():
        uniq = set(pc) - any_obs
        on = {k for k in uniq if k[0] in merged_genes}
        off = uniq - on
        detail[f"model_unique_{m}_gene_covered"] = rate(on, merged)
        detail[f"model_unique_{m}_gene_uncovered"] = rate(off, merged)
        rows.append(dict(kind="model-unique, gene HAS merged calls", source=m, **rate(on, merged)))
        # The rate in this second stratum is 0 BY CONSTRUCTION -- recovery requires a merged call on
        # the gene, and this stratum is defined by their absence. Its informative quantity is n (how
        # many extras land on a gene with no detectable translation in 28 pooled libraries), never
        # the rate. Labelled so the 0% is not quoted as a measurement.
        rows.append(dict(kind="model-unique, gene has NO merged call*",
                         source=m, **rate(off, merged)))

    hdr = ["kind", "source", "n", "recovered", "rate"]
    with open(out / "overcall_recovery.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[h]:.4f}" if h == "rate" else str(r[h]) for h in hdr) + "\n")

    crows = []
    for grp, d in detail.items():
        for c, r in d.get("by_class", {}).items():   # gene-coverage strata carry no class split
            crows.append(dict(group=grp, orf_class=c, **r))
    with open(out / "overcall_recovery_by_class.tsv", "w") as fh:
        fh.write("group\torf_class\tn\trecovered\trate\n")
        for r in crows:
            fh.write(f"{r['group']}\t{r['orf_class']}\t{r['n']}\t{r['recovered']}\t{r['rate']:.4f}\n")

    (out / "overcall_meta.json").write_text(json.dumps(dict(
        merged_reference=a.merged, merged_calls_in_space=len(merged),
        n_libraries_merged=28, arm=ARM, model_rna=MODEL_RNA,
        transcript_space=dict(n_tx=len(keep_tx), n_genes=len(keep_genes), source=str(base)),
        filtering="compare_dropin_calls.build_loader (pval<=0.05, ORF>=90nt, enrichment>=0.5x)",
        key="genomic (gene_id, ORF_gstop)",
        observed_calls={ds: len(obs[ds]) for ds in DATASETS},
        model_calls={m: len(p) for m, p in pred.items()},
        rows=rows, detail=detail), indent=2))

    print(f"\n{'kind':<40}{'source':<16}{'n':>8}{'recovered':>11}{'rate':>8}")
    for r in rows:
        print(f"{r['kind']:<40}{r['source']:<16}{r['n']:>8,}{r['recovered']:>11,}{r['rate']:>8.1%}")
    print("  * rate is 0 BY CONSTRUCTION in this stratum (recovery needs a merged call on the gene,\n"
          "    and the stratum is defined by their absence). Read n, never the rate.")
    ctrl = [r for r in rows if r["kind"].startswith("experiment-unique")]
    q = [r for r in rows if r["kind"].startswith("model-unique (")]
    if ctrl and q:
        c = sum(r["recovered"] for r in ctrl) / max(1, sum(r["n"] for r in ctrl))
        print(f"\nCONTROL experiment-unique recovery {c:.1%}  vs  model-unique "
              + ", ".join(f"{r['source']} {r['rate']:.1%}" for r in q))
    print(f"\nwrote 2 TSVs + overcall_meta.json to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
