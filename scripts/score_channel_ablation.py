#!/usr/bin/env python3
"""Score the ORF-track channel ablation: what does each input channel actually buy?

Each arm re-runs prediction with ONE channel of the ORF track zeroed (ch0/ch1/ch2 = frame-0/1/2
occupancy, ch3 = start_ext propensity, ch4 = is_stop), plus `f012` which zeroes all three frame
channels together. The readout is the change in ORF-calling F1 against the observed reference,
per ORF class -- a channel that matters is one whose removal costs calls.

Signs are the point, not magnitudes: a POSITIVE delta means the model calls that class BETTER
without the channel, i.e. the channel is actively misleading the model for that class.

Reference and base come from the same matched cell of the mouse-liver 3x3 (Ribo and RNA both
GSE243134), so the only difference between base and arm is the zeroed channel.

Filtering is imported from `compare_dropin_calls.build_loader` -- the same code path behind every
other drop-in number in this project -- and each arm's predicted calls are filtered with that ARM's
own pred_flat, since the enrichment filter must see the predictions it is judging. Previously this
scoring existed only as an inline shell heredoc, which is exactly what docs/PIPELINE_POLICY.md
forbids: a number in the manuscript needs a script someone can rerun.

Keying is GENOMIC (gene_id, ORF_gstop), matching score_liver3x3.py.
"""
import argparse
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_dropin_calls import build_loader  # noqa: E402
from score_liver3x3 import CLASSES, prf  # noqa: E402

ARMS = ["ch0", "ch1", "ch2", "ch3", "ch4", "f012"]
MODELS = ["attn", "mamba4"]
VARIANTS = ["pred_obsdepth", "pred_preddepth"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ablation-dir", required=True,
                    help="dir holding <model>_<arm>/ subdirs from orf_channel_ablation.sbatch")
    ap.add_argument("--base-dir", required=True,
                    help="matched-cell dir with the FULL-track dump (<model>_ribo-X_rna-X)")
    ap.add_argument("--base-pattern", default="{model}_ribo-gse243134_rna-gse243134")
    # build_loader defaults to the HUMAN tx2biotype; this ablation is mouse. Default to the mouse
    # map here rather than requiring the flag: build_loader errors loudly on a species mismatch, but
    # only after the expensive load, and the wrong default is the easier mistake to make.
    ap.add_argument("--tx2gene",
                    default=str(pathlib.Path(__file__).resolve().parent.parent
                                / "data" / "tx2biotype_mouse.tsv"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    abl = pathlib.Path(a.ablation_dir)
    base_root = pathlib.Path(a.base_dir)
    rows = []

    for model in MODELS:
        bdir = base_root / a.base_pattern.format(model=model)
        bprof = bdir / "pred_profiles.npz"
        if not bprof.exists():
            print(f"  SKIP {model}: no base dump at {bprof}", file=sys.stderr)
            continue
        blc, bkeep, _ = build_loader(str(bprof), key="genomic", tx2gene=a.tx2gene)
        # The reference is the OBSERVED call set of the same cell, loaded under the base dump's
        # filtering so base and every arm are scored against one identical reference.
        ref = blc(str(bdir / "real_collapsed.txt"), is_pred=False)

        for variant in VARIANTS:
            base_calls = blc(str(bdir / f"{variant}_collapsed.txt"), is_pred=True)
            base_sc = {c: prf(base_calls, ref, subset=[c]) for c in CLASSES}
            base_sc["ALL"] = prf(base_calls, ref)

            for arm in ARMS:
                adir = abl / f"{model}_{arm}"
                aprof = adir / "pred_profiles.npz"
                acall = adir / f"{variant}_collapsed.txt"
                if not (aprof.exists() and acall.exists()):
                    print(f"  SKIP {model}/{arm}/{variant}: missing dump or calls", file=sys.stderr)
                    continue
                alc, akeep, _ = build_loader(str(aprof), key="genomic", tx2gene=a.tx2gene)
                # An arm that dumped a different transcript set is not comparable to the base
                # (feedback_dump_silent_skip_onehot_fasta): the delta would mix a channel effect
                # with a change of denominator.
                if akeep != bkeep:
                    print(f"  WARN {model}/{arm}: test tx differ from base "
                          f"({len(akeep):,} vs {len(bkeep):,}) -- delta not comparable",
                          file=sys.stderr)
                arm_calls = alc(str(acall), is_pred=True)
                for c in list(CLASSES) + ["ALL"]:
                    sc = prf(arm_calls, ref, subset=None if c == "ALL" else [c])
                    b = base_sc[c]
                    rows.append(dict(
                        model=model, arm=arm, dropin=variant, orf_class=c,
                        f1=f"{sc['f1']:.4f}", base_f1=f"{b['f1']:.4f}",
                        delta=f"{sc['f1'] - b['f1']:+.4f}",
                        precision=f"{sc['precision']:.4f}", recall=f"{sc['recall']:.4f}",
                        n_pred=sc["n_pred"], n_ref=sc["n_ref"]))

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
