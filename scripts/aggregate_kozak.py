#!/usr/bin/env python3
"""Consolidate the 4 Kozak arms (V0 heuristic / V1 no-Kozak / V2 empirical PWM / V3 learned) into one
comparison table that answers the 3 questions from KOZAK_PLAN.md:

  Q1  How do metrics hold if the hand-picked Kozak heuristic is REMOVED?    -> V1 vs V0
  Q2  Does an EMPIRICALLY-fit Kozak (PWM from annotated starts) beat it?     -> V2 vs V0
  Q3  Can the model LEARN the context weights itself?                        -> V3 vs V0/V1
                                                                                + learned-kernel r

Kozak only modulates the START channel, so the metrics that should move are non-AUG / uORF detection
(the CTG drop-in and the uORF/non-canonical localization AUROC), NOT bulk CDS prediction. Bulk metrics
(annotated-ORF AUROC, whole-pc profile Pearson) are reported as a SANITY check -- they should be flat
across arms; if they move, the arm broke something unrelated.

Reads per arm (all under results/kozak/<arm>_onehot_fib2hep/):
  extra_metrics.json                 (profile Pearson sanity; written by train_kozak_arms.sbatch)
  localization_metrics.json          (uORF / non-canonical / annotated AUROC + fidelity)
  dropin_ctg/dropin_ctg_metrics.json (ATG vs CTG strict precision/recall/F1 + locus recall)
plus figures/kozak/learned_vs_empirical_vs_heuristic.json (V3 learned-vs-empirical kernel Pearson r).

Missing JSONs (an arm not yet evaluated) render as "--" so this can be run mid-chain. Writes
results/kozak/kozak_summary.md + .json. cas12a env (stdlib only)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ARMS = [("v0_heuristic", "V0 heuristic"), ("v1_nokozak", "V1 no-Kozak"),
        ("v2_pwm", "V2 empirical"), ("v3_learned", "V3 learned")]
# The CTG drop-in reports two predicted-depth variants; obsdepth = predicted SHAPE at the OBSERVED
# total depth, the fairest "does the profile recover the ORF in the caller" test (depth held fixed).
DEPTH = "pred_obsdepth"


def load(p):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else None


def g(d, *keys, default=None):
    """Nested .get; returns default if any level is missing/None."""
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return default if d is None else d


def fmt(v, nd=3):
    if v is None:
        return "--"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kozak_dir", default=str(NEW / "results/kozak"))
    ap.add_argument("--kernel_json",
                    default=str(NEW / "figures/kozak/learned_vs_empirical_vs_heuristic.json"))
    ap.add_argument("--out_md", default=str(NEW / "results/kozak/kozak_summary.md"))
    ap.add_argument("--out_json", default=str(NEW / "results/kozak/kozak_summary.json"))
    args = ap.parse_args()

    kd = Path(args.kozak_dir)
    per = {}
    for slug, _ in ARMS:
        run = kd / f"{slug}_onehot_fib2hep"
        per[slug] = {
            "extra": load(run / "extra_metrics.json"),
            "loc": load(run / "localization_metrics.json"),
            "ctg": load(run / "dropin_ctg" / "dropin_ctg_metrics.json"),
        }
    kernel = load(args.kernel_json)

    def ctg(slug, codon, field):
        block = g(per[slug]["ctg"], f"{DEPTH}_vs_real_{codon}")
        if field in ("precision", "recall", "f1"):
            return g(block, "strict", field)
        return g(block, field)

    def loc_auroc(slug, group, score="auroc_pred_frame0_lengthctrl"):
        return g(per[slug]["loc"], "discrimination", group, score)

    def loc_n(slug, group):
        return g(per[slug]["loc"], "discrimination", group, "n_pos")

    def loc_fid(slug, group, field):
        return g(per[slug]["loc"], "fidelity", group, field)

    def extra(slug, block, field="profile_pearson_median"):
        return g(per[slug]["extra"], block, field)

    # (row label, extractor) grouped into sections. Length-controlled AUROC is the headline
    # discrimination number (removes the ORF-length confound); raw pred_frame0 kept alongside.
    sections = [
        ("CTG non-AUG drop-in (predicted profile @ observed depth) -- THE sharpest Kozak test", [
            ("CTG strict F1", lambda s: ctg(s, "CTG", "f1")),
            ("CTG strict recall", lambda s: ctg(s, "CTG", "recall")),
            ("CTG strict precision", lambda s: ctg(s, "CTG", "precision")),
            ("CTG locus recall (real CTG loci pred-called, any codon)",
             lambda s: ctg(s, "CTG", "locus_recall_of_real_codon")),
            ("n real CTG loci", lambda s: ctg(s, "CTG", "n_real_codon")),
            ("n pred CTG loci", lambda s: ctg(s, "CTG", "n_pred_codon")),
        ]),
        ("uORF / non-canonical localization AUROC (length-controlled pred_frame0)", [
            ("non-canonical AUROC", lambda s: loc_auroc(s, "noncanonical")),
            ("  (n non-canonical positives)", lambda s: loc_n(s, "noncanonical")),
            ("uORF AUROC", lambda s: loc_auroc(s, "uORF")),
            ("  (n uORF positives)", lambda s: loc_n(s, "uORF")),
            ("non-canonical density Pearson (fidelity, +)",
             lambda s: loc_fid(s, "noncanonical", "density_pearson")),
            ("uORF frame0 Pearson (fidelity, +)",
             lambda s: loc_fid(s, "uORF", "frame0_pearson")),
        ]),
        ("SANITY -- should be ~flat across arms (Kozak only touches the start channel)", [
            ("ATG strict F1 (drop-in)", lambda s: ctg(s, "ATG", "f1")),
            ("annotated-ORF AUROC (length-controlled)", lambda s: loc_auroc(s, "annotated")),
            ("all-ORF AUROC (length-controlled)", lambda s: loc_auroc(s, "all")),
            ("pc whole-tx profile Pearson (median)", lambda s: extra(s, "protein_coding")),
            ("pc 5'UTR profile Pearson (median)", lambda s: extra(s, "uorf_5utr")),
        ]),
    ]

    slugs = [s for s, _ in ARMS]
    names = [n for _, n in ARMS]
    lines = []
    lines.append("## Kozak start-context ablation -- Fibroblast -> Hepatocytes (one-hot, orf_v2_attn)\n")
    lines.append("Arms differ ONLY in the ORF-track start-propensity channel (channel 3). "
                 "Occupancy + stop channels are identical across all four.\n")
    lines.append("| metric | " + " | ".join(names) + " |")
    lines.append("|" + "---|" * (len(names) + 1))
    json_rows = {}
    for title, rows in sections:
        lines.append(f"| **{title}** | " + " | ".join([""] * len(names)) + " |")
        for label, fn in rows:
            vals = [fn(s) for s in slugs]
            nd = 0 if label.strip().startswith("(") or label.startswith("n ") else 3
            lines.append(f"| {label} | " + " | ".join(fmt(v, nd) for v in vals) + " |")
            json_rows[label.strip()] = dict(zip(slugs, vals))

    # Q3 mechanistic: the learned kernel vs the empirical PWM.
    r = g(kernel, "learned_vs_empirical_pearson_r")
    lines.append("")
    lines.append("### Q3 mechanistic: what did V3 learn?\n")
    if r is not None:
        lines.append(f"- Learned 4x10 start-context kernel vs empirical PWM: **Pearson r = {r:.3f}** "
                     f"over the 7 context positions (per-position centered). See "
                     f"`figures/kozak/learned_vs_empirical_vs_heuristic.png`.")
    else:
        lines.append("- Learned-kernel comparison pending (V3 not yet plotted).")

    # Auto-interpretation of the 3 questions from the numbers actually present.
    def d(a, b):
        return (a - b) if (a is not None and b is not None) else None

    ctg_f1 = {s: ctg(s, "CTG", "f1") for s in slugs}
    uorf_au = {s: loc_auroc(s, "uORF") for s in slugs}
    lines.append("\n### Answers\n")
    q1 = d(ctg_f1["v1_nokozak"], ctg_f1["v0_heuristic"])
    q2 = d(ctg_f1["v2_pwm"], ctg_f1["v0_heuristic"])
    q3 = d(ctg_f1["v3_learned"], ctg_f1["v0_heuristic"])
    lines.append(f"- **Q1 (remove heuristic):** V1-V0 CTG-F1 delta = {fmt(q1)}; "
                 f"uORF AUROC {fmt(uorf_au['v1_nokozak'])} vs {fmt(uorf_au['v0_heuristic'])}. "
                 f"(>=0 => the heuristic was not helping.)")
    lines.append(f"- **Q2 (empirical PWM):** V2-V0 CTG-F1 delta = {fmt(q2)}. "
                 f"(<=0 corroborates the pre-registered prediction that a CDS-fit PWM is a "
                 f"miscalibrated prior for uORF/CTG starts.)")
    lines.append(f"- **Q3 (learned):** V3-V0 CTG-F1 delta = {fmt(q3)}; "
                 f"learned-vs-empirical kernel r = {fmt(r)}. "
                 f"(The one-hot backbone can reconstruct Kozak from sequence; r near the V0/V2 "
                 f"agreement means it rediscovered the empirical matrix without supervision.)")

    md = "\n".join(lines) + "\n"
    Path(args.out_md).write_text(md)
    Path(args.out_json).write_text(json.dumps(
        {"arms": dict(ARMS), "depth_variant": DEPTH, "rows": json_rows,
         "learned_vs_empirical_pearson_r": r,
         "deltas_ctg_f1_vs_v0": {"v1_nokozak": q1, "v2_pwm": q2, "v3_learned": q3}}, indent=2))
    print(md)
    print(f"wrote {args.out_md}\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
