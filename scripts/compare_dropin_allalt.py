#!/usr/bin/env python3
"""All-near-cognate drop-in comparison: generalizes compare_dropin_ctg.py from ATG-vs-CTG to the full
set of RiboCode alternative starts (CTG,GTG,TTG,ACG,ATA,ATT,ATC,AAG,AGG). After the drop-in is re-run
with --alt_start_codons <all 9> (ribocode_dropin_ctg.sbatch OUTSUB=dropin_allalt), this scores
real-vs-predicted ORF calling stratified BY START CODON, so it can be read per codon (does the model
help call GTG / ACG / ... ORFs, or only ATG + CTG?).

Same genomic-locus key + filters as compare_dropin_ctg.py (min_len 90, raw pval<=0.05, pred
enrichment>=0.5); reuses that module's load_seqs / load_calls_codon / strat_report / type_breakdown and
compare_dropin_calls.load_tx2gene, so only the codon LOOP is new. For every codon present in real OR
predicted it reports strict precision/recall/F1 (both sides called that codon) plus locus recall (of real
<codon> loci, fraction the model called at all, any codon). Writes <dir>/dropin_allalt_metrics.json +
a per-codon table. cas12a env (numpy)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_dropin_ctg import (  # noqa: E402
    DEFAULT_FASTA,
    load_calls_codon,
    load_seqs,
    strat_report,
    type_breakdown,
)
from compare_dropin_calls import load_tx2gene  # noqa: E402

# ATG first (canonical anchor), then the 9 near-cognates in RiboCode start-efficiency order.
CODON_ORDER = ["ATG", "CTG", "GTG", "TTG", "ACG", "ATA", "ATT", "ATC", "AAG", "AGG"]


def codon_counts(calls):
    cc = {}
    for v in calls.values():
        cc[v["codon"]] = cc.get(v["codon"], 0) + 1
    return dict(sorted(cc.items(), key=lambda kv: -kv[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dropin_dir", required=True, help="dropin_allalt dir (all-near-cognate drop-in)")
    ap.add_argument("--profiles", required=True, help="pred_profiles.npz (for the enrichment filter)")
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA),
                    help="RiboCode annotation transcripts_sequence.fa for start-codon lookup")
    ap.add_argument("--tx2gene", default=None, help="tx_id->gene_id map (default = compare's human)")
    ap.add_argument("--max_pval", type=float, default=0.05)
    ap.add_argument("--min_len", type=int, default=90)
    ap.add_argument("--min_enrichment", type=float, default=0.5)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    d = Path(args.dropin_dir)
    out = Path(args.out) if args.out else d
    out.mkdir(parents=True, exist_ok=True)

    npz = np.load(args.profiles, allow_pickle=False)
    keep_tx = set(npz["tx_ids"].tolist())
    lengths = npz["lengths"]
    pflat = npz["pred_flat"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    idx = {str(t): i for i, t in enumerate(npz["tx_ids"])}

    def enrich(tx, a, e):
        i = idx.get(str(tx))
        if i is None:
            return float("nan")
        p = pflat[off[i]:off[i + 1]]
        if a < 0 or e > len(p) or e <= a:
            return float("nan")
        return float(p[a:e].mean() * len(p))

    t2g = load_tx2gene(args.tx2gene) if args.tx2gene else load_tx2gene()
    keep_genes = {t2g[t] for t in keep_tx if t in t2g}
    print(f"test tx: {len(keep_tx):,}; test genes: {len(keep_genes):,}", file=sys.stderr)
    if keep_tx and not keep_genes:
        raise SystemExit("ERROR: 0 test tx matched tx2gene map (wrong assembly?)")

    files = {v: d / f"{v}_collapsed.txt" for v in ("real", "pred_obsdepth", "pred_preddepth")}
    for name, p in files.items():
        if not p.exists():
            raise SystemExit(f"ERROR: missing {p} -- has the all-alt drop-in finished?")

    need = set()
    for p in files.values():
        with open(p) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            ti = h.index("transcript_id")
            for ln in fh:
                need.add(ln.rstrip("\n").split("\t")[ti])
    print(f"resolving start codons for {len(need):,} transcripts ...", file=sys.stderr)
    seqs = load_seqs(args.fasta, need)

    real = load_calls_codon(files["real"], keep_genes, seqs, args.max_pval, args.min_len)
    pred_obs = load_calls_codon(files["pred_obsdepth"], keep_genes, seqs, args.max_pval, args.min_len,
                                enrich, args.min_enrichment)
    pred_pred = load_calls_codon(files["pred_preddepth"], keep_genes, seqs, args.max_pval,
                                 args.min_len, enrich, args.min_enrichment)

    present = [c for c in CODON_ORDER
               if any(c in codon_counts(x) for x in (real, pred_obs, pred_pred))]

    res = {"dropin_dir": str(d), "min_len": args.min_len, "max_pval": args.max_pval,
           "min_enrichment": args.min_enrichment, "codons": present,
           "codon_counts": {"real": codon_counts(real), "pred_obsdepth": codon_counts(pred_obs),
                            "pred_preddepth": codon_counts(pred_pred)}}
    for codon in present:
        res[f"pred_obsdepth_vs_real_{codon}"] = strat_report(pred_obs, real, codon,
                                                             f"pred_obsdepth_vs_real_{codon}")
        res[f"pred_preddepth_vs_real_{codon}"] = strat_report(pred_pred, real, codon,
                                                              f"pred_preddepth_vs_real_{codon}")
        res[f"real_{codon}_type_breakdown"] = type_breakdown(real, codon)
    (out / "dropin_allalt_metrics.json").write_text(json.dumps(res, indent=2))

    print("\n=== all-near-cognate drop-in: per-start-codon ORF calling (held-out) ===")
    print("start-codon composition (genomic loci, after filters):")
    for name, c in (("real", real), ("pred_obsdepth", pred_obs), ("pred_preddepth", pred_pred)):
        parts = ", ".join(f"{k}:{v}" for k, v in codon_counts(c).items())
        print(f"  {name:<15} n={len(c):<7,} {parts}")
    print(f"\n{'codon':>6} {'nRealC':>7} {'nPredC':>7} {'precS':>7} {'recallS':>8} "
          f"{'F1strict':>9} {'locusRec':>9}  (pred_obsdepth vs real)")
    for codon in present:
        r = res[f"pred_obsdepth_vs_real_{codon}"]; s = r["strict"]
        print(f"{codon:>6} {r['n_real_codon']:>7} {r['n_pred_codon']:>7} {s['precision']:>7.3f} "
              f"{s['recall']:>8.3f} {s['f1']:>9.3f} {r['locus_recall_of_real_codon']:>9.3f}")
    print(f"\nwrote {out / 'dropin_allalt_metrics.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
