#!/usr/bin/env python3
"""Start-codon-stratified RiboCode drop-in comparison (ATG vs CTG).

After the drop-in is re-run with --alt_start_codons CTG (ribocode_dropin_ctg.sbatch -> dropin_ctg/),
this splits the real / pred_obsdepth / pred_preddepth ORF calls by their ACTUAL start codon (the 3 nt
at ORF_tstart, looked up in the RiboCode annotation FASTA) and reports predicted-vs-real ORF-calling
precision/recall/F1 SEPARATELY for ATG-initiated and CTG-initiated ORFs. It uses the SAME genomic-locus
key (gene_id, ORF_gstop) and the SAME filters (min_len, raw pval_combined, pred enrichment) as
compare_dropin_calls.py, so the ATG rows should reproduce the ATG-only drop-in and the CTG rows are the
new result: does the predicted profile recover / propose non-AUG ORFs, or only canonical ATG ones?

For CTG two readouts are given because the collapse representative isoform (hence the recorded codon)
can differ between the real and predicted call at the same genomic stop:
  * strict   prf(pred CTG-calls, real CTG-calls)         -- both sides called it CTG
  * locus    of real CTG loci, fraction pred called at all (any codon) = recall; and of pred CTG loci,
             fraction real called at all = precision
cas12a env (numpy). Writes <dropin_ctg>/dropin_ctg_metrics.json + a printed table.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
from compare_dropin_calls import NONCANON, load_tx2gene, prf  # noqa: E402

ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/"
           "expression_context_human")
DEFAULT_FASTA = ECH / "data/ribocode_annot/transcripts_sequence.fa"


def load_seqs(fasta, need):
    """Stream the annotation FASTA once, keep only sequences for transcripts in `need`."""
    seq = {}
    cur = None
    keep = False
    buf = []
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None and keep:
                    seq[cur] = "".join(buf)
                cur = ln[1:].split()[0]
                keep = cur in need
                buf = []
            elif keep:
                buf.append(ln.strip())
        if cur is not None and keep:
            seq[cur] = "".join(buf)
    return seq


def start_codon(seq, tx, tstart):
    """3-nt start codon at 1-based ORF_tstart, uppercased, U->T. None if unresolved."""
    s = seq.get(tx)
    if s is None or tstart < 1 or tstart + 2 > len(s):
        return None
    return s[tstart - 1:tstart + 2].upper().replace("U", "T")


def load_calls_codon(path, keep_genes, seqs, max_pval, min_len, enrich=None, min_enrichment=0.0):
    """Parse a *_collapsed.txt into {genomic_key: {codon,type,tx,tstart}} with the same filters
    compare_dropin_calls.load_calls applies (genomic key = (gene_id, ORF_gstop), restricted to
    keep_genes; raw pval_combined <= max_pval; ORF_length >= min_len; pred enrichment >= threshold)."""
    calls = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tx = f[ci["transcript_id"]]
            gene = f[ci["gene_id"]] if "gene_id" in ci else "NA"
            if keep_genes is not None and gene not in keep_genes:
                continue
            try:
                tstart = int(f[ci["ORF_tstart"]])
                length = int(f[ci["ORF_length"]])
            except (ValueError, KeyError):
                continue
            if length < min_len:
                continue
            try:
                pc = float(f[ci["pval_combined"]])
            except (ValueError, KeyError):
                pc = float("nan")
            if max_pval is not None and not (pc <= max_pval):
                continue
            if enrich is not None and min_enrichment > 0:
                try:
                    enr = enrich(tx, tstart - 1, int(f[ci["ORF_tstop"]]))
                except (ValueError, KeyError):
                    enr = float("nan")
                if not (enr >= min_enrichment):
                    continue
            key = (gene, f[ci["ORF_gstop"]])
            calls[key] = {"codon": start_codon(seqs, tx, tstart),
                          "type": f[ci["ORF_type"]] if "ORF_type" in ci else "NA",
                          "tx": tx, "tstart": tstart}
    return calls


def subset(calls, codon):
    return {k for k, v in calls.items() if v["codon"] == codon}


def strat_report(pred, real, codon, label):
    """prf on the codon-subset of both sides (strict), plus locus-level recall/precision where the
    other side may have called the same genomic locus with a different codon."""
    pred_c = subset(pred, codon)
    real_c = subset(real, codon)
    pred_all = set(pred)
    real_all = set(real)
    strict = prf(pred_c, real_c)
    locus_recall = (len(real_c & pred_all) / len(real_c)) if real_c else float("nan")
    locus_prec = (len(pred_c & real_all) / len(pred_c)) if pred_c else float("nan")
    return {"label": label, "codon": codon, "strict": strict,
            "n_pred_codon": len(pred_c), "n_real_codon": len(real_c),
            "locus_recall_of_real_codon": locus_recall,
            "locus_precision_of_pred_codon": locus_prec}


def type_breakdown(calls, codon):
    out = {}
    for k, v in calls.items():
        if v["codon"] == codon:
            out[v["type"]] = out.get(v["type"], 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dropin_dir", required=True, help="dropin_ctg dir (ATG+CTG drop-in output)")
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

    files = {"real": d / "real_collapsed.txt",
             "pred_obsdepth": d / "pred_obsdepth_collapsed.txt",
             "pred_preddepth": d / "pred_preddepth_collapsed.txt"}
    for name, p in files.items():
        if not p.exists():
            raise SystemExit(f"ERROR: missing {p} -- has ribocode_dropin_ctg.sbatch finished?")

    # collect all transcripts appearing across the three files (pre-filter) for FASTA streaming
    need = set()
    for p in files.values():
        with open(p) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            ti = h.index("transcript_id")
            for ln in fh:
                need.add(ln.rstrip("\n").split("\t")[ti])
    print(f"resolving start codons for {len(need):,} transcripts from {args.fasta} ...",
          file=sys.stderr)
    seqs = load_seqs(args.fasta, need)

    real = load_calls_codon(files["real"], keep_genes, seqs, args.max_pval, args.min_len)
    pred_obs = load_calls_codon(files["pred_obsdepth"], keep_genes, seqs, args.max_pval, args.min_len,
                                enrich, args.min_enrichment)
    pred_pred = load_calls_codon(files["pred_preddepth"], keep_genes, seqs, args.max_pval,
                                 args.min_len, enrich, args.min_enrichment)

    def codon_counts(c):
        cc = {}
        for v in c.values():
            cc[v["codon"]] = cc.get(v["codon"], 0) + 1
        return dict(sorted(cc.items(), key=lambda kv: -kv[1]))

    res = {"dropin_dir": str(d), "min_len": args.min_len, "max_pval": args.max_pval,
           "min_enrichment": args.min_enrichment,
           "codon_counts": {"real": codon_counts(real), "pred_obsdepth": codon_counts(pred_obs),
                            "pred_preddepth": codon_counts(pred_pred)},
           "real_CTG_type_breakdown": type_breakdown(real, "CTG"),
           "real_ATG_type_breakdown": type_breakdown(real, "ATG")}
    for codon in ("ATG", "CTG"):
        res[f"pred_obsdepth_vs_real_{codon}"] = strat_report(pred_obs, real, codon,
                                                             f"pred_obsdepth_vs_real_{codon}")
        res[f"pred_preddepth_vs_real_{codon}"] = strat_report(pred_pred, real, codon,
                                                              f"pred_preddepth_vs_real_{codon}")
    (out / "dropin_ctg_metrics.json").write_text(json.dumps(res, indent=2))

    # ---- printed report ----
    print("\n=== ATG+CTG drop-in: start-codon-stratified ORF calling (held-out Hepatocytes) ===")
    print("start-codon composition of each call set (genomic loci, after filters):")
    for name, c in (("real", real), ("pred_obsdepth", pred_obs), ("pred_preddepth", pred_pred)):
        tot = len(c)
        parts = ", ".join(f"{k}:{v} ({100*v/tot:.1f}%)" for k, v in codon_counts(c).items())
        print(f"  {name:<15} n={tot:<7,} {parts}")
    print("\nreal CTG-initiated ORFs by type:", res["real_CTG_type_breakdown"])
    print("\n%-32s %7s %7s %7s %7s %7s %8s %9s" %
          ("comparison", "nPredC", "nRealC", "match", "precS", "recallS", "F1strict", "locusRec"))
    for codon in ("ATG", "CTG"):
        for base in ("pred_obsdepth", "pred_preddepth"):
            r = res[f"{base}_vs_real_{codon}"]
            s = r["strict"]
            print("%-32s %7d %7d %7d %7.3f %7.3f %8.3f %9.3f" %
                  (f"{base}_vs_real [{codon}]", r["n_pred_codon"], r["n_real_codon"], s["n_match"],
                   s["precision"], s["recall"], s["f1"], r["locus_recall_of_real_codon"]))
    print(f"\nwrote {out / 'dropin_ctg_metrics.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
