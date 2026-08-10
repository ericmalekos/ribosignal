#!/usr/bin/env python3
"""Poly(A) gate for RNA-seq arms, read off salmon quants.

The failure this exists to catch: a Ribo-Zero TOTAL-RNA library keeps tRNA- and 7SL-derived
transcripts, many of which GENCODE types as lncRNA, so a biotype filter cannot remove them. They
soak up the TPM mass and the universe built from that quant collapses (GSE243134 liver: 3,321
transcripts). SRA `library_selection` does not distinguish this -- it reads `cDNA` either way.

What is measured, per sample:
  * TPM mass by transcript biotype class.
  * Concentration: what fraction of all TPM sits in the top 20 transcripts. A poly(A) mRNA library
    spreads its mass; a contaminated one piles it into a handful of repeat-derived loci.
  * The actual top transcripts by TPM, with gene names, so a 7SL/tRNA pile-up is visible by name
    rather than inferred.
  * salmon mapping rate against the pc+lncRNA decoy-aware index.

This gate judges whether a library may DEFINE A UNIVERSE. It does not judge whether the library is
usable as model coverage input -- salmon EM plus length normalisation and STAR mm1 plus ncRNA
filtering disagree strongly on multi-copy loci, and a TPM-contaminated library can still produce
clean pack coverage.
"""
import argparse
import json
import pathlib
import sys
from collections import defaultdict

# GENCODE types that a poly(A) selection should have largely removed. tRNA/7SL-derived loci are
# frequently typed lncRNA or misc_RNA, so the name-level table below matters as much as this set.
SUSPECT = {"rRNA", "Mt_rRNA", "rRNA_pseudogene", "misc_RNA", "snRNA", "snoRNA", "scaRNA",
           "vault_RNA", "sRNA", "scRNA", "Mt_tRNA", "ribozyme"}
SUSPECT_NAME_PREFIX = ("RN7SL", "RN7SK", "RNU", "RMRP", "RPPH1", "SNOR", "SCARNA", "TRNA", "RNY")


def load_biotype(path):
    m = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ti = hdr.index("tx_id")
        bi = hdr.index("transcript_type")
        gi = hdr.index("gene_name")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ti, bi, gi):
                m[f[ti]] = (f[bi], f[gi])
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salmon-dir", required=True)
    ap.add_argument("--tx2biotype", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", help="dataset label per sample (same order)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()

    bt = load_biotype(a.tx2biotype)
    labels = a.labels or [""] * len(a.samples)
    rows, lines = [], []

    for srr, lab in zip(a.samples, labels):
        q = pathlib.Path(a.salmon_dir) / srr / "quant.sf"
        if not q.exists():
            print(f"missing {q}", file=sys.stderr)
            continue
        by_class = defaultdict(float)
        tot = 0.0
        entries = []
        with open(q) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                tpm = float(f[3])
                if tpm <= 0:
                    continue
                tx = f[0].split("|")[0]
                b, g = bt.get(tx, ("unknown", f[0].split("|")[5] if "|" in f[0] else tx))
                by_class[b] += tpm
                tot += tpm
                entries.append((tpm, g, b))
        entries.sort(reverse=True)
        topn = entries[:a.top]
        top_frac = sum(e[0] for e in topn) / tot if tot else 0.0
        susp = sum(v for k, v in by_class.items() if k in SUSPECT)
        susp_name = sum(t for t, g, b in entries if g.upper().startswith(SUSPECT_NAME_PREFIX))
        mi = pathlib.Path(a.salmon_dir) / srr / "aux_info" / "meta_info.json"
        rate = json.load(open(mi))["percent_mapped"] if mi.exists() else float("nan")

        rows.append(dict(sample=srr, dataset=lab, mapping_rate=rate,
                         pct_protein_coding=100 * by_class.get("protein_coding", 0) / tot,
                         pct_lncRNA=100 * by_class.get("lncRNA", 0) / tot,
                         pct_suspect_biotype=100 * susp / tot,
                         pct_suspect_by_name=100 * susp_name / tot,
                         pct_top20=100 * top_frac))
        lines.append(f"\n## {srr}  ({lab})")
        lines.append(f"  salmon mapping rate      {rate:.1f}%")
        lines.append(f"  protein_coding TPM       {100*by_class.get('protein_coding',0)/tot:.1f}%")
        lines.append(f"  lncRNA TPM               {100*by_class.get('lncRNA',0)/tot:.1f}%")
        lines.append(f"  suspect-biotype TPM      {100*susp/tot:.1f}%")
        lines.append(f"  suspect-by-name TPM      {100*susp_name/tot:.1f}%   "
                     f"(RN7SL/RN7SK/RNU/RMRP/RPPH1/SNOR/...)")
        lines.append(f"  top-{a.top} concentration    {100*top_frac:.1f}% of all TPM")
        lines.append(f"  top transcripts: " + ", ".join(f"{g}({b},{t:.0f})" for t, g, b in topn[:8]))

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    hdr = list(rows[0].keys())
    with open(out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.3f}" if isinstance(r[k], float) else str(r[k])
                               for k in hdr) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    sys.exit(main())
