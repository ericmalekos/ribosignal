#!/usr/bin/env python3
"""O3 Phase 1: extract the FP-filter feature vector per candidate ORF, from the model's SMOOTH
predicted per-nt density (where the uniformity/shape signal lives -- NOT the Poisson sample),
plus a real/spurious label.

Features (all numpy on the predicted profile; no BAM):
  RibORF-derived: frac_f0 (in-frame fraction), f1max (fraction of codons where frame-0 dominates),
                  pme (percent-of-max-entropy of per-codon in-frame signal = uniformity, 1=uniform)
  over-smoothing: gini, cv, max_to_mean, frac_zero (of per-codon in-frame values; a real ORF is
                  spiky/ramped, a spurious smooth ORF flattens), ramp5 (5' start enrichment),
                  drop3 (3'-of-stop signal ratio)
  context:        len_codons
Label: 1 if the ORF (coord key) is in the deep observed calls (--truth), else 0.

Usage: orf_features.py --profiles pred_profiles.npz --candidates <collapsed.txt>
                       --truth <real_collapsed.txt> --out features.tsv
"""
import argparse

import numpy as np


def coord_key(orf_id):
    return "_".join(orf_id.rsplit("_", 3)[-3:])   # gstart_gstop_len (annotation-version robust)


def read_collapsed(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(f)
    return rows, ci


def features(prof):
    """prof = 1-D predicted density over the ORF (start codon .. stop, in-frame at position 0)."""
    L = (len(prof) // 3) * 3
    if L < 9:
        return None
    cod = prof[:L].reshape(-1, 3).astype(np.float64)   # (n_cod, 3): col0 = in-frame
    n = cod.shape[0]
    fsum = cod.sum(0)                                   # [f0,f1,f2]
    tot = fsum.sum() + 1e-9
    c = cod[:, 0]                                       # per-codon in-frame signal
    csum = c.sum() + 1e-9
    p = c / csum
    nz = p[p > 0]
    pme = float(-(nz * np.log(nz)).sum() / np.log(n)) if n > 1 else 0.0   # 1=uniform
    f1max = float((cod[:, 0] >= cod.max(1)).mean())    # codons where in-frame is (tied-)max
    # Gini of c
    cs = np.sort(c)
    gini = float((2 * np.arange(1, n + 1) - n - 1).dot(cs) / (n * cs.sum() + 1e-9))
    mean_c = c.mean() + 1e-9
    return {
        "frac_f0": float(fsum[0] / tot),
        "f1max": f1max,
        "pme": pme,
        "gini": gini,
        "cv": float(c.std() / mean_c),
        "max_to_mean": float(c.max() / mean_c),
        "frac_zero": float((c == 0).mean()),
        "ramp5": float(c[:5].sum() / csum),
        "drop3": float(c[-3:].mean() / mean_c),         # <1 = ribosome release drop at stop
        "len_codons": n,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    npz = np.load(a.profiles, allow_pickle=True)
    tx = np.array([str(x) for x in npz["tx_ids"]])
    lengths = npz["lengths"].astype(np.int64)
    off = np.concatenate([[0], np.cumsum(lengths)])
    pred = npz["pred_flat"]
    row_of = {t: i for i, t in enumerate(tx)}

    truth_rows, tci = read_collapsed(a.truth)
    truth = {coord_key(r[tci["ORF_ID"]]) for r in truth_rows}

    cand, ci = read_collapsed(a.candidates)
    feat_names = ["frac_f0", "f1max", "pme", "gini", "cv", "max_to_mean", "frac_zero",
                  "ramp5", "drop3", "len_codons"]
    n_out = n_skip = 0
    with open(a.out, "w") as o:
        o.write("ORF_ID\tORF_type\t" + "\t".join(feat_names) + "\tlabel\n")
        for r in cand:
            txid = r[ci["transcript_id"]]
            i = row_of.get(txid)
            if i is None:
                n_skip += 1
                continue
            ts = int(r[ci["ORF_tstart"]])
            olen = int(r[ci["ORF_length"]])
            s = off[i] + (ts - 1)
            prof = np.asarray(pred[s:s + olen], dtype=np.float64)
            if prof.size < 9:
                n_skip += 1
                continue
            fv = features(prof)
            if fv is None:
                n_skip += 1
                continue
            lab = 1 if coord_key(r[ci["ORF_ID"]]) in truth else 0
            o.write(f"{r[ci['ORF_ID']]}\t{r[ci['ORF_type']]}\t"
                    + "\t".join(f"{fv[k]:.5g}" for k in feat_names) + f"\t{lab}\n")
            n_out += 1
    print(f"wrote {a.out}: {n_out} ORFs featurized, {n_skip} skipped (tx missing / too short)")


if __name__ == "__main__":
    main()
