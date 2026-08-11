#!/usr/bin/env python3
"""Per-codon ribosome occupancy from OBSERVED and PREDICTED profiles, on the same transcripts.

Every `pred_profiles.npz` carries `obs_flat` and `pred_flat` over an identical transcript axis, so the
empirical codon statistic and the model's are the SAME computation applied to two arrays. That is the
whole design: no separate pipeline, no alignment step, nothing that could differ between the two
sides except the numbers themselves.

METHOD
  For each transcript with an annotated CDS and at least --min-cds-counts in the CDS:
    * take the per-nt profile over the CDS only,
    * divide by that transcript's CDS mean  -> removes expression, so a highly translated transcript
      does not dominate; every transcript contributes a relative profile,
    * walk the CDS in codons and accumulate the normalised value onto a codon identity.
  Averaging over all transcripts gives one value per sense codon: >1 means the ribosome is enriched
  there relative to that transcript's own average.

P-SITE vs A-SITE, reported separately and never silently
  RiboCode assigns each read to its P-site. The conventional dwell proxy is the A-site codon, which
  is the NEXT codon along. Both are emitted:
    psite[C] = mean normalised occupancy at codons whose own identity is C
    asite[C] = mean normalised occupancy at codon k, attributed to the identity of codon k+1
  Reporting only one, without saying which, is a common way for codon-dwell numbers to become
  non-comparable between papers.

WHAT IT IS FOR
  The model never receives codon identity, tRNA abundance, or amino acid as an input -- only one-hot
  sequence, an ORF track and RNA coverage. So any codon structure recovered from `pred_flat` was
  learned from sequence. Whether that structure MATCHES the observed one, and whether it survives on
  transcripts the model never saw (mouse), is the actual question.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections import defaultdict

import numpy as np

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [c for c in CODONS if c not in STOPS]


def load_cds(path):
    cds = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                try:
                    cds[f[0]] = (int(f[1]), int(f[2]))
                except ValueError:
                    pass
    return cds


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = ln[1:].split("|")[0].split()[0]
                buf = []
            else:
                buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def occupancy(npz, cds, seqs, min_cds_counts, min_codons):
    """Return {source: {'psite': {codon: (sum, n)}, 'asite': ...}} plus per-source tx counts."""
    z = np.load(npz, allow_pickle=True)
    tx = [t.decode() if isinstance(t, bytes) else str(t) for t in z["tx_ids"]]
    lens = np.asarray(z["lengths"], dtype=np.int64)
    off = np.concatenate([[0], np.cumsum(lens)])
    out = {s: {k: defaultdict(lambda: [0.0, 0]) for k in ("psite", "asite")}
           for s in ("observed", "predicted")}
    used = {"observed": 0, "predicted": 0}
    arrays = {"observed": np.asarray(z["obs_flat"]), "predicted": np.asarray(z["pred_flat"])}

    for i, t in enumerate(tx):
        c = cds.get(t)
        s = seqs.get(t)
        if c is None or s is None:
            continue
        cds0, cds1 = c[0] - 1, c[1]          # transcripts_cds.txt is 1-based inclusive
        n_cod = (cds1 - cds0) // 3
        if n_cod < min_codons:
            continue
        if cds1 > len(s) or cds1 > lens[i]:
            continue
        codons = [s[cds0 + 3 * k: cds0 + 3 * k + 3].upper() for k in range(n_cod)]
        a, b = off[i] + cds0, off[i] + cds0 + 3 * n_cod

        # The signal gate applies to OBSERVED COUNTS ONLY. `pred_flat` is a per-transcript
        # probability distribution summing to exactly 1.0 (scale lives separately in `pred_total`),
        # so gating it on a count threshold rejected every transcript. Gating on observed also makes
        # both sides average over the IDENTICAL transcript set, without which the two 61-vectors
        # would not be comparable -- which matters more than the bug that surfaced it.
        if np.asarray(arrays["observed"][a:b], dtype=np.float64).sum() < min_cds_counts:
            continue

        for src, flat in arrays.items():
            v = np.asarray(flat[a:b], dtype=np.float64)
            per = v.reshape(n_cod, 3).sum(axis=1)     # per-codon signal
            m = per.mean()
            if m <= 0:
                continue
            per = per / m                              # normalise by THIS transcript's CDS mean
            used[src] += 1
            for k in range(n_cod):
                cod = codons[k]
                if len(cod) == 3 and cod in SENSE:
                    e = out[src]["psite"][cod]
                    e[0] += per[k]; e[1] += 1
                if k + 1 < n_cod:
                    nxt = codons[k + 1]
                    if len(nxt) == 3 and nxt in SENSE:
                        e = out[src]["asite"][nxt]
                        e[0] += per[k]; e[1] += 1
    return out, used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cds-counts", type=float, default=100.0)
    ap.add_argument("--min-codons", type=int, default=50)
    a = ap.parse_args()

    cds = load_cds(a.cds)
    seqs = load_fasta(a.fasta)
    occ, used = occupancy(pathlib.Path(a.npz), cds, seqs, a.min_cds_counts, a.min_codons)
    rows = []
    for src in ("observed", "predicted"):
        for site in ("psite", "asite"):
            for cod in SENSE:
                s, n = occ[src][site].get(cod, [0.0, 0])
                rows.append(dict(label=a.label, source=src, site=site, codon=cod,
                                 occupancy=(s / n if n else float("nan")), n=n))
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("label\tsource\tsite\tcodon\toccupancy\tn\n")
        for r in rows:
            fh.write(f"{r['label']}\t{r['source']}\t{r['site']}\t{r['codon']}\t"
                     f"{r['occupancy']:.6f}\t{r['n']}\n")
    print(f"  {a.label}: transcripts used obs={used['observed']:,} pred={used['predicted']:,} "
          f"-> {out.name}")
    # A quick within-dump sanity signal: obs vs pred correlation over the 61 A-site codons.
    o = np.array([occ["observed"]["asite"].get(c, [np.nan, 0])[0] /
                  max(occ["observed"]["asite"].get(c, [0, 1])[1], 1) for c in SENSE])
    p = np.array([occ["predicted"]["asite"].get(c, [np.nan, 0])[0] /
                  max(occ["predicted"]["asite"].get(c, [0, 1])[1], 1) for c in SENSE])
    ok = np.isfinite(o) & np.isfinite(p)
    if ok.sum() > 10:
        print(f"    A-site obs-vs-pred Pearson r = {np.corrcoef(o[ok], p[ok])[0,1]:.3f} "
              f"over {ok.sum()} codons")
    return 0


if __name__ == "__main__":
    sys.exit(main())
