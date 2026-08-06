#!/usr/bin/env python3
"""Proteogenomics UC-A step: enumerate candidate ORFs on A549's expressed transcripts and score each by
the MODEL's predicted per-nt profile (not RiboCode -- Task 19b showed its caller collapses to CTG). For
each transcript this enumerates every AUG..in-frame-stop ORF (same construction as eval_localization.py /
the orf_track), scores it by pred_frame0 (in-frame fraction of the predicted density -- the model's
translation-localization signal, AUROC ~0.94 on held-out truth) + predicted density, classifies it vs the
annotated CDS (canonical / uORF / dORF / N-ext / internal / lncRNA-ORF), and translates it to protein.

Output feeds the DB build: the NULL db = every candidate; the MODEL db = candidates the model scores as
translated (pred_frame0 above threshold). So both DBs are built from the same candidate set; the model's
value is the SELECTION. cas12a env (numpy only -- no torch; reads the pred_profiles.npz dump)."""
from __future__ import annotations

import argparse
from pathlib import Path

import sys
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
CODON = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _c1 in enumerate(_bases):
    for _j, _c2 in enumerate(_bases):
        for _k, _c3 in enumerate(_bases):
            CODON[_c1 + _c2 + _c3] = _aas[_i * 16 + _j * 4 + _k]


def translate(s):
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


# Near-cognate starts = the nine single-nucleotide variants of ATG, the set with documented
# mammalian initiation activity. The model's ORF track already SCORES these (CUG 0.5, GUG/ACG 0.35,
# UUG 0.3, AUA 0.25, AUU/AUC 0.2, AAG/AGG 0.15); until 2026-08-01 the enumerator emitted AUG only, so
# no non-AUG ORF could ever reach a search database.
NEAR_COGNATE = ("CTG", "GTG", "TTG", "ACG", "AGG", "AAG", "ATA", "ATT", "ATC")
ALL_STARTS = ("ATG",) + NEAR_COGNATE


def candidate_orfs(seq, min_nt, starts=("ATG",)):
    """All <start>..in-frame-stop ORFs -> (start0, end0_excl, start_codon); end0_excl past the stop.

    One ORF per (start position, next in-frame stop), so nested ORFs sharing a stop but opening at
    different starts are all emitted -- that is the point for non-AUG: a CUG upstream of the first
    in-frame AUG yields an N-terminally EXTENDED protein whose upstream tryptic peptides exist in no
    AUG-only database.
    """
    s = seq.upper().replace("U", "T")
    L = len(s)
    if L < 6:
        return []
    b = np.frombuffer(s.encode("ascii", "replace"), dtype=np.uint8)
    T_, A_, G_ = 84, 65, 71
    orfs = []
    for f in range(3):
        n = (L - f) // 3
        if n < 2:
            continue
        idx = f + 3 * np.arange(n)
        c0, c1, c2 = b[idx], b[idx + 1], b[idx + 2]
        is_stop = (((c0 == T_) & (c1 == A_) & (c2 == A_)) | ((c0 == T_) & (c1 == A_) & (c2 == G_))
                   | ((c0 == T_) & (c1 == G_) & (c2 == A_)))
        stop_pos = idx[is_stop]
        if stop_pos.size == 0:
            continue
        for sc in starts:
            k0, k1, k2 = (ord(sc[0]), ord(sc[1]), ord(sc[2]))
            m = (c0 == k0) & (c1 == k1) & (c2 == k2)
            for a in idx[m]:
                j = int(np.searchsorted(stop_pos, a, side="right"))
                if j >= stop_pos.size:
                    continue
                e = int(stop_pos[j]) + 3
                if e - a >= min_nt:
                    orfs.append((int(a), e, sc))
    return orfs


def cds_relationship(a, e, cds):
    if cds is None:
        return "no_cds"
    cs, ce = cds
    inframe = (a - cs) % 3 == 0
    if inframe and a >= ce:
        return "dorf_inframe"
    if inframe and a > cs:
        return "cds_inframe"
    if inframe and a == cs:
        return "cds_canonical"
    if max(a, cs) < min(e, ce):
        return "ext_inframe" if inframe else "cds_offframe"
    if e <= cs:
        return "utr5"
    if a >= ce:
        return "utr3"
    return "other"


# map (cds_relationship, biotype) -> a clean ORF class for the DB
def orf_class(crel, biotype):
    if crel == "cds_canonical":
        return "canonical"
    if biotype == "lncRNA" or crel == "no_cds":
        return "lncRNA_orf"
    if crel == "utr5":
        return "uORF"
    if crel in ("utr3", "dorf_inframe"):
        return "dORF"
    if crel == "ext_inframe":
        return "nterm_ext"
    if crel in ("cds_offframe", "cds_inframe"):
        return "internal"
    return "other"


def read_fasta(path, need):
    seqs, name, buf, keep = {}, None, [], False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None and keep:
                    seqs[name] = "".join(buf)
                name = ln[1:].split()[0].split("|")[0]; keep = name in need; buf = []
            elif keep:
                buf.append(ln.strip())
        if name is not None and keep:
            seqs[name] = "".join(buf)
    return seqs


def load_cds(path):
    out = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len")}
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); u5, cl = int(f[ix["utr5_len"]]), int(f[ix["cds_len"]])
            if cl > 0:
                out[f[ix["tx_id"]]] = (u5, u5 + cl)
    return out


def load_biotype(path):
    out = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ti, gi = h.index("tx_id"), h.index("gene_type")
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); out[f[ti]] = f[gi]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", default=str(NEW / "proteogenomics/data/A549_pilot/pred/pred_profiles.npz"))
    ap.add_argument("--fasta", default=str(NEW / "data/fibroblast_universe.fa"))
    ap.add_argument("--min_aa", type=int, default=10, help="min protein length (aa); ORF >= (min_aa+1)*3 nt")
    ap.add_argument("--out_dir", default=str(NEW / "proteogenomics/data/A549_pilot/orfs"))
    ap.add_argument("--tx2cds", default=str(NEW / "data/tx2cds.tsv"),
                    help="CDS table (default human; mouse: data/mouse_tx2cds.tsv)")
    ap.add_argument("--tx2biotype", default=str(NEW / "data/tx2biotype.tsv"),
                    help="biotype table (default human; mouse: data/mouse_tx2biotype.tsv)")
    ap.add_argument("--start_codons", default="ATG",
                    help="comma-separated start codons, or 'near_cognate' for ATG + the nine "
                         "single-nt variants (CTG,GTG,TTG,ACG,AGG,AAG,ATA,ATT,ATC). Default ATG "
                         "reproduces the pre-2026-08-01 behaviour exactly. Non-AUG ORFs are what make "
                         "the ncStart column live: they extend a protein N-terminally past the first "
                         "in-frame AUG, so their upstream peptides exist in no AUG-only database.")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    min_nt = (args.min_aa + 1) * 3
    starts = (ALL_STARTS if args.start_codons.strip().lower() == "near_cognate"
              else tuple(c.strip().upper().replace("U", "T") for c in args.start_codons.split(",") if c.strip()))
    bad = [c for c in starts if len(c) != 3 or set(c) - set("ACGT")]
    if bad:
        raise SystemExit(f"bad start codon(s): {bad}")
    print(f"start codons ({len(starts)}): {','.join(starts)}", file=sys.stderr)

    z = np.load(args.profiles, allow_pickle=False)
    tx_ids = [str(t) for t in z["tx_ids"]]
    lengths = z["lengths"]; pred = z["pred_flat"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    seqs = read_fasta(args.fasta, set(tx_ids))
    cds = load_cds(args.tx2cds)
    bt = load_biotype(args.tx2biotype)

    tsv = open(out / "candidates.tsv", "w")
    tsv.write("orf_id\ttx\tstart0\tend0\taa_len\torf_class\tbiotype\tpred_frame0\tpred_density"
              "\tstart_codon\n")
    faa = open(out / "candidates.faa", "w")
    from collections import Counter
    cls_n = Counter(); start_n = Counter(); n_orf = 0; n_tx = 0
    for i, tx in enumerate(tx_ids):
        seq = seqs.get(tx)
        L = int(lengths[i])
        if seq is None or len(seq) != L:
            continue
        p = pred[off[i]:off[i + 1]]
        cand = candidate_orfs(seq, min_nt, starts)
        if not cand:
            continue
        n_tx += 1
        biotype = bt.get(tx, "unknown")
        cbounds = cds.get(tx)
        for a, e, sc in cand:
            if e > L:
                continue
            pm = p[a:e]; psum = float(pm.sum())
            f0 = float(pm[0::3].sum() / psum) if psum > 0 else 0.0
            dens = float(pm.mean() * L)          # density scaled by tx length (comparable to eval_loc)
            klass = orf_class(cds_relationship(a, e, cbounds), biotype)
            prot = translate(seq[a:e - 3].upper().replace("U", "T"))   # drop the stop codon
            if len(prot) < args.min_aa or "*" in prot:   # drop truncated / internal-stop (shouldn't occur)
                continue
            oid = f"{tx}|{a}|{e}|{klass}"
            tsv.write(f"{oid}\t{tx}\t{a}\t{e}\t{len(prot)}\t{klass}\t{biotype}\t{f0:.4f}\t{dens:.3f}"
                      f"\t{sc}\n")
            # start= is consumed by build_a549_dbs to keep the NULL db AUG-only while the MODEL db may
            # select non-AUG ORFs -- the asymmetry is deliberate (the null is the naive baseline).
            faa.write(f">{oid} f0={f0:.3f} dens={dens:.2f} start={sc}\n{prot}\n")
            cls_n[klass] += 1; n_orf += 1; start_n[sc] += 1
    tsv.close(); faa.close()
    print(f"transcripts with candidates: {n_tx:,}; candidate ORFs written: {n_orf:,}")
    print("by class:", dict(cls_n.most_common()))
    print(f"-> {out}/candidates.tsv + candidates.faa")


if __name__ == "__main__":
    main()
