#!/usr/bin/env python3
"""Expressed-transcript universe for the `pgx` pipeline.

Same rule as the validated `build_line_universe.py`, lifted here so the model arm and BOTH null
arms are built from ONE universe object, so the null contrast compares like with like:

    TPM >= --min-tpm  AND  transcript_type in {protein_coding, lncRNA}
    AND  chrom != chrM  AND  mature length <= --max-len

Several salmon quants are merged by MEAN TPM across replicates (as `build_macro_universes.py` does),
so a population's reps collapse to one universe. Writes <prefix>_universe_tx.txt (versioned ENST,
one per line) and <prefix>_universe.fa (headers = versioned ENST).

  python -m pgx.universe --salmon a/quant.sf b/quant.sf --species mouse --out-prefix <dir>/BMDM
"""
from __future__ import annotations

import argparse
from pathlib import Path

from .refs import require
from .seqtools import iter_fasta, read_fasta


def mean_tpm(salmon_paths) -> dict:
    """{versioned_tx: mean TPM across the given quant.sf files}. Missing in a rep counts as 0."""
    tot, n = {}, 0
    for q in salmon_paths:
        q = Path(q)
        if not q.exists():
            raise SystemExit(f"missing salmon quant: {q}")
        n += 1
        with open(q) as fh:
            fh.readline()
            for ln in fh:
                p = ln.split("\t")
                tot[p[0].split("|")[0]] = tot.get(p[0].split("|")[0], 0.0) + float(p[3])
    return {tx: v / max(n, 1) for tx, v in tot.items()}


def load_tx_meta(tx2biotype, field="transcript_type") -> dict:
    """{tx_id: (chrom, <field>, length)} from the tx2biotype table."""
    meta = {}
    with open(tx2biotype) as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        if field not in ci:
            raise SystemExit(f"{tx2biotype} has no column {field!r}; has {sorted(ci)}")
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            meta[p[ci["tx_id"]]] = (p[ci["chrom"]], p[ci[field]], int(p[ci["length"]]))
    return meta


KEEP_BT = {"protein_coding", "lncRNA"}
BIOTYPE_FIELDS = ("transcript_type", "gene_type")


def build(salmon_paths, species, min_tpm=1.0, max_len=10000,
          biotype_field="transcript_type", **ref_overrides):
    """-> (universe:set[tx], tpm:dict[tx->float], meta:dict) applying the documented filters.

    biotype_field picks WHICH biotype the {protein_coding, lncRNA} filter applies to:

      transcript_type (default)  only canonical coding isoforms. Matches the validated
                                 `build_line_universe.py` rule and every existing pack/universe in
                                 this project, so it is the continuity-preserving choice.
      gene_type                  every isoform of a protein-coding or lncRNA GENE, which also
                                 admits nonsense_mediated_decay, retained_intron, TEC and
                                 non_stop_decay transcripts.

    The two differ materially and the difference is proteogenomically relevant: on A549 the
    gene_type universe is 45,166 vs 41,139 transcripts, and the 4,027 extra are 3,747 NMD, 106
    retained_intron, 52 TEC and 21 non_stop_decay. Those isoforms are exactly where non-canonical
    ORFs are enriched, so gene_type is the more inclusive choice for discovery; transcript_type is
    the more conservative one and keeps comparability with prior results.
    """
    refs = require(species, **ref_overrides)
    tpm = mean_tpm(salmon_paths)
    meta = load_tx_meta(refs["tx2biotype"], biotype_field)
    universe = {tx for tx, v in tpm.items()
                if v >= min_tpm and (m := meta.get(tx)) and m[1] in KEEP_BT
                and m[0] != "chrM" and m[2] <= max_len}
    return universe, tpm, meta


def write(universe, species, out_prefix, **ref_overrides):
    """Write <prefix>_universe_tx.txt + <prefix>_universe.fa. -> (tx_path, fa_path, n_fa)."""
    refs = require(species, **ref_overrides)
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tx_path = Path(f"{out_prefix}_universe_tx.txt")
    tx_path.write_text("\n".join(sorted(universe)) + "\n")
    fa_path = Path(f"{out_prefix}_universe.fa")
    n_fa = 0
    with open(fa_path, "w") as w:
        for fa in refs["tx_fastas"]:
            for name, seq in iter_fasta(fa):
                if name in universe:
                    w.write(f">{name}\n{seq}\n"); n_fa += 1
    return tx_path, fa_path, n_fa


def read_universe_fasta(fa_path, keep=None) -> dict:
    """{tx: sequence} from a universe FASTA (optionally restricted to `keep`)."""
    return read_fasta(fa_path, keep)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--salmon", nargs="+", required=True, help="one or more quant.sf (mean-TPM merged)")
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out-prefix", required=True, help="e.g. <dir>/BMDM -> _universe_tx.txt / _universe.fa")
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--max-len", type=int, default=10000)
    ap.add_argument("--biotype-field", default="transcript_type", choices=BIOTYPE_FIELDS,
                    help="which biotype the {protein_coding, lncRNA} filter applies to. "
                         "transcript_type (default) matches build_line_universe.py and every "
                         "existing universe in this project; gene_type additionally admits NMD / "
                         "retained_intron / TEC isoforms of coding genes, where non-canonical ORFs "
                         "are enriched (A549: 45,166 vs 41,139 transcripts). See build().")
    ap.add_argument("--tx2biotype", default=None, help="override the species tx2biotype table")
    a = ap.parse_args()

    ov = {k: v for k, v in (("tx2biotype", a.tx2biotype),) if v}
    universe, tpm, _ = build(a.salmon, a.species, a.min_tpm, a.max_len, a.biotype_field, **ov)
    tx_path, fa_path, n_fa = write(universe, a.species, a.out_prefix, **ov)
    n_expr = sum(1 for v in tpm.values() if v >= a.min_tpm)
    print(f"quants merged: {len(a.salmon)}  expressed(TPM>={a.min_tpm}): {n_expr:,}  "
          f"universe(pc/lncRNA, <={a.max_len} nt, no chrM): {len(universe):,}  FASTA: {n_fa:,}")
    print(f"wrote {tx_path}\nwrote {fa_path}")
    if n_fa != len(universe):
        print(f"WARNING: {len(universe) - n_fa} universe tx missing from FASTA (id mismatch?)")


if __name__ == "__main__":
    main()
