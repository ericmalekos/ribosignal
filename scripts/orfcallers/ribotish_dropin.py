#!/usr/bin/env python3
"""Ribo-TISH drop-in: run `ribotish predict` on the model's PREDICTED per-nt profiles.

Same contract as scripts/ribocode_dropin.py, so the two callers can be compared on identical
density arrays. Three variants, identical caller and parameters, ONLY the density differs:

    real            observed P-sites (truth)
    pred_obsdepth   predicted shape scaled to the transcript's OBSERVED total
    pred_preddepth  predicted shape scaled to the count head's pred_total * --pred_scale
                    (standalone; --pred_poisson draws counts instead of rounding)

WHY THIS WORKS WITHOUT A BAM. `ribotish predict --inprofile` reads a per-transcript P-site
profile instead of opening bam files (predict.py:212-220). Its Ribo object is indexed 0-based
over the mature transcript (`self.length = trans.cdna_length()`, ribo.py:56-66) with
`nhead = ntail = 0`, which is the same axis as the pack, so a density array drops straight in.
The BAM path is skipped only when the profile covers EVERY transcript of a gene
(predict.py:381-385, `load = False`), so the GTF must be restricted to exactly the transcripts
written here -- pass --gtf a subset, not the full annotation.

FILE FORMAT (predict.py:217, ribo.py:338-352):
    Gid \\t Tid \\t Symbol \\t {pos:count, ...} \\t {pos:count, ...}
column 3 is the TIS profile (written empty: no TIS/LTM data here), column 4 the RPF profile.
Both are sparse python dict literals over 0-based cDNA positions and are `eval`ed.

LENGTH CHECK IS NOT OPTIONAL. Ribo-TISH derives cDNA length from the GTF exons; the pack derives
it from the sequence it was built on. If they disagree the profile is silently shifted, so every
transcript's GTF length is compared against the pack length and any mismatch is dropped and
counted rather than written.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

NEW = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(NEW / "scripts"))
from ribocode_dropin import (  # noqa: E402  same density code path as RiboCode
    build_density,
    validate_density_flags,
)

ATTR = re.compile(r'(\S+) "([^"]*)"')


def gtf_index(gtf: Path):
    """tx_id -> (gene_id, symbol, cdna_length) from exon lines."""
    exlen, meta = defaultdict(int), {}
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            t = a.get("transcript_id")
            if not t:
                continue
            exlen[t] += int(f[4]) - int(f[3]) + 1
            if t not in meta:
                meta[t] = (a.get("gene_id", ""), a.get("gene_name", a.get("gene_id", "")))
    return {t: (meta[t][0], meta[t][1], exlen[t]) for t in exlen}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True,
                    help="pred_profiles.npz from dump_pred_profiles.py")
    ap.add_argument("--gtf", required=True, help="GTF RESTRICTED to the transcripts being scored")
    ap.add_argument("--genome", required=True, help="genome fasta (ribotish -f)")
    ap.add_argument("--variant", required=True,
                    choices=["real", "pred_obsdepth", "pred_preddepth"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--pred_scale", type=float, default=1.0)
    ap.add_argument("--pred_poisson", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minaalen", type=int, default=5,
                    help="ribotish --minaalen (RiboCode uses min_aa 5)")
    ap.add_argument("--fpth", type=float, default=0.05, help="ribotish frame p threshold")
    ap.add_argument("--numproc", type=int, default=4)
    ap.add_argument("--ribotish", default="ribotish",
                    help="ribotish binary; the default resolves on $PATH. Ribo-TISH does not "
                         "have to share this python env -- point at its own env's binary.")
    ap.add_argument("--emit-gtf", default=None,
                    help="write a GTF restricted to exactly the transcripts profiled here, and "
                         "use it instead of --gtf. Ribo-TISH only skips its bam path when the "
                         "profile covers EVERY transcript of a gene, so a GTF filtered on any "
                         "other list (the pack, say) leaves genes with unprofiled transcripts and "
                         "Ribo-TISH tries to open the bam. Filtering here removes that class of "
                         "error, because this is the only place that knows what was written.")
    ap.add_argument("--longest", action="store_true",
                    help="ribotish --longest: one ORF per stop codon, the analogue of RiboCode's "
                         "*_collapsed.txt keying on (gene, ORF_gstop). Without it Ribo-TISH also "
                         "reports every in-frame downstream ATG as a separate Truncated ORF.")
    ap.add_argument("--profile_only", action="store_true", help="write the inprofile, do not run")
    a = ap.parse_args()
    try:
        validate_density_flags(a.variant, pred_scale=a.pred_scale, poisson=a.pred_poisson)
    except ValueError as e:
        ap.error(str(e))

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tag = a.variant
    if a.pred_scale != 1.0:
        tag += f"_depth{a.pred_scale:g}"
    if a.pred_poisson:
        tag += f"_poisson{a.seed}"

    z = np.load(a.profiles, allow_pickle=True)
    ids = list(z["tx_ids"]); lens = z["lengths"]
    off = np.concatenate([[0], np.cumsum(lens)])
    pf, ob, pt = z["pred_flat"], z["obs_flat"], z["pred_total"]
    rng = np.random.default_rng(a.seed) if a.pred_poisson else None

    gi = gtf_index(Path(a.gtf))
    print(f"GTF carries {len(gi):,} transcripts; npz carries {len(ids):,}", file=sys.stderr)

    n_written = n_absent = n_lenmismatch = n_empty = 0
    written = []
    inprof = out / f"{tag}_inprofile.txt"
    with inprof.open("w") as fh:
        fh.write("Gid\tTid\tSymbol\tTISProf\tRiboProf\n")
        for j, t in enumerate(ids):
            g = gi.get(t)
            if g is None:
                n_absent += 1
                continue
            gid, sym, glen = g
            aa, bb = off[j], off[j + 1]
            if glen != int(lens[j]):
                n_lenmismatch += 1
                continue
            dens = build_density(a.variant, pf[aa:bb], ob[aa:bb], pt[j],
                                 pred_scale=a.pred_scale, rng=rng, poisson=a.pred_poisson)
            nz = np.nonzero(dens)[0]
            if nz.size == 0:
                n_empty += 1
            body = ", ".join(f"{int(i)}:{int(dens[i])}" for i in nz)
            fh.write(f"{gid}\t{t}\t{sym}\t{{}}\t{{{body}}}\n")
            written.append(t)
            n_written += 1
    print(f"wrote {inprof}  tx={n_written:,}  absent_from_gtf={n_absent:,}  "
          f"length_mismatch={n_lenmismatch:,}  all_zero={n_empty:,}", file=sys.stderr)
    if n_lenmismatch:
        print(f"WARNING: {n_lenmismatch} transcripts dropped on cDNA-length mismatch between the "
              f"GTF and the pack. A silent mismatch would shift the profile.", file=sys.stderr)
    if a.profile_only:
        return 0

    gtf_used = a.gtf
    if a.emit_gtf:
        # Keep only transcripts that were actually profiled, and only genes ALL of whose
        # transcripts survived that filter. A gene with even one unprofiled transcript must go,
        # or Ribo-TISH will read a bam for it.
        keep_tx = set(written)
        gene_tx, lines = defaultdict(set), []
        with Path(a.gtf).open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                at = dict(ATTR.findall(f[8]))
                g, t = at.get("gene_id"), at.get("transcript_id")
                if g and t:
                    gene_tx[g].add(t)
                    lines.append((g, t, line))
        whole = {g for g, ts in gene_tx.items() if ts <= keep_tx}
        with Path(a.emit_gtf).open("w") as out:
            n = sum(bool(out.write(ln)) for g, t, ln in lines if g in whole and t in keep_tx)
        n_tx = len({t for g, t, _ in lines if g in whole and t in keep_tx})
        print(f"wrote {a.emit_gtf}: {len(whole):,} of {len(gene_tx):,} genes fully profiled, "
              f"{n_tx:,} transcripts, {n:,} lines", file=sys.stderr)
        if not whole:
            print("FATAL: no gene is fully covered by the profile; Ribo-TISH would read a bam "
                  "for every one of them.", file=sys.stderr)
            return 2
        gtf_used = a.emit_gtf

    # predict.py:82-84 hard-exits unless -t or -b is non-empty, even though --inprofile
    # supplies everything. The path is never opened: multiRiboGene is inside `if load:`
    # (predict.py:385-388) and load is False once the profile covers every transcript of the
    # gene, and find_offset only probes for a sibling .para.py. A DELIBERATELY NONEXISTENT path
    # is used rather than an empty BAM: if `load` ever became True this crashes loudly instead
    # of silently scoring every transcript at zero counts.
    cmd = [a.ribotish, "predict", "-g", str(gtf_used), "-f", str(a.genome),
           "-b", "/nonexistent/inprofile_only__no_bam_should_be_read.bam",
           "--inprofile", str(inprof), "-o", str(out / f"{tag}.txt"),
           "--minaalen", str(a.minaalen), "--fpth", str(a.fpth),
           "-p", str(a.numproc), "-v"]
    if a.longest:
        cmd.append("--longest")
    print("RUN " + " ".join(cmd), file=sys.stderr)
    r = subprocess.run(cmd)
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
