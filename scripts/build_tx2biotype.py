#!/usr/bin/env python3
"""Build a transcript -> (gene, biotype, chrom, length) table from the GENCODE v49 GTF.

One row per transcript. Length is the mature (spliced) transcript length = sum of
exon lengths, which matches the transcriptome-BAM @SQ LN and the per-nt P-site array
length used as the model target axis.

Output TSV columns (tab-separated, header included):
  tx_id  gene_id  gene_name  chrom  strand  transcript_type  gene_type  length

Reused by:
  - verify_target_inputs.py   (gene_name resolution for positive controls; biotype of ncRNA)
  - build_fibroblast_psite_target.py  (summary TSV biotype/gene join)
  - later held-out-chromosome split   (gene -> chrom)
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths  # noqa: E402

# The GTF is an input you supply; only the output has a sensible in-repo default.
DEFAULT_GTF = paths._env_path("RIBO_GTF",
                              paths._env_path("RIBO_ANNOT_DIR", paths.data_dir() / "annotations")
                              / "gencode.v49.annotation.gtf")
DEFAULT_OUT = paths.tx2biotype()

_ATTR = {k: re.compile(k + r' "([^"]+)"') for k in
         ("gene_id", "transcript_id", "gene_name", "gene_type", "transcript_type")}


def _get(pat, s):
    m = pat.search(s)
    return m.group(1) if m else "NA"


def main():
    ap = argparse.ArgumentParser(
        description="Build a tx -> (gene, biotype, chrom, length) table from a GENCODE GTF. "
                    "Defaults to human v49; pass --gtf/--out for other assemblies (e.g. mouse vM38).")
    ap.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    gtf, out_path = args.gtf, args.out
    if not gtf.is_file():
        sys.exit(paths.missing(gtf, "annotation GTF", "RIBO_GTF", "gtf"))
    out_path.parent.mkdir(parents=True, exist_ok=True)

    meta = {}          # tx_id -> [gene_id, gene_name, chrom, strand, transcript_type, gene_type]
    length = {}        # tx_id -> summed exon length
    n_tx = n_exon = 0

    with gtf.open() as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat == "transcript":
                attr = f[8]
                tx = _get(_ATTR["transcript_id"], attr)
                meta[tx] = [
                    _get(_ATTR["gene_id"], attr),
                    _get(_ATTR["gene_name"], attr),
                    f[0],                       # chrom
                    f[6],                       # strand
                    _get(_ATTR["transcript_type"], attr),
                    _get(_ATTR["gene_type"], attr),
                ]
                length.setdefault(tx, 0)
                n_tx += 1
            elif feat == "exon":
                tx = _get(_ATTR["transcript_id"], f[8])
                # GTF exon start/end are 1-based inclusive
                length[tx] = length.get(tx, 0) + (int(f[4]) - int(f[3]) + 1)
                n_exon += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        out.write("tx_id\tgene_id\tgene_name\tchrom\tstrand\t"
                  "transcript_type\tgene_type\tlength\n")
        for tx in sorted(meta):
            gid, gname, chrom, strand, ttype, gtype = meta[tx]
            out.write(f"{tx}\t{gid}\t{gname}\t{chrom}\t{strand}\t"
                      f"{ttype}\t{gtype}\t{length.get(tx, 0)}\n")

    print(f"transcripts: {n_tx:,}  exon-lines: {n_exon:,}", file=sys.stderr)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
