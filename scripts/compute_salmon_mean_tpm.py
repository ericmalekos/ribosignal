#!/usr/bin/env python3
"""Mean salmon TPM / NumReads per transcript across the 32 Fibroblast RNAseq samples.

Independent abundance reference for validating the posture-A RNAseq per-nt coverage
(per-transcript total coverage should correlate strongly with salmon TPM) and for later
defining the expressed training universe. Uses the existing decoy-aware salmon quant;
no realignment.

Output: data/tpm/fibroblast_salmon_mean_tpm.tsv  (tx_id, mean_tpm, mean_numreads, n_samples)
"""

# ============================================================================
# CANNOT RUN AS-IS -- input lost 2026-08-15.
#
# SAL/SALMON below points at expression_context_human/data/salmon_quant_chothani_decoy, which went
# with the deleted biotype_probe tree. There is NO surviving copy: no salmon quant anywhere under
# data/ holds the Chothani SRR15513* accessions (data/tpm/salmon_human_rna is a different cohort).
# It was deliberately NOT repointed, because there is nothing to point at.
#
# This script's OUTPUT survives, so nothing downstream is blocked today. What is lost is the
# ability to REBUILD that output from source. Re-deriving it means re-running salmon on the 69
# Chothani RNA libraries; the BAMs are still in data/rnaseq_bam_mm1/.
# ============================================================================
import sys
from collections import defaultdict
from pathlib import Path

ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/biotype_probe/expression_context_human")
NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
SALMON = ECH / "data" / "salmon_quant_chothani_decoy"
SRR_LIST = NEW / "data" / "fibroblast_rnaseq_srr.txt"
OUT = NEW / "data" / "tpm" / "fibroblast_salmon_mean_tpm.tsv"


def main():
    srrs = [s.strip() for s in SRR_LIST.read_text().splitlines() if s.strip()]
    tpm = defaultdict(float)
    nreads = defaultdict(float)
    n_used = 0
    for srr in srrs:
        qf = SALMON / srr / "quant.sf"
        if not qf.exists():
            print(f"WARN missing {qf}", file=sys.stderr)
            continue
        n_used += 1
        with qf.open() as fh:
            fh.readline()                       # header: Name Length EffectiveLength TPM NumReads
            for line in fh:
                f = line.rstrip("\n").split("\t")
                tx = f[0].split("|")[0]         # versioned ENST from the pipe header
                tpm[tx] += float(f[3])
                nreads[tx] += float(f[4])
    print(f"samples used: {n_used}/{len(srrs)}  transcripts: {len(tpm):,}", file=sys.stderr)

    with OUT.open("w") as o:
        o.write("tx_id\tmean_tpm\tmean_numreads\tn_samples\n")
        for tx in sorted(tpm):
            o.write(f"{tx}\t{tpm[tx] / n_used:.6f}\t{nreads[tx] / n_used:.3f}\t{n_used}\n")
    print(f"wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
