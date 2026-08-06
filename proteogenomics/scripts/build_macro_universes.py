#!/usr/bin/env python3
"""Per-population expressed universe + FASTA for the 12 macrophage populations (per-population
DB build). Per population: mean-TPM-merge its 3 salmon quants, then build_line_universe.py --species
mouse (TPM>=1, pc+lncRNA, <=10k nt, no chrM). Writes universes/<pop>_universe_tx.txt + _universe.fa."""
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
PD = NEW / "proteogenomics" / "data" / "macrophage_tissue"
BLU = NEW / "proteogenomics" / "scripts" / "build_line_universe.py"
CASPY = "/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3"
OUT = PD / "universes"
OUT.mkdir(exist_ok=True)

SHORT = {
    "BMDM": "BMDM", "Kupffer": "Kupffer", "Large intestinal macrophages": "LargeIntestinal",
    "Liver_recruited macrophages": "LiverRecruited", "Lung_recruited macrophages": "LungRecruited",
    "Lung_resident macrophages": "LungResident", "Microglia": "Microglia",
    "Peritoneal macrophages": "Peritoneal", "Raw264.7": "RAW264",
    "Small intestinal macrophages": "SmallIntestinal",
    "Spleen_recruited macrophages": "SpleenRecruited",
    "Spleen_resident macrophages": "SpleenResident",
}

rows = list(csv.DictReader(open(PD / "rnaseq_manifest.tsv"), delimiter="\t"))
pops = defaultdict(list)
for r in rows:
    pops[r["population"]].append(r["run"])


def merge_salmon(srrs, out):
    """mean-TPM merge of the population's salmon quants (TPM=mean over reps, NumReads=sum)."""
    data = {}
    n = 0
    for srr in srrs:
        q = PD / "salmon" / srr / "quant.sf"
        if not q.exists():
            print(f"    WARN missing {q}")
            continue
        n += 1
        with open(q) as f:
            f.readline()
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                name, L, eff, tpm, nr = p[0], p[1], p[2], float(p[3]), float(p[4])
                d = data.setdefault(name, [L, eff, 0.0, 0.0])
                d[2] += tpm
                d[3] += nr
    with open(out, "w") as w:
        w.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for name, (L, eff, tsum, nr) in data.items():
            w.write(f"{name}\t{L}\t{eff}\t{tsum / max(n, 1):.6f}\t{nr:.3f}\n")
    return len(data)


for pop, srrs in sorted(pops.items()):
    short = SHORT[pop]
    merged = OUT / f"{short}_merged_quant.sf"
    ntx = merge_salmon(srrs, merged)
    print(f"[{short}] merged {len(srrs)} reps ({ntx} tx) -> build universe")
    subprocess.run([CASPY, str(BLU), "--salmon", str(merged), "--out_prefix", str(OUT / short),
                    "--species", "mouse", "--min_tpm", "1.0", "--max_len", "10000"], check=True)
print("ALL 12 population universes built ->", OUT)
