#!/usr/bin/env bash
#SBATCH --job-name=screen_polya
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/liver_polya/screen_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/liver_polya/screen_%j.err
set -uo pipefail
# GATE for the GSE302188 poly(A) liver RNA arm: run the exact test GSE243134's totalRNA failed, on the
# first replicate, before committing to the remaining downloads and the full pack build.
#
# Metadata cannot settle this (SRA reports library_selection=cDNA for both poly(A) AND Ribo-Zero
# libraries), so the library is judged on its own mapping statistics.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
IDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
FQ=$NEW/data/external/mouse_liver_polya_rnaseq/fastq
QD=$NEW/data/salmon_quant_gse302188_liver
RUN=${SCREEN_RUN:-SRR34449471}
mkdir -p "$QD"

R1=$FQ/${RUN}_1.fastq.gz; R2=$FQ/${RUN}_2.fastq.gz
[ -s "$R1" ] && [ -s "$R2" ] || { echo "MISSING mates for $RUN ($R1 / $R2)"; exit 1; }

echo "=== adapter check (full local file) ==="
A=$(zcat "$R1" | head -400000 | awk 'NR%4==2' | grep -c AGATCGGAAGAGC)
L=$(zcat "$R1" | head -4000 | awk 'NR%4==2{print length($0)}' | sort -n | awk '{v[NR]=$1}END{print v[int(NR/2)]}')
echo "  $RUN read_len=$L  adapter=$(echo "$A" | awk '{printf "%.1f%%",$1/1000}')  ($A/100000 reads)"

if [ ! -s "$QD/$RUN/quant.sf" ]; then
  T=/data/tmp/emalekos/gse302188_trim/$RUN; mkdir -p "$T"
  "$ENV/bin/cutadapt" -j 16 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 20 \
    -o "$T/1.fq.gz" -p "$T/2.fq.gz" "$R1" "$R2" > "$T/cut.log" 2>&1
  grep -E "Reads written|Total read pairs" "$T/cut.log" | sed 's/^/  /'
  "$ENV/bin/salmon" quant -i "$IDX" -l A -1 "$T/1.fq.gz" -2 "$T/2.fq.gz" -p 16 \
    --seqBias --gcBias -o "$QD/$RUN" > "$QD/$RUN.log" 2>&1
  rm -rf "$T"
fi

echo
echo "===== SCREEN VERDICT ====="
SCREEN_RUN=$RUN "$PY" - <<'PY'
import json, os
NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
QD  = f"{NEW}/data/salmon_quant_gse302188_liver"
RUN = os.environ.get("SCREEN_RUN", "SRR34449471")
name = {}
with open(f"{NEW}/data/mouse_tx2biotype.tsv") as f:
    h = f.readline().rstrip("\n").split("\t")
    ti, gi, bi = h.index("tx_id"), h.index("gene_name"), h.index("gene_type")
    for ln in f:
        p = ln.rstrip("\n").split("\t"); name[p[ti]] = (p[gi], p[bi])

rows = []
with open(f"{QD}/{RUN}/quant.sf") as f:
    f.readline()
    for ln in f:
        p = ln.rstrip("\n").split("\t"); rows.append((p[0].split("|")[0], float(p[3])))
rows.sort(key=lambda r: -r[1]); tot = sum(t for _, t in rows) or 1.0
rate = json.load(open(f"{QD}/{RUN}/aux_info/meta_info.json")).get("percent_mapped", 0.0)
n1   = sum(1 for _, t in rows if t >= 1)
top10 = sum(t for _, t in rows[:10]) / tot * 100

print(f"{'library':<40} {'map%':>7} {'TPM>=1':>9} {'top10%':>8}   verdict")
print(f"{'-'*40} {'-'*7} {'-'*9} {'-'*8}   {'-'*7}")
print(f"{'Wang P42_Liver_mRNA  poly(A) REFERENCE':<40} {92.0:>7.1f} {14331:>9,} {32.9:>8.1f}   pass")
print(f"{'GSE243134 totalRNA   REJECTED':<40} {17.7:>7.1f} {3534:>9,} {93.6:>8.1f}   FAIL")
ok = rate >= 70 and top10 <= 50 and n1 >= 10000
print(f"{'GSE302188 ' + RUN + '  poly(A)':<40} {rate:>7.1f} {n1:>9,} {top10:>8.1f}   "
      f"{'**PASS**' if ok else '**FAIL**'}")
print(f"\n  pass bar: map>=70  TPM>=1 >=10,000  top10<=50")
print(f"\n  top 10 transcripts ({RUN}):")
for t, v in rows[:10]:
    g, b = name.get(t, (t.split(".")[0], "?"))
    print(f"    {g:<14} {b:<16} TPM={v:>10,.0f}")
print(f"\n  GATE: {'PROCEED to full download + pack build' if ok else 'STOP -- do not build a universe from this library'}")
PY
