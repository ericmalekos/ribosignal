#!/usr/bin/env bash
# Submit each Chothani run to the decoy-aware salmon array as soon as BOTH its mates are
# md5-verified, instead of waiting for the whole 67 GB fetch. Download and compute overlap.
#
# Idempotent: a run is submitted once (tracked in submitted.txt) and skipped if its quant.sf
# already exists. Exits when all 25 are quantified or the fetch has stopped and nothing is left.
set -uo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
D=$NEW/data/external/chothani_rna
SUB=$D/logs/submitted.txt; touch "$SUB"
PY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
while true; do
  rows=$(PYTHONNOUSERSITE=1 $PY - <<'PYEOF'
import os
R="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
D=f"{R}/data/external/chothani_rna"
sub=set(x.strip() for x in open(f"{D}/logs/submitted.txt") if x.strip())
out=[]
for i,l in enumerate(open(f"{D}/manifest_25.tsv")):
    if i==0: continue
    f=l.rstrip("\n").split("\t")
    if len(f)<5 or f[0] in sub: continue
    if os.path.exists(f"{R}/data/tpm/chothani_decoy/{f[1]}/{f[0]}/quant.sf"): continue
    if all(os.path.exists(f"{D}/fastq/"+os.path.basename(p)) and
           os.path.getsize(f"{D}/fastq/"+os.path.basename(p))>=int(b)
           for p,b in zip(f[2].split(";"), f[4].split(";"))):
        out.append((i,f[0]))
print(";".join(f"{i}:{s}" for i,s in out))
PYEOF
)
  if [ -n "$rows" ]; then
    idx=$(echo "$rows" | tr ';' '\n' | cut -d: -f1 | paste -sd,)
    jid=$(sbatch --parsable --array="$idx" $NEW/scripts/salmon_chothani_decoy.sbatch 2>/dev/null)
    if [ -n "$jid" ]; then
      echo "[$(date +%T)] submitted $jid for rows $idx"
      echo "$rows" | tr ';' '\n' | cut -d: -f2 >> "$SUB"
    fi
  fi
  n=$(find $NEW/data/tpm/chothani_decoy -name quant.sf 2>/dev/null | wc -l)
  [ "$n" -ge 25 ] && { echo "[$(date +%T)] ALL 25 QUANTIFIED"; break; }
  # stop if the fetch is gone AND nothing further can become ready
  if ! ps -u emalekos -o comm= --no-headers | grep -qx wget; then
    [ -z "$rows" ] && [ "$(wc -l < "$SUB")" -ge 25 ] && { echo "[$(date +%T)] fetch done, all submitted"; break; }
  fi
  sleep 300
done
