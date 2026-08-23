#!/usr/bin/env bash
# Re-fetch the 57 Chothani RNA-seq runs (~193 GB) whose salmon quant died with the deleted
# biotype_probe tree. HEAD NODE, SERIAL, SINGLE CONNECTION, md5-verified -- per the standing rules:
# compute nodes throttle a single stream to ~0.4 MB/s and parallel array downloads trip ENA's
# per-IP limit; multi-connection (aria2c -x4) silently corrupts large files at the correct size.
#
# Resume policy: a file that already matches its ENA md5 is skipped. A file that FAILS md5 is
# DELETED and refetched from scratch -- `wget -c` would keep the bad bytes and append to them.
#
#   bash scripts/fetch_chothani_rna.sh
set -uo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
D=$NEW/data/external/chothani_rna
MAN=${MAN:-$D/manifest_25.tsv}   # default: the 25 runs whose TPM was lost (Fibroblast survived)
OUT=$D/fastq
mkdir -p "$OUT" "$D/logs"
LOG=$D/logs/fetch.log
ok=0; skip=0; fail=0
while IFS=$'\t' read -r run tissue ftp md5s bytes; do
  [ "$run" = "run" ] && continue
  IFS=';' read -ra F <<< "$ftp"; IFS=';' read -ra M <<< "$md5s"
  for i in "${!F[@]}"; do
    url="https://${F[$i]}"; want="${M[$i]}"; f="$OUT/$(basename "${F[$i]}")"
    if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ]; then
      skip=$((skip+1)); continue
    fi
    IFS=';' read -ra B <<< "$bytes"; want_b="${B[$i]}"
    if [ -e "$f" ] && [ "$(stat -c %s "$f")" -ge "$want_b" ]; then
      rm -f "$f"                       # complete but md5-bad: resume would keep the bad bytes
      echo "  [refetch] $(basename "$f") was complete but failed md5" | tee -a "$LOG"
    fi
    echo "[$(date +%T)] $run $tissue $(basename "$f") ($((want_b/1000000)) MB)" | tee -a "$LOG"
    wget -q -c --tries=5 --timeout=60 --waitretry=10 -O "$f" "$url" 2>>"$LOG"
    got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)
    if [ "$got" = "$want" ]; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); rm -f "$f"
      echo "  MD5 MISMATCH $run $(basename "$f") want=$want got=${got:-none}" | tee -a "$LOG"
    fi
  done
done < "$MAN"
echo "[$(date +%T)] DONE downloaded=$ok already_ok=$skip failed=$fail" | tee -a "$LOG"
df_free=$(getfattr -n ceph.dir.rbytes --only-values /private/groups/carpenterlab 2>/dev/null)
echo "  quota now: $(python3 -c "print(f'{$df_free/1e12:.2f} TB ({100*$df_free/15e12:.1f}%)')")" | tee -a "$LOG"
[ "$fail" -eq 0 ]
