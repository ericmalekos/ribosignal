#!/usr/bin/env bash
# Identify the 3' adapter and read-length profile of every downloaded FASTQ, PER RUN,
# by measuring rather than assuming, and write the answer to a TSV the samplesheet
# generator reads.
#
# WHY A PANEL AND NOT JUST TruSeq. The standing cluster rule names AGATCGGAAGAGC,
# and every Ribo-seq dataset this project had processed before used it. The worm arm
# does not: GSE52905 is 83.2% `TGGAATTCTCGGGTGCCAAGG` (Illumina small-RNA) and
# **0.0%** TruSeq. Assuming TruSeq there would have trimmed nothing, left ~20 nt of
# adapter on every 51 nt read, and then `-m 20 -M 40` would have DISCARDED those
# reads for being too long -- producing a near-empty BAM rather than an error. That
# is the same failure shape as the Janich library that once mapped 0.03%.
#
# WHY PER RUN AND NOT PER DATASET. PRJEB65856 alone carries two trimming regimes
# inside one study, split by read length: the 35 nt human and macaque-iPSC-CM Ribo
# runs are pre-clipped and N-padded, the 51 nt runs still carry their adapter.
#
# WHY `zcat | head` AND NOT A BYTE RANGE. The rule is to probe the COMPLETE file. A
# streamed byte-range slice once reported 0% where the real file was 93%, because a
# mid-file gzip offset does not decompress to record boundaries. `zcat | head`
# decompresses from the start and head closes the pipe, so it costs seconds.
#
# Usage: bash scripts/xspecies/probe_adapters.sh [--force] [label ...]
set -uo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
FQROOT="$NEW/data/external/xspecies/fastq"
MANIFEST="$NEW/data/external/xspecies/download_manifest.tsv"
OUT="$NEW/data/external/xspecies/adapter_probe.tsv"
NREADS=100000

# name:sequence. All are scored; the winner is the LONGEST that clears MIN_PCT (see below).
PANEL=(
  "truseq:AGATCGGAAGAGC"                 # Illumina TruSeq / NEBNext
  "smallrna:TGGAATTCTCGGGTGCCAAGG"       # Illumina TruSeq Small RNA; the worm arm
  "ingolia:CTGTAGGCACCATCAAT"            # classic Ingolia ribosome-profiling linker
  "ingolia_short:CTGTAGGCACC"            # conserved prefix of the Ingolia linker. The zebrafish
                                         # arm (GSE46512) uses a variant that diverges after ~11 nt,
                                         # so the full 17-mer scored 0.0% by exact match here while
                                         # cutadapt -- which tolerates mismatches -- trimmed 74% of
                                         # reads with it. Keep both.
  "smallrna_short:TGGAATTCTCGGG"         # 13 nt prefix of the small-RNA adapter
  "polyA:AAAAAAAAAAAA"
)
MIN_PCT=20        # below this, no adapter is called present

[ -s "$MANIFEST" ] || { echo "missing $MANIFEST" >&2; exit 2; }
FORCE=0
[ "${1:-}" = "--force" ] && { FORCE=1; shift; }
# --force with labels re-probes ONLY those labels. It used to delete the whole table, which
# silently discarded 117 already-probed runs and made the samplesheet generator refuse
# everything. Scope the reset to what was asked for.
if [ "$FORCE" = 1 ]; then
  if [ "$#" -gt 0 ] && [ -s "$OUT" ]; then
    keep=$(mktemp); head -1 "$OUT" > "$keep"
    awk -F'\t' -v labs="$*" 'BEGIN{n=split(labs,a," "); for(i=1;i<=n;i++) drop[a[i]]=1}
                              NR>1 && !($2 in drop)' "$OUT" >> "$keep"
    mv "$keep" "$OUT"
    echo "re-probing only: $* (kept $(( $(wc -l < "$OUT") - 1 )) other rows)"
  else
    rm -f "$OUT"
  fi
fi
if [ ! -s "$OUT" ]; then
  printf 'run\tlabel\tspecies\tassay\tfile\tn_probed\tbest_adapter\tbest_seq\tbest_pct\tall_pcts\tmodal_len\tmin_len\tmax_len\tnpad_pct\tverdict\tcutadapt_arg\n' > "$OUT"
fi

WANT=("$@")
want_label() {
  [ "${#WANT[@]}" -eq 0 ] && return 0
  for w in "${WANT[@]}"; do [ "$1" = "$w" ] && return 0; done
  return 1
}

while IFS=$'\t' read -r label species assay run; do
  want_label "$label" || continue
  for f in "$FQROOT/$label/$run"*.fastq.gz; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    # R2 shares the library chemistry with R1; probing one mate is enough.
    case "$base" in *_2.fastq.gz) continue;; esac
    awk -F'\t' -v r="$run" -v b="$base" 'NR>1 && $1==r && $5==b{found=1} END{exit !found}' "$OUT" && continue

    SEQF="$(mktemp)"
    zcat "$f" 2>/dev/null | head -n $((NREADS*4)) | awk 'NR%4==2' > "$SEQF"
    n_tot=$(wc -l < "$SEQF")
    [ "${n_tot:-0}" -gt 0 ] || { echo "  WARN $base: no reads" >&2; rm -f "$SEQF"; continue; }

    # Score every panel member, then pick the LONGEST sequence that clears MIN_PCT --
    # not the highest-scoring one. A 13 nt prefix of an adapter always matches at least
    # as often as the full 21 nt version (it also matches reads where the adapter runs
    # off the read end), so scoring by percentage alone would hand cutadapt a truncated
    # adapter. Measured on the worm arm: smallrna_short 95.9% vs full smallrna 83.2%,
    # same adapter. cutadapt matches partial 3' adapters itself, so the full sequence
    # is strictly the better argument.
    best_name=none; best_seq=""; best_pct=0; best_len=0; all=""
    for entry in "${PANEL[@]}"; do
      nm=${entry%%:*}; sq=${entry#*:}
      c=$(grep -c "$sq" "$SEQF" || true)
      p=$(awk -v c="$c" -v t="$n_tot" 'BEGIN{printf "%.1f",100*c/t}')
      all="${all}${all:+,}${nm}=${p}"
      if awk -v a="$p" -v m="$MIN_PCT" 'BEGIN{exit !(a>m)}'; then
        L=${#sq}
        if [ "$L" -gt "$best_len" ]; then
          best_len=$L; best_pct=$p; best_name=$nm; best_seq=$sq
        fi
      fi
    done

    read -r modal mn mx npad <<< "$(awk '
      { L=length($0); len[L]++; if(mn==""||L<mn)mn=L; if(L>mx)mx=L; if($0 ~ /NN$/) p++ }
      END { best=0; bl=0; for (L in len) if (len[L]>best){best=len[L]; bl=L}
            printf "%d %d %d %.1f", bl, mn, mx, 100*p/NR }' "$SEQF")"
    rm -f "$SEQF"

    # Verdict, in priority order:
    #   <adapter>  a panel member is present in >MIN_PCT of reads -> cutadapt -a <seq>
    #   trim-n     little adapter but the reads are N-padded, i.e. the submitter
    #              already clipped them (the Chothani/Ruiz-Orera format). Keyed on
    #              N-padding ALONE: the human arm is 99.7% padded but its lengths
    #              run 35-51, so requiring a fixed length misses it.
    #   none       little adapter, no padding -> already trimmed upstream
    if [ "$best_name" != "none" ]; then
      verdict=$best_name; cut_arg=$best_seq
    elif awk -v p="$npad" 'BEGIN{exit !(p>5)}'; then
      verdict=trim-n;     cut_arg="--trim-n"; best_name=none; best_seq=""
    else
      verdict=none;       cut_arg="--trim-n"; best_name=none; best_seq=""
    fi

    # HARD GUARD. This probe greps for EXACT adapter substrings; cutadapt matches with an
    # error tolerance and will find adapters this probe misses. So "no adapter detected" is
    # not proof there is none. The dangerous combination is a FIXED-LENGTH read longer than
    # the Ribo length window (-M 40): if the verdict is trim-n/none, cutadapt discards EVERY
    # read and the arm silently vanishes. Measured on zebrafish GSE46512: fixed 44 nt,
    # exact-match adapter 0.0%, and trimming with --trim-n alone left 75 of 200,000 reads.
    if [ "$assay" = "ribo" ] && [ "$mn" = "$mx" ] && [ "$mn" -gt 40 ] && [ "$verdict" != "truseq" ] \
       && [ "$best_name" = "none" ]; then
      echo "    *** WARNING $base: fixed ${mn} nt reads, longer than the -M 40 Ribo window, and no" >&2
      echo "        adapter matched exactly. cutadapt WILL discard essentially every read. Find the" >&2
      echo "        real adapter from the 3' end of the reads before aligning this run. ***" >&2
      verdict="NEEDS_ADAPTER"; cut_arg="UNRESOLVED"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n' \
      "$run" "$label" "$species" "$assay" "$base" "$n_tot" "$best_name" "$best_seq" \
      "$best_pct" "$all" "$modal" "$mn" "$mx" "$npad" "$verdict" "$cut_arg" >> "$OUT"
    printf '  %-12s %-22s %-14s %5s%%  len %3d-%-3d (modal %3d) Npad=%5s%%  -> %s\n' \
      "$run" "$label" "$best_name" "$best_pct" "$mn" "$mx" "$modal" "$npad" "$verdict"
  done
done < <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3"\t"$6}' "$MANIFEST")

echo
echo "=== verdict by label ==="
awk -F'\t' 'NR>1{c[$2"\t"$4"\t"$15]++} END{for(k in c) printf "%-30s %3d\n", k, c[k]}' "$OUT" | sort
echo
echo "=== any label with MIXED verdicts (a per-dataset setting would be wrong there) ==="
awk -F'\t' 'NR>1{v[$2]=v[$2]" "$15} END{for(l in v){n=split(v[l],a," "); u=""; for(i=1;i<=n;i++) if(index(u,a[i])==0) u=u" "a[i]; if(split(u,b," ")>1) print "  "l": "u}}' "$OUT"
echo "wrote $OUT"
