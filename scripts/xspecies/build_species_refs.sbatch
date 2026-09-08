#!/usr/bin/env bash
#SBATCH --job-name=xsp_refs
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=logs/xspecies/refs_%A_%a.out
#SBATCH --error=logs/xspecies/refs_%A_%a.err
# Log paths are relative to the submission directory; mkdir -p logs/xspecies first, or
# SLURM drops the job before it starts.
#
# Build every per-species reference artifact the cross-species pipeline needs.
# One array task per species:
#
#   1. normalize_annotation.py --write-gtf   -> normalized GTF + ncRNA drop list +
#                                               tx2gene + tx2biotype + 2 BEDs + report
#   2. gffread -w                            -> transcriptome FASTA
#   3. STAR --runMode genomeGenerate         -> shared STAR_indexes/
#   4. salmon index (decoy-aware)            -> shared Salmon_indexes/
#   5. RiboCode prepare_transcripts          -> data/ribocode_annot_<species>/
#
# WALLTIME IS A SCHEDULING DECISION, NOT A SAFETY MARGIN. This was first submitted
# with --time=23:45:00, and the jobs then sat PENDING for over two hours at
# Reason=Priority while phoenix-04/05 had 192 idle CPUs and 1.9 TB of free memory.
# SLURM backfill will only start a job in a gap it can FINISH inside, and essentially
# no gap is 24 hours long, so an oversized walltime makes a job unschedulable no
# matter how small it actually is. Resubmitting the identical job with --time=02:00:00
# started it in under 90 seconds. Measured usage, for sizing:
#
#     worm  0:01:43 /  1 GB      macaque  1:13:38 / 73 GB
#     fly   0:02:02 /  4 GB      chimp    0:46:57 / 85 GB
#     fish  0:27:33 / 41 GB      gorilla  0:50:57 / 74 GB
#
# The defaults below cover the worst case (gorilla/chimp). For the small genomes,
# override on the command line and they will backfill almost instantly:
#     sbatch --array=0-2 --time=00:20:00 --mem=16G  ...   # yeast, worm, fly
#
# Launch:  sbatch --array=0-6 scripts/xspecies/build_species_refs.sbatch
#          sbatch --array=0   scripts/xspecies/build_species_refs.sbatch    # yeast only, smoke test
#
# Every step has an idempotent skip guard, so a re-run after a timeout resumes.
#
# WHY THE TRANSCRIPTOME COMES FROM gffread AND NOT NCBI's _rna.fna:
# NCBI's rna_from_genomic headers are `lcl|NC_001133.9_mrna_NM_001180043.1_1`, not bare
# accessions, so they do not join to the GTF or tx2gene without rewriting -- a silent
# failure mode for salmon. Deriving from the NORMALIZED GTF guarantees the transcript IDs
# are byte-identical to the ones STAR, RiboCode and the tx2gene map use.
set -euo pipefail

# ---- locations ---------------------------------------------------------------------
# Nothing below is specific to one cluster. RIBO_REPO defaults to the directory sbatch
# was launched from, which is the repo root if you follow the launch line above; the
# tool prefix defaults to the conda environment that is active at submit time. Set
# RIBO_ENV / RIBO_REPO explicitly when either assumption does not hold.
REPO=${RIBO_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}
ENV=${RIBO_ENV:-${CONDA_PREFIX:?set RIBO_ENV=<prefix of the env holding STAR/salmon/samtools/cutadapt>}}
DATA=${RIBO_DATA_DIR:-$REPO/data}
REFS=${RIBO_SPECIES_REFS:?set RIBO_SPECIES_REFS=<dir holding <species>/genome.fna and genomic.gtf>}
STARIDX=${RIBO_STAR_INDEXES:?set RIBO_STAR_INDEXES=<dir to write STAR indexes into>}
SALIDX=${RIBO_SALMON_INDEXES:?set RIBO_SALMON_INDEXES=<dir to write salmon indexes into>}
GFFREAD=${RIBO_GFFREAD:-gffread}
RIBOCODE=${RIBO_RIBOCODE_BIN:-$ENV/bin}
T=${SLURM_CPUS_PER_TASK:-16}

# The per-species parameter table is scripts/xspecies/species.tsv, read here rather than
# duplicated. An earlier version hardcoded the same seven rows in this file AND in that
# one, which is exactly how an annotation version drifts between the index that gets
# built and the index that gets used.
#
# Array index = row order in species.tsv (comments and header skipped):
#   0 yeast  1 celegans  2 fly  3 zebrafish  4 macaque  5 chimp  6 gorilla
#   7 human_refseq  8 mouse_refseq
# List them with:  awk -F'\t' '!/^#/ && NR>1{print NR-2, $1}' scripts/xspecies/species.tsv
SPTSV=$REPO/scripts/xspecies/species.tsv
[ -s "$SPTSV" ] || { echo "missing $SPTSV" >&2; exit 2; }
mapfile -t SPECS < <(awk -F'\t' '!/^#/ && $1!="species" && NF>1 {print $1"\t"$2"\t"$7"\t"$8}' "$SPTSV")

IDX=${SLURM_ARRAY_TASK_ID:?set --array}
[ "$IDX" -lt "${#SPECS[@]}" ] || { echo "array index $IDX exceeds ${#SPECS[@]} species in $SPTSV" >&2; exit 2; }
IFS=$'\t' read -r SP ANN OVERHANG SAN <<< "${SPECS[$IDX]}"

GENOME="$REFS/$SP/genome.fna"
RAWGTF="$REFS/$SP/genomic.gtf"
ANNOTD="$DATA/annot/$SP"
NORMGTF="$ANNOTD/$SP.normalized.gtf"
TXFA="$ANNOTD/transcripts.fa"
SIDX="$STARIDX/star_index_${SP}_${ANN}"
QIDX="$SALIDX/salmon_index_decoy_${SP}_${ANN}"
RCANNOT="$DATA/ribocode_annot_$SP"

SCR=${RIBO_TMPDIR:-${TMPDIR:-/tmp}}/xsp_refs_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR" "$ANNOTD"
trap 'rm -rf "$SCR"' EXIT

echo "=== $SP  ann=$ANN  overhang=$OVERHANG  SAn=$SAN"
echo "    host=$(hostname) cpus=$T start=$(date -Iseconds)"
[ -s "$GENOME" ] || { echo "missing genome $GENOME" >&2; exit 2; }
[ -s "$RAWGTF" ] || { echo "missing gtf $RAWGTF" >&2; exit 2; }

# ---- 1. normalize -----------------------------------------------------------
if [ -s "$NORMGTF" ] && [ -s "$ANNOTD/ncrna_tx.txt" ]; then
  echo "[1/5] normalize: SKIP (exists, $(wc -l < "$ANNOTD/ncrna_tx.txt") ncRNA tx)"
else
  echo "[1/5] normalize $(date +%H:%M:%S)"
  "$ENV/bin/python" "$REPO/scripts/normalize_annotation.py" \
    --species "$SP" --gtf "$RAWGTF" --genome "$GENOME" --write-gtf --out-dir "$ANNOTD"
fi
NNC=$(wc -l < "$ANNOTD/ncrna_tx.txt")
[ "$NNC" -ge 10 ] || { echo "ABORT: ncrna_tx.txt has only $NNC entries" >&2; exit 3; }

# ---- 2. transcriptome FASTA -------------------------------------------------
if [ -s "$TXFA" ]; then
  echo "[2/5] gffread: SKIP ($(grep -c '^>' "$TXFA") transcripts)"
else
  echo "[2/5] gffread $(date +%H:%M:%S)"
  "$GFFREAD" -w "$SCR/tx.fa" -g "$GENOME" "$NORMGTF"
  # gffread appends ' gene=NAME' etc to the header; keep only the transcript id so
  # salmon's target names join to tx2gene and the GTF without any rewriting.
  awk '/^>/{print $1; next}{print}' "$SCR/tx.fa" > "$TXFA"
  rm -f "$SCR/tx.fa"
fi
NTX=$(grep -c '^>' "$TXFA")
echo "      transcriptome: $NTX sequences"
# Every transcriptome sequence must appear in tx2gene, or salmon output cannot be
# aggregated and the universe build will silently lose transcripts.
MISS=$("$ENV/bin/python" - "$TXFA" "$ANNOTD/tx_to_gene.tsv" <<'PY'
import sys
fa={l[1:].split()[0] for l in open(sys.argv[1]) if l.startswith('>')}
t2g={l.split('\t')[0] for l in open(sys.argv[2])}
print(len(fa-t2g))
PY
)
[ "$MISS" = "0" ] || { echo "ABORT: $MISS transcriptome IDs are absent from tx_to_gene.tsv" >&2; exit 4; }
echo "      all $NTX transcript IDs join to tx_to_gene.tsv"

# ---- 3. STAR index ----------------------------------------------------------
if [ -s "$SIDX/SA" ] && [ -s "$SIDX/genomeParameters.txt" ]; then
  echo "[3/5] STAR index: SKIP ($SIDX)"
else
  echo "[3/5] STAR genomeGenerate $(date +%H:%M:%S)"
  mkdir -p "$SIDX"
  # cd into node-local scratch so STAR's _STARtmp lands there, not on ceph.
  # --limitGenomeGenerateRAM default is 31 GB and is too low for a 3.5 Gb genome.
  # --limitSjdbInsertNsj default is 1e6; RefSeq mammalian annotation exceeds it.
  ( cd "$SCR" && "$ENV/bin/STAR" --runMode genomeGenerate --runThreadN "$T" \
      --genomeDir "$SIDX" --genomeFastaFiles "$GENOME" --sjdbGTFfile "$NORMGTF" \
      --sjdbOverhang "$OVERHANG" --genomeSAindexNbases "$SAN" \
      --limitGenomeGenerateRAM 90000000000 --limitSjdbInsertNsj 4000000 )
  [ -s "$SIDX/SA" ] || { echo "ABORT: STAR index has no SA" >&2; exit 5; }
fi
grep -E 'sjdbOverhang|genomeSAindexNbases' "$SIDX/genomeParameters.txt" | sed 's/^/      /'

# ---- 4. salmon decoy-aware index --------------------------------------------
if [ -s "$QIDX/info.json" ]; then
  echo "[4/5] salmon index: SKIP ($QIDX)"
else
  echo "[4/5] salmon index $(date +%H:%M:%S)"
  # The gentrome and decoys list were NOT kept for the human/mouse indexes, so the
  # build could not be reproduced. Keep them beside the index this time.
  mkdir -p "$QIDX"
  grep '^>' "$GENOME" | cut -d' ' -f1 | sed 's/^>//' > "$QIDX/decoys.txt"
  cat "$TXFA" "$GENOME" > "$SCR/gentrome.fa"
  "$ENV/bin/salmon" index -t "$SCR/gentrome.fa" -d "$QIDX/decoys.txt" \
      -i "$QIDX" -k 31 -p "$T" --tmpdir "$SCR/saltmp"
  [ -s "$QIDX/info.json" ] || { echo "ABORT: salmon index incomplete" >&2; exit 6; }
fi
"$ENV/bin/python" -c "import json,sys;d=json.load(open('$QIDX/info.json'));\
print('      salmon: num_refs=%s num_decoys=%s k=%s'%(d.get('num_refs'),d.get('num_decoys'),d.get('k')))"

# ---- 5. RiboCode annotation --------------------------------------------------
if [ -s "$RCANNOT/transcripts_sequence.fa" ]; then
  echo "[5/5] RiboCode annot: SKIP ($RCANNOT)"
else
  echo "[5/5] RiboCode prepare_transcripts $(date +%H:%M:%S)"
  "$RIBOCODE/prepare_transcripts" -g "$NORMGTF" -f "$GENOME" -o "$RCANNOT"
  [ -s "$RCANNOT/transcripts_sequence.fa" ] || { echo "ABORT: RiboCode annot empty" >&2; exit 7; }
fi
NRC=$(grep -c '^>' "$RCANNOT/transcripts_sequence.fa" || echo 0)
echo "      RiboCode annot: $NRC transcripts (GTF has $NTX)"
# A large shortfall means the GTF shape confused prepare_transcripts -- the FlyBase
# `mRNA`-instead-of-`transcript` problem the normalizer exists to fix.
"$ENV/bin/python" -c "
n,g=$NRC,$NTX
r=n/g if g else 0
print('      ratio %.3f'%r)
raise SystemExit(0 if r>=0.5 else 'ABORT: RiboCode kept only %.1f%% of transcripts'%(100*r))"

echo "=== $SP DONE=$(date -Iseconds)"
du -sh "$SIDX" "$QIDX" "$RCANNOT" "$ANNOTD" 2>/dev/null | sed 's/^/    /'
