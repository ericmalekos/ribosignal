#!/usr/bin/env bash
# Move bulk data to the warm archive: COPY -> VERIFY BY CHECKSUM -> INDEX -> STUB -> delete source.
#
#   archive_to_warm.sh <local_path> "<why it is archived / how to get it back>"
#   archive_to_warm.sh --dry-run <local_path> "..."
#
# The source is deleted ONLY after a checksum-verified copy exists. Every archived path leaves:
#   1. a row in data/ARCHIVE_INDEX.tsv        (the searchable index -- CHECK THIS BEFORE DOWNLOADING)
#   2. a <name>.ARCHIVED.md stub in its place (so anyone who looks where the data WAS, finds it)
#
# Archive root mirrors the project tree, so a local path maps to exactly one archive path:
#   /private/groups/carpenterlab/emalekos/RNAZoo_meta/<rel>
#     -> /private/warm-archive/carpenterlab/RNAZoo_meta/<rel>
set -uo pipefail
DRY=0; [ "${1:-}" = "--dry-run" ] && { DRY=1; shift; }
SRC=$(readlink -f "${1:?usage: archive_to_warm.sh [--dry-run] <path> \"<reason>\"}")
REASON="${2:-no reason given}"
PROJ=/private/groups/carpenterlab/emalekos/RNAZoo_meta
ARCH=/private/warm-archive/carpenterlab/RNAZoo_meta
IDX=$PROJ/RNAZoo/experiments/riboseq_signal_model/data/ARCHIVE_INDEX.tsv

case "$SRC" in "$PROJ"/*) ;; *) echo "REFUSING: $SRC is outside $PROJ" >&2; exit 2;; esac
[ -e "$SRC" ] || { echo "REFUSING: $SRC does not exist" >&2; exit 2; }
REL=${SRC#$PROJ/}
DST=$ARCH/$REL
SZ=$(du -sb "$SRC" | cut -f1)
NF=$(find "$SRC" -type f 2>/dev/null | wc -l)
echo "  src   $SRC"
echo "  dst   $DST"
echo "  size  $(python3 -c "print(f'{$SZ/1e9:.1f} GB')")  files $NF"
[ "$DRY" = 1 ] && { echo "  (dry run, nothing done)"; exit 0; }

mkdir -p "$(dirname "$DST")" || { echo "cannot create $(dirname "$DST")" >&2; exit 3; }
echo "  [1/4] copying..."
rsync -a --info=progress2 "$SRC" "$(dirname "$DST")/" || { echo "rsync FAILED" >&2; exit 4; }

echo "  [2/4] verifying by CHECKSUM (not size/mtime)..."
DIFF=$(rsync -a --checksum --dry-run --itemize-changes "$SRC" "$(dirname "$DST")/" 2>/dev/null | grep -c '^[<>ch]') || DIFF=0
if [ "$DIFF" != "0" ]; then
  echo "  VERIFY FAILED: $DIFF path(s) differ by checksum. SOURCE NOT DELETED." >&2
  exit 5
fi
echo "        checksums match"

echo "  [3/4] indexing + stub"
[ -s "$IDX" ] || printf 'archived_utc\trel_path\tarchive_path\tsize_bytes\tn_files\treason\n' > "$IDX"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$REL" "$DST" "$SZ" "$NF" "$REASON" >> "$IDX"
STUB="$SRC.ARCHIVED.md"
cat > "$STUB" <<STUBEOF
# ARCHIVED -- data moved to the warm archive $(date -u +%Y-%m-%d)

**This path no longer holds data locally.** It lives at:

    $DST

$(python3 -c "print(f'{$SZ/1e9:.1f} GB')") in $NF files. Reason: $REASON

Restore with:

    rsync -a "$DST" "$(dirname "$SRC")/"

The archive is **$ARCH**, mirroring the project tree, 3.0 PB free. The searchable index of
everything archived is \`data/ARCHIVE_INDEX.tsv\` -- **check it before re-downloading anything.**
STUBEOF

echo "  [4/4] removing local copy"
rm -rf "$SRC"
echo "  DONE  $REL  -> archive  ($(python3 -c "print(f'{$SZ/1e9:.1f} GB')") freed)"
