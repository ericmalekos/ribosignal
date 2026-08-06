#!/usr/bin/env bash
#SBATCH --job-name=macro_enum
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=0-11
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_enum_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_enum_%A_%a.err
set -euo pipefail
# Enumerate + model-score candidate ORFs per macrophage population (numpy-only, cas12a env; no GPU).
# Reads the deployed-model dump (pred_profiles.npz) + the population's mouse universe FASTA, anchors each
# candidate AUG..stop ORF against the MOUSE CDS/biotype tables, and writes candidates.tsv + candidates.faa.
# candidates.faa is the NULL search DB (all candidates); the MODEL DB is a pred_frame0 filter applied next.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
# TAG namespaces artefacts per checkpoint (see macro_dumps.sbatch). Empty = original single-model layout.
TAG=${TAG:-}
mkdir -p "$NEW/proteogenomics/logs"

mapfile -t POPS < <(ls "$PD"/universes/*_universe_tx.txt | xargs -n1 basename | sed 's/_universe_tx.txt//' | sort)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"
OUT="$PD/orfs${TAG}/$POP"
mkdir -p "$OUT"
echo "enumerate population=$POP tag=${TAG:-<none>} host=$(hostname) start=$(date -Iseconds)"

"$PY" "$NEW/proteogenomics/scripts/enumerate_score_orfs.py" \
  --profiles "$PD/pred${TAG}/$POP/pred_profiles.npz" \
  --fasta "$PD/universes/${POP}_universe.fa" \
  --tx2cds "$NEW/data/mouse_tx2cds.tsv" \
  --tx2biotype "$NEW/data/mouse_tx2biotype.tsv" \
  --min_aa 10 --out_dir "$OUT"

# guard: canonical ORFs must be a non-trivial fraction, else the mouse CDS anchor silently failed
"$PY" - "$OUT/candidates.tsv" <<'PY'
import sys
from collections import Counter
c = Counter(); n = 0
for ln in open(sys.argv[1]):
    if ln.startswith("orf_id"):
        continue
    c[ln.split("\t")[5]] += 1; n += 1
canon = c.get("canonical", 0)
print(f"{n:,} candidates; canonical={canon:,} ({canon/max(n,1):.1%}); classes={dict(c.most_common())}")
sys.exit(1 if canon < 1000 else 0)
PY
echo "DONE enumerate $POP $(date -Iseconds)"
