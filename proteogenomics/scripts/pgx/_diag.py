from collections import Counter
from pathlib import Path
from pgx.report import load_rank1, cut, canonical_at_global_fdr
X = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/macrophage_tissue/pgx_xsubtype")
for pop in ["Microglia", "RAW264", "BMDM"]:
    print(f"=== {pop}", flush=True)
    for a in ["gencode", "model_ribocode_bmdm_nt"]:
        rows = load_rank1(X / f"search/{pop}/{a}")
        k = Counter(r[2] for r in rows)
        g = cut(rows, {"canon_t", "novel_t"}, {"novel_d", "other_d"})
        c, _ = canonical_at_global_fdr(rows)
        print(f"  {a:24s} rank1={len(rows):8,d} canon_t={k['canon_t']:8,d} "
              f"other_d={k['other_d']:8,d} novel_t={k['novel_t']:6,d} cut={g:6.2f} canon@1%={c:8,d}", flush=True)
