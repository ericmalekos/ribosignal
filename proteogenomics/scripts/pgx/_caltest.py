from collections import Counter
from pathlib import Path
from pgx.report import load_rank1, cut, canonical_at_global_fdr
X = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/macrophage_tissue/pgx_xsubtype")
for lab, base in [("calibrate_mass=2 (default)", X/"search/Microglia"),
                  ("calibrate_mass=0 (test)",    X/"calibtest/nocal")]:
    print(f"=== {lab}", flush=True)
    res = {}
    for a in ["gencode", "model_ribocode_bmdm_nt"]:
        rows = load_rank1(base/a)
        k = Counter(r[2] for r in rows)
        c, _ = canonical_at_global_fdr(rows)
        g = cut(rows, {"canon_t","novel_t"}, {"novel_d","other_d"})
        res[a] = (len(rows), k["canon_t"], g, c)
        print(f"  {a:24s} rank1={len(rows):9,d} canon_t={k['canon_t']:9,d} cut={g:6.2f} canon@1%={c:9,d}", flush=True)
    d_rank1 = res["model_ribocode_bmdm_nt"][0] - res["gencode"][0]
    d_psm   = res["model_ribocode_bmdm_nt"][3] - res["gencode"][3]
    print(f"  --> delta rank1 = {d_rank1:+,}   dPSM = {d_psm:+,}", flush=True)
