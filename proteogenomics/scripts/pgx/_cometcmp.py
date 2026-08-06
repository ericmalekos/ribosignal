from collections import Counter
from pathlib import Path
from pgx.report import load_rank1, cut, canonical_at_global_fdr, novel_at_class_fdr
X = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/data/macrophage_tissue/pgx_xsubtype")
for eng, base in [("MSFragger (calibrate_mass=2)", X/"search/Microglia"),
                  ("MSFragger (calibrate_mass=0)", X/"calibtest/nocal"),
                  ("Comet (single pass)",          X/"comet/Microglia")]:
    print(f"=== {eng}", flush=True)
    res = {}
    for a in ["gencode", "model_ribocode_bmdm_nt"]:
        d = base/a
        if not (d/"rank1.tsv.gz").exists() and not list(d.glob("*.tsv")):
            print(f"  {a}: missing"); continue
        rows = load_rank1(d)
        k = Counter(r[2] for r in rows)
        c, _ = canonical_at_global_fdr(rows)
        g = cut(rows, {"canon_t","novel_t"}, {"novel_d","other_d"})
        np_, peps, _, _ = novel_at_class_fdr(rows)
        res[a] = (len(rows), c, len(peps))
        print(f"  {a:24s} rank1={len(rows):9,d} canon_t={k['canon_t']:9,d} "
              f"cut={g:7.3f} canon@1%={c:8,d} novel_pep={len(peps):4d}", flush=True)
    if len(res) == 2:
        g_, m_ = res["gencode"], res["model_ribocode_bmdm_nt"]
        print(f"  --> delta rank1 = {m_[0]-g_[0]:+,}    dPSM = {m_[1]-g_[1]:+,}", flush=True)
