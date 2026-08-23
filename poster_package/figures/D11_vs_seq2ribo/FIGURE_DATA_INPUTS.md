# FIGURE D11 -- vs seq2ribo: data inputs
**What it shows:** this model vs seq2ribo on WITHIN-transcript profile shape (per-codon CDS Shape r). This
model 0.526 >> seq2ribo published 0.05-0.19; ceiling ~0.95.
**Generator:** `make_vs_seq2ribo.py` (cas12a). THIS_MODEL=0.526 (orf_v2_attn fold-0, results.md per-codon
section); SEQ2RIBO=0.05-0.19 (Kaynar & Kingsford 2026 published within-transcript Shape r); ceiling from
`results/replicate_concordance.json` (pc_cds_codon, Spearman-Brown). Key point: seq2ribo's advertised 0.92 is
cross-transcript 'Elemwise r' (~ expression), not within-transcript shape -- so it is NOT the tool that
closes the shape gap. UPDATE THIS_MODEL on mm1 retrain.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/D11_vs_seq2ribo
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_vs_seq2ribo.py
```
