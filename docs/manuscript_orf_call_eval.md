# RiboCode ORF-call evaluation of model variants (manuscript results)

Per-model, standing-rule ORF-call metrics on the **Chothani Hepatocytes hold-out** for every architecture /
embedding / universe variant. This is the manuscript's "how do the actual ORF calls change" section, and every
input path + count is recorded so the numbers are fully reproducible. Populated automatically by
`scripts/heldout/orf_call_metrics.py` (one `###` block per model), driven by the eval orchestrators.

## Inputs (what each model is + where its calls come from)
- **Model families** (all: 2 transformer layers unless noted, mm1 unique-mapper RNA coverage, no-Kozak ORF
  track, trained on 7 tissues holding out Hepatocytes; seed 0):
  - *Architecture sweep* (one-hot input): `onehot 2-attn` (deployed recipe), `attn4` (4 transformer layers),
    `mamba` (2 bidirectional Mamba layers).
  - *FM-embedding sweep* (input embedding swapped for one-hot): `rinalmo` (1280-d), `orthrus` (512-d),
    `hydrarna` (1024-d).
  - *Universe*: `fibroblast` = the 36,668-tx Fibroblast-expressed universe (salmon TPM>=1 in Fibroblast);
    `union` = the 84,472-tx union universe (TPM>=1 in ANY of the 8 tissues; adds 47,817 tx incl. 2,158
    Hepatocytes-exclusive). Each model is scored against ITS OWN universe's observed calls (noted per block).
- **Reference (truth)** for each model = its `dropin/real_collapsed.txt` -- RiboCode `detectORF` run on the
  REAL Hepatocytes P-site density (the same pack target the model is trained against), on that model's universe.
- **Predicted call sets** per model (`scripts/dump_pred_profiles.py` -> `scripts/ribocode_dropin.py`):
  - `pred_obsdepth` -- predicted per-nt profile scaled to the REAL per-tx P-site total (isolates SHAPE, i.e.
    whether the ORF distribution is right, independent of the count head).
  - `pred_preddepth` -- fully standalone / de novo: predicted profile x predicted count-head total, theta=1.
  - Poisson CDS-anchored -- `pred_preddepth` re-called after `--pred_poisson --pred_scale theta` over a theta
    sweep; the reported point is the theta whose CDS recall is nearest 0.90 (the two-arm partner of theta=1).

## Methods
- RiboCode `detectORF` on a SUBSTITUTED per-nt density (`ribocode_dropin.py`), `--min_aa 5 --pval 0.05`,
  annotation `expression_context_human/data/ribocode_annot` (GENCODE v49). ORFs matched between predicted and
  observed on the **genomic coord key** `gstart_gstop_len` (annotation-version robust; `orf_call_metrics.py`).
- **Minimum ORF length.** Calls are made at the permissive floor `--min_aa 5` (a 5-aa / ~18-nt ORF) to include
  the shortest micro-ORFs; RiboCode's own signal gate (frame-0 sum >=5 and >=5 non-zero frame-0 codons) removes
  signal-less tiny candidates. Because 5 aa is well below RiboCode's default (20 aa) and the tiny-ORF tail is
  uORF-heavy and the most over-call-prone, every model additionally reports a **length sweep**: the standalone
  theta=1 calls are re-scored after restricting BOTH the predicted and observed sets to ORFs >= {5,10,15,20,30,40}
  aa (post-hoc length stratification, `orf_call_metrics.py --min-aa-sweep`; aa = ORF_length_nt//3 - 1). This
  holds the calling procedure fixed and shows how precision/recall/F1 depend on the length floor (it is NOT a
  re-run of detectORF at each floor, so the FDR pool is constant -- a cleaner length-dependence measure).
- **Primary floor = 20 aa** (RiboCode's own default; `--primary-min-aa 20`). The headline P/R/F1 tables and the
  Poisson arm are computed at >= 20 aa; the sweep below each block shows the full length-dependence. The primary
  is a post-hoc >= 20 aa cut of the min_aa=5 calls, which was verified against a NATIVE `detectORF --min_AA 20`
  run on the baseline: Jaccard 0.99 (real 0.9925, pred 0.9887), all per-class P/R/F1 within 0.003, and the
  post-hoc set is a strict subset of the native set (~1% fewer calls -- the larger min_aa=5 FDR pool makes BH
  marginally stricter, so the post-hoc primary is if anything slightly conservative). Hence one min_aa=5 call
  set feeds both the primary (post-hoc >= 20 aa) and the sweep, with no separate 20-aa calling run.
- **Two-arm protocol** (standing rule): every model reports the standalone `theta=1` deterministic arm AND the
  Poisson CDS-anchored arm. theta=1 exposes de novo over-calling; the Poisson arm is the calibrated operating
  point (trades non-canonical recall for precision at a fixed CDS-recall target).
- Per ORF class: **annotated** = canonical CDS (trusted), **uORF**, **novel**, **dORF**, plus the
  **non-canonical** aggregate (everything but CDS). Metrics: precision (frac of predicted calls in the observed
  set), recall (frac of observed recovered), F1, and the **over-call ratio** n_predicted / n_observed.
- Profile fidelity (`val_pearson`, per-nt Pearson on the held-out val fold) is reported alongside for context;
  it is NOT the ORF-call metric -- a model can have lower val_pearson yet call ORFs as well or better.

## How to read the tables
`CDS P/R/F1` = precision/recall/F1 for canonical ORFs; `uORF/novel/dORF P/R` = precision/recall per class;
`non-canon F1` = F1 over all non-CDS ORFs. Higher precision at a fixed CDS recall (Poisson arm) = fewer
spurious non-canonical calls. The over-call line gives the raw n_predicted vs n_observed inflation at theta=1.

---

## Results (per model, appended as each finishes)


### onehot 2-attn (baseline)  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5268 (e6)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.960/0.974/0.967 | 0.444/0.877 | 0.371/0.868 | 0.072/0.533 | 0.480 |
| pred_preddepth (standalone theta=1) | 0.936/0.954/0.945 | 0.509/0.818 | 0.397/0.824 | 0.122/0.378 | 0.535 |
- standalone theta=1 over-call: total 15,630 calls vs 13,942 real (all 1.121x; novel 2.075x, dORF 3.1x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.933/0.924/0.929 | 0.719/0.214 | 0.605/0.260 | 0.292/0.117 | 0.294 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,872/15,415 | 0.945 | 0.447/0.847/0.585 | 0.347/0.812/0.486 | 0.107/0.361/0.165 | 0.410/0.708/0.519 |
| 10 | 17,596/14,988 | 0.945 | 0.487/0.836/0.616 | 0.364/0.813/0.503 | 0.114/0.368/0.175 | 0.436/0.694/0.535 |
| 15 | 16,541/14,463 | 0.945 | 0.498/0.821/0.620 | 0.386/0.826/0.526 | 0.116/0.367/0.176 | 0.442/0.677/0.535 |
| 20 | 15,630/13,942 | 0.945 | 0.509/0.818/0.627 | 0.397/0.824/0.536 | 0.122/0.378/0.184 | 0.447/0.667/0.535 |
| 30 | 14,252/13,110 | 0.945 | 0.525/0.815/0.639 | 0.427/0.852/0.569 | 0.127/0.384/0.191 | 0.455/0.651/0.536 |
| 40 | 13,326/12,492 | 0.946 | 0.531/0.791/0.636 | 0.435/0.872/0.581 | 0.107/0.346/0.164 | 0.446/0.627/0.521 |

### onehot mamba  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5171 (e4)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.962/0.974/0.968 | 0.424/0.890 | 0.372/0.863 | 0.071/0.500 | 0.480 |
| pred_preddepth (standalone theta=1) | 0.932/0.944/0.938 | 0.507/0.793 | 0.401/0.811 | 0.146/0.289 | 0.541 |
- standalone theta=1 over-call: total 15,224 calls vs 13,942 real (all 1.092x; novel 2.02x, dORF 1.972x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.90) | 0.929/0.901/0.915 | 0.756/0.184 | 0.581/0.229 | 0.342/0.072 | 0.251 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,333/15,415 | 0.938 | 0.454/0.821/0.585 | 0.347/0.810/0.486 | 0.121/0.296/0.172 | 0.425/0.688/0.525 |
| 10 | 17,157/14,988 | 0.938 | 0.489/0.812/0.611 | 0.364/0.811/0.502 | 0.134/0.297/0.184 | 0.451/0.676/0.541 |
| 15 | 16,125/14,463 | 0.938 | 0.498/0.797/0.613 | 0.385/0.820/0.524 | 0.144/0.296/0.194 | 0.460/0.658/0.541 |
| 20 | 15,224/13,942 | 0.938 | 0.507/0.793/0.618 | 0.401/0.811/0.537 | 0.146/0.289/0.194 | 0.467/0.643/0.541 |
| 30 | 13,887/13,110 | 0.938 | 0.524/0.783/0.628 | 0.434/0.833/0.571 | 0.152/0.281/0.198 | 0.482/0.622/0.543 |
| 40 | 13,040/12,492 | 0.938 | 0.512/0.761/0.612 | 0.445/0.843/0.582 | 0.149/0.280/0.195 | 0.477/0.601/0.532 |

### orthrus  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_orthrus_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5109 (e4)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.959/0.971/0.965 | 0.415/0.892 | 0.345/0.881 | 0.054/0.522 | 0.443 |
| pred_preddepth (standalone theta=1) | 0.935/0.950/0.942 | 0.495/0.819 | 0.372/0.822 | 0.129/0.350 | 0.523 |
- standalone theta=1 over-call: total 15,577 calls vs 13,942 real (all 1.117x; novel 2.207x, dORF 2.722x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.90) | 0.930/0.904/0.917 | 0.753/0.163 | 0.505/0.220 | 0.365/0.106 | 0.231 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 19,220/15,415 | 0.942 | 0.422/0.850/0.564 | 0.311/0.816/0.450 | 0.111/0.370/0.170 | 0.385/0.699/0.497 |
| 10 | 17,792/14,988 | 0.942 | 0.466/0.843/0.600 | 0.331/0.817/0.471 | 0.122/0.378/0.185 | 0.416/0.687/0.518 |
| 15 | 16,584/14,463 | 0.942 | 0.483/0.827/0.610 | 0.353/0.824/0.494 | 0.130/0.378/0.194 | 0.429/0.668/0.523 |
| 20 | 15,577/13,942 | 0.942 | 0.495/0.819/0.617 | 0.372/0.822/0.512 | 0.129/0.350/0.188 | 0.437/0.651/0.523 |
| 30 | 14,050/13,110 | 0.943 | 0.519/0.792/0.627 | 0.406/0.841/0.548 | 0.139/0.329/0.196 | 0.455/0.618/0.524 |
| 40 | 13,125/12,492 | 0.943 | 0.522/0.779/0.625 | 0.420/0.853/0.563 | 0.135/0.299/0.185 | 0.457/0.595/0.517 |

### hydrarna  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_hydrarna_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5102 (e6)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.951/0.964/0.957 | 0.439/0.870 | 0.364/0.857 | 0.065/0.494 | 0.471 |
| pred_preddepth (standalone theta=1) | 0.931/0.948/0.940 | 0.477/0.825 | 0.380/0.815 | 0.102/0.367 | 0.516 |
- standalone theta=1 over-call: total 15,977 calls vs 13,942 real (all 1.146x; novel 2.148x, dORF 3.589x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.930/0.916/0.923 | 0.717/0.239 | 0.579/0.273 | 0.270/0.133 | 0.312 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 19,529/15,415 | 0.940 | 0.417/0.848/0.559 | 0.323/0.795/0.459 | 0.094/0.380/0.151 | 0.379/0.711/0.494 |
| 10 | 18,119/14,988 | 0.940 | 0.457/0.843/0.593 | 0.342/0.794/0.478 | 0.102/0.388/0.162 | 0.406/0.702/0.515 |
| 15 | 16,968/14,463 | 0.940 | 0.470/0.832/0.601 | 0.365/0.815/0.504 | 0.103/0.378/0.162 | 0.415/0.687/0.517 |
| 20 | 15,977/13,942 | 0.940 | 0.477/0.825/0.604 | 0.380/0.815/0.518 | 0.102/0.367/0.160 | 0.417/0.674/0.516 |
| 30 | 14,448/13,110 | 0.940 | 0.497/0.812/0.617 | 0.406/0.836/0.546 | 0.104/0.356/0.161 | 0.427/0.653/0.516 |
| 40 | 13,483/12,492 | 0.940 | 0.489/0.794/0.606 | 0.418/0.856/0.561 | 0.094/0.336/0.147 | 0.419/0.637/0.506 |

### onehot attn4  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5236 (e12)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.961/0.976/0.968 | 0.493/0.895 | 0.402/0.881 | 0.095/0.567 | 0.535 |
| pred_preddepth (standalone theta=1) | 0.934/0.948/0.941 | 0.573/0.823 | 0.437/0.802 | 0.166/0.350 | 0.593 |
- standalone theta=1 over-call: total 15,196 calls vs 13,942 real (all 1.09x; novel 1.835x, dORF 2.106x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.89) | 0.934/0.886/0.909 | 0.806/0.256 | 0.686/0.293 | 0.460/0.094 | 0.352 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,213/15,415 | 0.941 | 0.495/0.850/0.626 | 0.376/0.792/0.510 | 0.159/0.375/0.223 | 0.462/0.735/0.568 |
| 10 | 17,049/14,988 | 0.941 | 0.539/0.843/0.658 | 0.398/0.796/0.531 | 0.174/0.383/0.239 | 0.494/0.727/0.588 |
| 15 | 16,052/14,463 | 0.941 | 0.557/0.827/0.665 | 0.421/0.805/0.553 | 0.171/0.367/0.233 | 0.506/0.711/0.591 |
| 20 | 15,196/13,942 | 0.941 | 0.573/0.823/0.676 | 0.437/0.802/0.566 | 0.166/0.350/0.225 | 0.514/0.702/0.593 |
| 30 | 13,876/13,110 | 0.941 | 0.601/0.800/0.686 | 0.470/0.815/0.596 | 0.167/0.336/0.223 | 0.530/0.678/0.595 |
| 40 | 13,035/12,492 | 0.941 | 0.612/0.796/0.692 | 0.486/0.810/0.608 | 0.150/0.327/0.206 | 0.532/0.664/0.590 |

### rinalmo  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_rinalmo_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5000 (e6)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.951/0.962/0.957 | 0.408/0.914 | 0.369/0.903 | 0.061/0.533 | 0.460 |
| pred_preddepth (standalone theta=1) | 0.929/0.943/0.936 | 0.450/0.863 | 0.379/0.839 | 0.101/0.394 | 0.511 |
- standalone theta=1 over-call: total 16,472 calls vs 13,942 real (all 1.181x; novel 2.214x, dORF 3.889x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.928/0.916/0.922 | 0.643/0.245 | 0.566/0.247 | 0.244/0.106 | 0.307 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 20,384/15,415 | 0.936 | 0.390/0.876/0.540 | 0.327/0.820/0.467 | 0.089/0.412/0.146 | 0.357/0.740/0.482 |
| 10 | 18,841/14,988 | 0.936 | 0.428/0.872/0.575 | 0.340/0.819/0.481 | 0.096/0.411/0.156 | 0.383/0.731/0.502 |
| 15 | 17,556/14,463 | 0.936 | 0.442/0.864/0.585 | 0.363/0.840/0.507 | 0.101/0.408/0.161 | 0.393/0.720/0.509 |
| 20 | 16,472/13,942 | 0.936 | 0.450/0.863/0.592 | 0.379/0.839/0.522 | 0.101/0.394/0.161 | 0.399/0.712/0.511 |
| 30 | 14,734/13,110 | 0.936 | 0.471/0.845/0.605 | 0.406/0.863/0.552 | 0.110/0.397/0.172 | 0.412/0.688/0.515 |
| 40 | 13,649/12,492 | 0.936 | 0.468/0.833/0.600 | 0.424/0.866/0.570 | 0.091/0.355/0.145 | 0.406/0.664/0.504 |

### T-cell GSE155087 (mouse held-out)  (mouse-tcell universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout/mouse_gse155087_tcell/onehot`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 11,417 (CDS 9,159, uORF 1,107, novel 275, dORF 211)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.976/0.978/0.977 | 0.619/0.636 | 0.483/0.636 | 0.164/0.379 | 0.489 |
| pred_preddepth (standalone theta=1) | 0.961/0.970/0.966 | 0.553/0.758 | 0.412/0.782 | 0.133/0.408 | 0.495 |
- standalone theta=1 over-call: total 12,336 calls vs 11,417 real (all 1.08x; novel 1.898x, dORF 3.071x)
| Poisson CDS-anchored (th=0.02, CDSrec 0.92) | 0.959/0.919/0.939 | 0.733/0.070 | 0.592/0.269 | 0.257/0.043 | 0.145 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 14,890/12,792 | 0.966 | 0.491/0.770/0.600 | 0.346/0.750/0.473 | 0.106/0.394/0.167 | 0.406/0.631/0.494 |
| 10 | 13,869/12,372 | 0.966 | 0.533/0.764/0.628 | 0.367/0.753/0.493 | 0.116/0.392/0.179 | 0.428/0.616/0.505 |
| 15 | 13,057/11,887 | 0.966 | 0.548/0.761/0.637 | 0.391/0.763/0.517 | 0.123/0.390/0.187 | 0.432/0.603/0.503 |
| 20 | 12,336/11,417 | 0.966 | 0.553/0.758/0.639 | 0.412/0.782/0.539 | 0.133/0.408/0.200 | 0.429/0.587/0.495 |
| 30 | 11,346/10,726 | 0.966 | 0.559/0.756/0.643 | 0.457/0.817/0.586 | 0.145/0.418/0.215 | 0.421/0.564/0.482 |
| 40 | 10,643/10,192 | 0.966 | 0.549/0.762/0.638 | 0.474/0.834/0.604 | 0.143/0.436/0.216 | 0.400/0.541/0.460 |

### onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.6099 (e21)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.945/0.962/0.953 | 0.502/0.836 | 0.358/0.850 | 0.081/0.559 | 0.517 |
| pred_preddepth (standalone theta=1) | 0.913/0.928/0.920 | 0.576/0.725 | 0.381/0.767 | 0.124/0.397 | 0.547 |
- standalone theta=1 over-call: total 19,072 calls vs 17,284 real (all 1.103x; novel 2.012x, dORF 3.201x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.90) | 0.911/0.898/0.904 | 0.766/0.342 | 0.591/0.397 | 0.270/0.190 | 0.438 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 22,442/18,934 | 0.920 | 0.510/0.767/0.613 | 0.337/0.757/0.466 | 0.107/0.412/0.171 | 0.435/0.684/0.532 |
| 10 | 21,204/18,469 | 0.920 | 0.548/0.755/0.635 | 0.351/0.758/0.480 | 0.117/0.418/0.183 | 0.456/0.674/0.544 |
| 15 | 20,086/17,891 | 0.920 | 0.562/0.736/0.637 | 0.368/0.763/0.496 | 0.121/0.415/0.188 | 0.464/0.660/0.545 |
| 20 | 19,072/17,284 | 0.920 | 0.576/0.725/0.642 | 0.381/0.767/0.509 | 0.124/0.397/0.189 | 0.471/0.652/0.547 |
| 30 | 17,480/16,240 | 0.920 | 0.605/0.711/0.653 | 0.400/0.786/0.530 | 0.126/0.372/0.189 | 0.482/0.645/0.552 |
| 40 | 16,377/15,418 | 0.920 | 0.605/0.694/0.646 | 0.408/0.787/0.537 | 0.117/0.374/0.179 | 0.476/0.635/0.544 |

### mamba onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.6020 (e21)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.943/0.961/0.952 | 0.468/0.848 | 0.355/0.831 | 0.071/0.592 | 0.497 |
| pred_preddepth (standalone theta=1) | 0.911/0.926/0.918 | 0.554/0.741 | 0.398/0.767 | 0.125/0.386 | 0.548 |
- standalone theta=1 over-call: total 19,044 calls vs 17,284 real (all 1.102x; novel 1.927x, dORF 3.078x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.89) | 0.909/0.892/0.900 | 0.725/0.348 | 0.588/0.360 | 0.221/0.168 | 0.424 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 22,333/18,934 | 0.918 | 0.499/0.775/0.607 | 0.351/0.758/0.480 | 0.107/0.398/0.169 | 0.438/0.682/0.534 |
| 10 | 21,160/18,469 | 0.918 | 0.536/0.769/0.632 | 0.362/0.757/0.490 | 0.115/0.404/0.178 | 0.459/0.675/0.546 |
| 15 | 20,046/17,891 | 0.918 | 0.546/0.750/0.632 | 0.382/0.764/0.509 | 0.119/0.390/0.182 | 0.467/0.660/0.547 |
| 20 | 19,044/17,284 | 0.918 | 0.554/0.741/0.634 | 0.398/0.767/0.524 | 0.125/0.386/0.189 | 0.473/0.652/0.548 |
| 30 | 17,433/16,240 | 0.919 | 0.582/0.733/0.649 | 0.417/0.789/0.545 | 0.142/0.393/0.209 | 0.489/0.647/0.557 |
| 40 | 16,307/15,418 | 0.919 | 0.591/0.725/0.651 | 0.422/0.790/0.550 | 0.135/0.383/0.200 | 0.487/0.636/0.552 |

### onehot attn2 seed1  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed1`  |  val_pearson 0.5216 (e6)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.960/0.974/0.967 | 0.421/0.899 | 0.373/0.877 | 0.063/0.544 | 0.466 |
| pred_preddepth (standalone theta=1) | 0.940/0.959/0.949 | 0.462/0.871 | 0.384/0.835 | 0.089/0.417 | 0.513 |
- standalone theta=1 over-call: total 16,544 calls vs 13,942 real (all 1.187x; novel 2.174x, dORF 4.694x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.93) | 0.938/0.935/0.936 | 0.728/0.315 | 0.585/0.326 | 0.216/0.122 | 0.372 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 20,512/15,415 | 0.949 | 0.403/0.893/0.555 | 0.324/0.818/0.464 | 0.080/0.435/0.135 | 0.360/0.754/0.487 |
| 10 | 18,990/14,988 | 0.949 | 0.439/0.885/0.587 | 0.342/0.821/0.483 | 0.085/0.440/0.142 | 0.384/0.743/0.506 |
| 15 | 17,675/14,463 | 0.949 | 0.453/0.873/0.596 | 0.366/0.834/0.509 | 0.086/0.429/0.144 | 0.393/0.728/0.510 |
| 20 | 16,544/13,942 | 0.949 | 0.462/0.871/0.604 | 0.384/0.835/0.526 | 0.089/0.417/0.146 | 0.399/0.717/0.513 |
| 30 | 14,850/13,110 | 0.949 | 0.478/0.857/0.614 | 0.413/0.858/0.557 | 0.095/0.425/0.155 | 0.408/0.696/0.515 |
| 40 | 13,741/12,492 | 0.950 | 0.474/0.842/0.606 | 0.424/0.866/0.570 | 0.082/0.383/0.136 | 0.399/0.667/0.499 |

### attn4 onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.6095 (e23)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.943/0.960/0.951 | 0.489/0.857 | 0.372/0.844 | 0.094/0.514 | 0.530 |
| pred_preddepth (standalone theta=1) | 0.916/0.932/0.924 | 0.572/0.765 | 0.392/0.787 | 0.155/0.369 | 0.566 |
- standalone theta=1 over-call: total 19,103 calls vs 17,284 real (all 1.105x; novel 2.01x, dORF 2.38x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.90) | 0.915/0.901/0.908 | 0.751/0.360 | 0.565/0.380 | 0.276/0.151 | 0.442 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 22,622/18,934 | 0.924 | 0.505/0.806/0.621 | 0.342/0.782/0.475 | 0.137/0.403/0.204 | 0.443/0.711/0.546 |
| 10 | 21,321/18,469 | 0.924 | 0.545/0.795/0.647 | 0.358/0.786/0.492 | 0.150/0.404/0.219 | 0.468/0.702/0.562 |
| 15 | 20,171/17,891 | 0.924 | 0.558/0.778/0.650 | 0.376/0.788/0.509 | 0.156/0.390/0.223 | 0.478/0.689/0.564 |
| 20 | 19,103/17,284 | 0.924 | 0.572/0.765/0.655 | 0.392/0.787/0.523 | 0.155/0.369/0.218 | 0.486/0.677/0.566 |
| 30 | 17,441/16,240 | 0.924 | 0.602/0.746/0.667 | 0.413/0.796/0.544 | 0.170/0.372/0.234 | 0.501/0.663/0.571 |
| 40 | 16,294/15,418 | 0.924 | 0.605/0.729/0.661 | 0.420/0.790/0.549 | 0.156/0.346/0.215 | 0.497/0.646/0.562 |

### onehot attn2 seed2  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed2`  |  val_pearson 0.5248 (e12)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.959/0.974/0.966 | 0.511/0.882 | 0.390/0.894 | 0.091/0.511 | 0.534 |
| pred_preddepth (standalone theta=1) | 0.933/0.947/0.940 | 0.595/0.764 | 0.441/0.806 | 0.184/0.350 | 0.581 |
- standalone theta=1 over-call: total 14,845 calls vs 13,942 real (all 1.065x; novel 1.828x, dORF 1.9x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.88) | 0.931/0.877/0.904 | 0.766/0.211 | 0.676/0.216 | 0.447/0.094 | 0.294 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 17,599/15,415 | 0.940 | 0.522/0.807/0.634 | 0.382/0.790/0.515 | 0.168/0.366/0.231 | 0.478/0.694/0.566 |
| 10 | 16,541/14,988 | 0.940 | 0.566/0.793/0.660 | 0.399/0.792/0.531 | 0.181/0.373/0.244 | 0.506/0.681/0.581 |
| 15 | 15,626/14,463 | 0.940 | 0.581/0.775/0.664 | 0.423/0.803/0.554 | 0.180/0.357/0.239 | 0.517/0.664/0.582 |
| 20 | 14,845/13,942 | 0.940 | 0.595/0.764/0.669 | 0.441/0.806/0.570 | 0.184/0.350/0.241 | 0.523/0.653/0.581 |
| 30 | 13,666/13,110 | 0.940 | 0.619/0.741/0.674 | 0.466/0.825/0.596 | 0.190/0.343/0.244 | 0.535/0.633/0.580 |
| 40 | 12,913/12,492 | 0.940 | 0.631/0.740/0.681 | 0.474/0.820/0.601 | 0.182/0.364/0.243 | 0.534/0.625/0.576 |

### mamba4 onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.6190 (e23)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.946/0.963/0.955 | 0.517/0.834 | 0.387/0.833 | 0.105/0.514 | 0.548 |
| pred_preddepth (standalone theta=1) | 0.912/0.929/0.921 | 0.590/0.754 | 0.397/0.771 | 0.156/0.397 | 0.570 |
- standalone theta=1 over-call: total 18,922 calls vs 17,284 real (all 1.095x; novel 1.942x, dORF 2.536x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.90) | 0.913/0.901/0.907 | 0.756/0.370 | 0.552/0.388 | 0.305/0.162 | 0.443 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 22,180/18,934 | 0.921 | 0.525/0.788/0.631 | 0.348/0.766/0.478 | 0.137/0.426/0.207 | 0.457/0.698/0.553 |
| 10 | 21,025/18,469 | 0.921 | 0.562/0.781/0.654 | 0.363/0.767/0.493 | 0.148/0.433/0.220 | 0.479/0.690/0.565 |
| 15 | 19,918/17,891 | 0.921 | 0.581/0.765/0.660 | 0.382/0.774/0.511 | 0.152/0.415/0.222 | 0.490/0.678/0.569 |
| 20 | 18,922/17,284 | 0.921 | 0.590/0.754/0.662 | 0.397/0.771/0.524 | 0.156/0.397/0.224 | 0.497/0.668/0.570 |
| 30 | 17,398/16,240 | 0.921 | 0.615/0.740/0.671 | 0.418/0.792/0.547 | 0.158/0.379/0.223 | 0.505/0.659/0.572 |
| 40 | 16,288/15,418 | 0.921 | 0.623/0.729/0.672 | 0.427/0.796/0.556 | 0.153/0.383/0.219 | 0.503/0.649/0.567 |

### onehot attn4 seed2  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed2`  |  val_pearson 0.5300 (e10)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.959/0.975/0.967 | 0.521/0.860 | 0.381/0.868 | 0.079/0.561 | 0.515 |
| pred_preddepth (standalone theta=1) | 0.936/0.949/0.942 | 0.594/0.730 | 0.430/0.777 | 0.127/0.411 | 0.549 |
- standalone theta=1 over-call: total 14,867 calls vs 13,942 real (all 1.066x; novel 1.808x, dORF 3.233x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.88) | 0.934/0.879/0.906 | 0.765/0.185 | 0.656/0.260 | 0.246/0.094 | 0.284 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 17,608/15,415 | 0.942 | 0.523/0.769/0.623 | 0.364/0.765/0.494 | 0.119/0.440/0.188 | 0.454/0.660/0.538 |
| 10 | 16,542/14,988 | 0.942 | 0.567/0.756/0.648 | 0.389/0.766/0.516 | 0.129/0.445/0.200 | 0.481/0.648/0.552 |
| 15 | 15,643/14,463 | 0.942 | 0.581/0.740/0.651 | 0.409/0.776/0.536 | 0.133/0.439/0.204 | 0.489/0.631/0.551 |
| 20 | 14,867/13,942 | 0.942 | 0.594/0.730/0.655 | 0.430/0.777/0.554 | 0.127/0.411/0.194 | 0.493/0.619/0.549 |
| 30 | 13,716/13,110 | 0.942 | 0.613/0.717/0.661 | 0.461/0.798/0.584 | 0.133/0.411/0.201 | 0.499/0.602/0.546 |
| 40 | 12,965/12,492 | 0.943 | 0.626/0.707/0.664 | 0.469/0.803/0.593 | 0.121/0.421/0.188 | 0.488/0.589/0.534 |

### attn6 onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn6_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5808 (e9)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.950/0.963/0.956 | 0.472/0.841 | 0.325/0.850 | 0.075/0.520 | 0.492 |
| pred_preddepth (standalone theta=1) | 0.911/0.928/0.920 | 0.564/0.708 | 0.372/0.751 | 0.148/0.318 | 0.538 |
- standalone theta=1 over-call: total 18,782 calls vs 17,284 real (all 1.087x; novel 2.018x, dORF 2.151x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.91) | 0.911/0.912/0.911 | 0.731/0.303 | 0.560/0.351 | 0.250/0.123 | 0.394 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 22,214/18,934 | 0.920 | 0.491/0.754/0.595 | 0.322/0.752/0.451 | 0.132/0.343/0.191 | 0.431/0.661/0.522 |
| 10 | 20,903/18,469 | 0.920 | 0.534/0.743/0.622 | 0.339/0.750/0.467 | 0.145/0.346/0.204 | 0.458/0.650/0.538 |
| 15 | 19,770/17,891 | 0.920 | 0.552/0.723/0.626 | 0.359/0.750/0.485 | 0.145/0.333/0.202 | 0.469/0.634/0.539 |
| 20 | 18,782/17,284 | 0.920 | 0.564/0.708/0.628 | 0.372/0.751/0.498 | 0.148/0.318/0.202 | 0.475/0.622/0.538 |
| 30 | 17,256/16,240 | 0.920 | 0.584/0.681/0.629 | 0.392/0.770/0.519 | 0.150/0.310/0.203 | 0.482/0.607/0.538 |
| 40 | 16,191/15,418 | 0.920 | 0.587/0.669/0.625 | 0.396/0.771/0.523 | 0.140/0.290/0.189 | 0.479/0.596/0.531 |

### onehot attn4 seed1  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_attn4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed1`  |  val_pearson 0.5239 (e6)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.959/0.971/0.965 | 0.433/0.904 | 0.372/0.888 | 0.071/0.567 | 0.483 |
| pred_preddepth (standalone theta=1) | 0.933/0.947/0.940 | 0.517/0.817 | 0.418/0.797 | 0.145/0.306 | 0.557 |
- standalone theta=1 over-call: total 15,364 calls vs 13,942 real (all 1.102x; novel 1.907x, dORF 2.1x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.89) | 0.931/0.887/0.908 | 0.763/0.223 | 0.614/0.220 | 0.412/0.078 | 0.290 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,560/15,415 | 0.940 | 0.451/0.844/0.587 | 0.362/0.782/0.495 | 0.132/0.315/0.185 | 0.427/0.712/0.533 |
| 10 | 17,329/14,988 | 0.940 | 0.490/0.836/0.618 | 0.379/0.780/0.510 | 0.141/0.321/0.196 | 0.455/0.700/0.551 |
| 15 | 16,253/14,463 | 0.940 | 0.506/0.819/0.626 | 0.402/0.795/0.533 | 0.141/0.311/0.194 | 0.467/0.683/0.554 |
| 20 | 15,364/13,942 | 0.940 | 0.517/0.817/0.633 | 0.418/0.797/0.548 | 0.145/0.306/0.197 | 0.474/0.674/0.557 |
| 30 | 13,989/13,110 | 0.940 | 0.541/0.800/0.645 | 0.450/0.823/0.582 | 0.155/0.315/0.208 | 0.491/0.655/0.561 |
| 40 | 13,116/12,492 | 0.940 | 0.551/0.793/0.650 | 0.460/0.826/0.591 | 0.133/0.299/0.184 | 0.487/0.634/0.551 |

### onehot mamba seed1  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed1`  |  val_pearson 0.5140 (e4)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.962/0.973/0.967 | 0.428/0.892 | 0.379/0.874 | 0.090/0.444 | 0.501 |
| pred_preddepth (standalone theta=1) | 0.933/0.947/0.940 | 0.498/0.824 | 0.404/0.815 | 0.145/0.244 | 0.546 |
- standalone theta=1 over-call: total 15,351 calls vs 13,942 real (all 1.101x; novel 2.018x, dORF 1.689x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.91) | 0.930/0.911/0.920 | 0.735/0.179 | 0.571/0.267 | 0.333/0.089 | 0.258 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,597/15,415 | 0.940 | 0.445/0.848/0.584 | 0.344/0.803/0.481 | 0.127/0.282/0.175 | 0.419/0.703/0.525 |
| 10 | 17,353/14,988 | 0.940 | 0.480/0.840/0.611 | 0.363/0.807/0.501 | 0.141/0.282/0.188 | 0.447/0.691/0.543 |
| 15 | 16,280/14,463 | 0.940 | 0.490/0.829/0.616 | 0.387/0.817/0.525 | 0.146/0.260/0.187 | 0.458/0.673/0.545 |
| 20 | 15,351/13,942 | 0.940 | 0.498/0.824/0.621 | 0.404/0.815/0.540 | 0.145/0.244/0.182 | 0.466/0.659/0.546 |
| 30 | 13,969/13,110 | 0.940 | 0.512/0.814/0.629 | 0.435/0.839/0.573 | 0.158/0.240/0.191 | 0.481/0.636/0.548 |
| 40 | 13,094/12,492 | 0.940 | 0.500/0.796/0.614 | 0.449/0.849/0.587 | 0.141/0.215/0.170 | 0.476/0.613/0.536 |

### onehot mamba seed2  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed2`  |  val_pearson 0.5199 (e4)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.962/0.974/0.968 | 0.446/0.873 | 0.383/0.866 | 0.088/0.461 | 0.499 |
| pred_preddepth (standalone theta=1) | 0.934/0.949/0.941 | 0.508/0.806 | 0.396/0.819 | 0.143/0.289 | 0.539 |
- standalone theta=1 over-call: total 15,274 calls vs 13,942 real (all 1.096x; novel 2.07x, dORF 2.022x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.932/0.918/0.925 | 0.725/0.192 | 0.592/0.291 | 0.294/0.083 | 0.273 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,477/15,415 | 0.941 | 0.452/0.833/0.586 | 0.341/0.816/0.481 | 0.122/0.319/0.177 | 0.419/0.690/0.521 |
| 10 | 17,234/14,988 | 0.941 | 0.488/0.822/0.613 | 0.361/0.817/0.501 | 0.135/0.321/0.190 | 0.447/0.676/0.538 |
| 15 | 16,202/14,463 | 0.941 | 0.498/0.810/0.617 | 0.380/0.824/0.520 | 0.139/0.306/0.191 | 0.455/0.658/0.538 |
| 20 | 15,274/13,942 | 0.941 | 0.508/0.806/0.623 | 0.396/0.819/0.534 | 0.143/0.289/0.191 | 0.464/0.644/0.539 |
| 30 | 13,924/13,110 | 0.941 | 0.522/0.789/0.628 | 0.431/0.847/0.572 | 0.146/0.274/0.191 | 0.475/0.617/0.537 |
| 40 | 13,082/12,492 | 0.941 | 0.506/0.766/0.609 | 0.445/0.859/0.586 | 0.130/0.252/0.172 | 0.467/0.594/0.523 |

### onehot mamba4 seed2  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed2`  |  val_pearson 0.5086 (e4)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.961/0.974/0.967 | 0.459/0.889 | 0.379/0.885 | 0.076/0.528 | 0.497 |
| pred_preddepth (standalone theta=1) | 0.939/0.953/0.946 | 0.526/0.835 | 0.393/0.824 | 0.122/0.356 | 0.550 |
- standalone theta=1 over-call: total 15,503 calls vs 13,942 real (all 1.112x; novel 2.095x, dORF 2.906x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.936/0.922/0.929 | 0.775/0.222 | 0.564/0.302 | 0.238/0.111 | 0.305 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,917/15,415 | 0.946 | 0.454/0.857/0.593 | 0.342/0.818/0.482 | 0.098/0.384/0.157 | 0.410/0.716/0.522 |
| 10 | 17,558/14,988 | 0.946 | 0.497/0.850/0.627 | 0.363/0.823/0.504 | 0.114/0.388/0.176 | 0.443/0.705/0.544 |
| 15 | 16,461/14,463 | 0.946 | 0.513/0.839/0.637 | 0.382/0.830/0.523 | 0.117/0.372/0.178 | 0.454/0.690/0.548 |
| 20 | 15,503/13,942 | 0.946 | 0.526/0.835/0.645 | 0.393/0.824/0.532 | 0.122/0.356/0.182 | 0.463/0.677/0.550 |
| 30 | 14,067/13,110 | 0.946 | 0.540/0.818/0.651 | 0.424/0.836/0.562 | 0.137/0.343/0.196 | 0.477/0.650/0.550 |
| 40 | 13,155/12,492 | 0.946 | 0.538/0.800/0.644 | 0.443/0.856/0.584 | 0.126/0.318/0.180 | 0.477/0.630/0.543 |

### onehot mamba4 seed1  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes_seed1`  |  val_pearson 0.5103 (e3)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.959/0.969/0.964 | 0.407/0.907 | 0.351/0.872 | 0.055/0.506 | 0.451 |
| pred_preddepth (standalone theta=1) | 0.942/0.957/0.950 | 0.477/0.862 | 0.348/0.863 | 0.152/0.406 | 0.530 |
- standalone theta=1 over-call: total 16,048 calls vs 13,942 real (all 1.151x; novel 2.478x, dORF 2.661x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.936/0.923/0.930 | 0.704/0.185 | 0.562/0.311 | 0.275/0.078 | 0.275 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 19,895/15,415 | 0.950 | 0.412/0.881/0.561 | 0.289/0.854/0.432 | 0.123/0.412/0.189 | 0.376/0.739/0.498 |
| 10 | 18,370/14,988 | 0.950 | 0.453/0.875/0.597 | 0.311/0.860/0.457 | 0.138/0.416/0.208 | 0.406/0.729/0.522 |
| 15 | 17,133/14,463 | 0.950 | 0.467/0.865/0.606 | 0.330/0.867/0.478 | 0.145/0.408/0.214 | 0.417/0.714/0.527 |
| 20 | 16,048/13,942 | 0.950 | 0.477/0.862/0.614 | 0.348/0.863/0.496 | 0.152/0.406/0.222 | 0.426/0.702/0.530 |
| 30 | 14,421/13,110 | 0.950 | 0.494/0.846/0.624 | 0.383/0.882/0.534 | 0.165/0.384/0.231 | 0.442/0.675/0.534 |
| 40 | 13,421/12,492 | 0.950 | 0.486/0.833/0.614 | 0.397/0.888/0.549 | 0.158/0.355/0.218 | 0.440/0.656/0.526 |

### onehot mamba4  (fibroblast universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba4_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5087 (e6)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 13,942 (CDS 10,947, uORF 1,461, novel 454, dORF 180)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.961/0.975/0.968 | 0.435/0.916 | 0.389/0.881 | 0.085/0.611 | 0.499 |
| pred_preddepth (standalone theta=1) | 0.933/0.949/0.941 | 0.506/0.870 | 0.445/0.828 | 0.154/0.350 | 0.570 |
- standalone theta=1 over-call: total 15,590 calls vs 13,942 real (all 1.118x; novel 1.861x, dORF 2.267x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.91) | 0.931/0.907/0.919 | 0.748/0.255 | 0.621/0.238 | 0.357/0.056 | 0.326 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 19,040/15,415 | 0.941 | 0.440/0.886/0.588 | 0.383/0.803/0.519 | 0.132/0.380/0.196 | 0.422/0.747/0.539 |
| 10 | 17,709/14,988 | 0.941 | 0.477/0.880/0.619 | 0.402/0.805/0.536 | 0.147/0.378/0.212 | 0.452/0.736/0.560 |
| 15 | 16,579/14,463 | 0.941 | 0.491/0.868/0.627 | 0.426/0.828/0.562 | 0.152/0.367/0.215 | 0.464/0.720/0.564 |
| 20 | 15,590/13,942 | 0.941 | 0.506/0.870/0.640 | 0.445/0.828/0.579 | 0.154/0.350/0.214 | 0.476/0.710/0.570 |
| 30 | 14,171/13,110 | 0.941 | 0.526/0.870/0.656 | 0.473/0.863/0.611 | 0.173/0.356/0.233 | 0.496/0.698/0.580 |
| 40 | 13,224/12,492 | 0.941 | 0.535/0.865/0.661 | 0.485/0.875/0.625 | 0.156/0.327/0.211 | 0.500/0.680/0.577 |

### mamba ds32 onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba_ds32_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5982 (e21)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.946/0.964/0.955 | 0.481/0.841 | 0.368/0.834 | 0.070/0.598 | 0.502 |
| pred_preddepth (standalone theta=1) | 0.907/0.924/0.916 | 0.559/0.747 | 0.392/0.771 | 0.118/0.458 | 0.543 |
- standalone theta=1 over-call: total 19,220 calls vs 17,284 real (all 1.112x; novel 1.969x, dORF 3.883x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.89) | 0.907/0.893/0.900 | 0.740/0.370 | 0.564/0.361 | 0.197/0.168 | 0.438 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 22,623/18,934 | 0.916 | 0.499/0.782/0.609 | 0.347/0.760/0.477 | 0.105/0.495/0.173 | 0.427/0.686/0.527 |
| 10 | 21,383/18,469 | 0.916 | 0.535/0.776/0.634 | 0.360/0.760/0.489 | 0.114/0.495/0.185 | 0.449/0.678/0.540 |
| 15 | 20,244/17,891 | 0.916 | 0.550/0.760/0.638 | 0.378/0.766/0.506 | 0.115/0.472/0.185 | 0.458/0.666/0.543 |
| 20 | 19,220/17,284 | 0.916 | 0.559/0.747/0.639 | 0.392/0.771/0.520 | 0.118/0.458/0.188 | 0.463/0.656/0.543 |
| 30 | 17,600/16,240 | 0.916 | 0.580/0.732/0.647 | 0.413/0.797/0.544 | 0.128/0.455/0.199 | 0.473/0.649/0.547 |
| 40 | 16,434/15,418 | 0.916 | 0.590/0.727/0.652 | 0.420/0.799/0.551 | 0.134/0.477/0.209 | 0.475/0.644/0.547 |

### mamba6 onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba6_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.6251 (e20)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.945/0.962/0.953 | 0.528/0.818 | 0.398/0.815 | 0.106/0.520 | 0.551 |
| pred_preddepth (standalone theta=1) | 0.911/0.926/0.919 | 0.605/0.731 | 0.423/0.763 | 0.151/0.357 | 0.577 |
- standalone theta=1 over-call: total 18,626 calls vs 17,284 real (all 1.078x; novel 1.805x, dORF 2.363x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.89) | 0.909/0.894/0.902 | 0.778/0.362 | 0.609/0.380 | 0.263/0.140 | 0.446 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 21,653/18,934 | 0.919 | 0.539/0.765/0.632 | 0.375/0.750/0.500 | 0.135/0.389/0.200 | 0.475/0.683/0.561 |
| 10 | 20,549/18,469 | 0.919 | 0.579/0.757/0.656 | 0.389/0.750/0.512 | 0.147/0.394/0.215 | 0.499/0.675/0.574 |
| 15 | 19,531/17,891 | 0.919 | 0.598/0.744/0.663 | 0.409/0.759/0.532 | 0.147/0.374/0.212 | 0.511/0.666/0.578 |
| 20 | 18,626/17,284 | 0.919 | 0.605/0.731/0.662 | 0.423/0.763/0.544 | 0.151/0.357/0.213 | 0.515/0.656/0.577 |
| 30 | 17,165/16,240 | 0.919 | 0.628/0.719/0.670 | 0.439/0.773/0.560 | 0.157/0.338/0.214 | 0.525/0.646/0.579 |
| 40 | 16,106/15,418 | 0.919 | 0.638/0.703/0.669 | 0.445/0.774/0.565 | 0.139/0.308/0.192 | 0.522/0.632/0.572 |

### mamba ds64 onehot union  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba_ds64_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`  |  val_pearson 0.5904 (e16)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.948/0.964/0.956 | 0.445/0.854 | 0.351/0.845 | 0.080/0.542 | 0.501 |
| pred_preddepth (standalone theta=1) | 0.908/0.923/0.915 | 0.512/0.769 | 0.381/0.754 | 0.112/0.386 | 0.534 |
- standalone theta=1 over-call: total 19,538 calls vs 17,284 real (all 1.13x; novel 1.982x, dORF 3.436x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.89) | 0.909/0.893/0.901 | 0.725/0.399 | 0.541/0.361 | 0.218/0.184 | 0.446 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 23,303/18,934 | 0.915 | 0.454/0.799/0.579 | 0.331/0.746/0.458 | 0.100/0.421/0.161 | 0.404/0.697/0.511 |
| 10 | 21,923/18,469 | 0.915 | 0.487/0.791/0.603 | 0.346/0.745/0.472 | 0.110/0.428/0.175 | 0.426/0.689/0.526 |
| 15 | 20,654/17,891 | 0.915 | 0.504/0.781/0.613 | 0.368/0.751/0.494 | 0.111/0.405/0.174 | 0.439/0.678/0.533 |
| 20 | 19,538/17,284 | 0.915 | 0.512/0.769/0.615 | 0.381/0.754/0.506 | 0.112/0.386/0.174 | 0.445/0.668/0.534 |
| 30 | 17,771/16,240 | 0.915 | 0.531/0.753/0.623 | 0.407/0.781/0.535 | 0.126/0.386/0.190 | 0.459/0.658/0.541 |
| 40 | 16,493/15,418 | 0.915 | 0.545/0.742/0.628 | 0.415/0.776/0.541 | 0.123/0.383/0.187 | 0.464/0.644/0.539 |

### mamba4 onehot mm10cov armB excl  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm10cov_exclmm50_holdout_Hepatocytes`  |  val_pearson 0.6084 (e13)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.949/0.967/0.958 | 0.492/0.835 | 0.352/0.843 | 0.081/0.609 | 0.510 |
| pred_preddepth (standalone theta=1) | 0.915/0.933/0.924 | 0.524/0.791 | 0.345/0.797 | 0.110/0.503 | 0.527 |
- standalone theta=1 over-call: total 20,105 calls vs 17,284 real (all 1.163x; novel 2.307x, dORF 4.592x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.92) | 0.915/0.918/0.917 | 0.713/0.429 | 0.511/0.453 | 0.179/0.207 | 0.470 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 24,171/18,934 | 0.924 | 0.458/0.821/0.588 | 0.298/0.790/0.433 | 0.096/0.528/0.163 | 0.386/0.723/0.504 |
| 10 | 22,701/18,469 | 0.924 | 0.494/0.816/0.616 | 0.312/0.791/0.447 | 0.105/0.534/0.175 | 0.407/0.716/0.519 |
| 15 | 21,344/17,891 | 0.924 | 0.512/0.801/0.625 | 0.329/0.795/0.466 | 0.108/0.523/0.180 | 0.418/0.704/0.524 |
| 20 | 20,105/17,284 | 0.924 | 0.524/0.791/0.630 | 0.345/0.797/0.482 | 0.110/0.503/0.180 | 0.425/0.693/0.527 |
| 30 | 18,192/16,240 | 0.924 | 0.550/0.782/0.646 | 0.365/0.812/0.504 | 0.120/0.497/0.193 | 0.437/0.683/0.533 |
| 40 | 16,886/15,418 | 0.925 | 0.552/0.767/0.642 | 0.365/0.808/0.503 | 0.116/0.486/0.188 | 0.429/0.666/0.522 |

### mamba4 onehot mm10cov armA  (union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm10cov_holdout_Hepatocytes`  |  val_pearson 0.6155 (e26)
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 17,284 (CDS 13,223, uORF 1,977, novel 827, dORF 179)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.945/0.962/0.954 | 0.523/0.821 | 0.388/0.822 | 0.103/0.559 | 0.544 |
| pred_preddepth (standalone theta=1) | 0.911/0.927/0.919 | 0.598/0.743 | 0.396/0.753 | 0.155/0.419 | 0.569 |
- standalone theta=1 over-call: total 18,812 calls vs 17,284 real (all 1.088x; novel 1.902x, dORF 2.709x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.89) | 0.912/0.894/0.903 | 0.784/0.357 | 0.562/0.375 | 0.327/0.179 | 0.445 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 21,893/18,934 | 0.919 | 0.539/0.777/0.636 | 0.349/0.748/0.476 | 0.136/0.440/0.208 | 0.466/0.689/0.556 |
| 10 | 20,815/18,469 | 0.919 | 0.573/0.771/0.657 | 0.363/0.746/0.488 | 0.151/0.452/0.227 | 0.486/0.682/0.568 |
| 15 | 19,756/17,891 | 0.919 | 0.590/0.754/0.662 | 0.380/0.751/0.504 | 0.154/0.436/0.228 | 0.496/0.669/0.570 |
| 20 | 18,812/17,284 | 0.919 | 0.598/0.743/0.662 | 0.396/0.753/0.519 | 0.155/0.419/0.226 | 0.500/0.660/0.569 |
| 30 | 17,303/16,240 | 0.919 | 0.619/0.727/0.668 | 0.422/0.784/0.549 | 0.161/0.400/0.230 | 0.513/0.654/0.575 |
| 40 | 16,213/15,418 | 0.919 | 0.632/0.720/0.673 | 0.428/0.788/0.555 | 0.148/0.374/0.212 | 0.513/0.646/0.572 |
