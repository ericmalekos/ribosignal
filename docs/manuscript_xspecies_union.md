
### Janich liver (mamba4_union xspecies)  (xsp-mamba4_union universe)

> **SUPERSEDED (2026-07-30).** This Janich block was computed on the broken 2,434-tx universe (salmon run on untrimmed ARTseq reads, 0.03% mapping). Quarantined under `data/_quarantine/janich_stale_universe_2026-07-30/`. See the re-run further down this file, on the rebuilt 10,431-tx universe.
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_janich_liver`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 2,050 (CDS 1,656, uORF 143, novel 51, dORF 101)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.990/1.000/0.995 | 0.769/0.790 | 0.550/0.647 | 0.213/0.257 | 0.525 |
| pred_preddepth (standalone theta=1) | 0.990/0.979/0.984 | 0.800/0.364 | 0.522/0.471 | 0.227/0.050 | 0.370 |
- standalone theta=1 over-call: total 1,812 calls vs 2,050 real (all 0.884x; novel 0.902x, dORF 0.218x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.93) | 0.989/0.933/0.960 | 0.833/0.035 | 0.458/0.216 | 0.000/0.000 | 0.096 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 1,952/2,313 | 0.984 | 0.762/0.369/0.497 | 0.460/0.439/0.450 | 0.167/0.040/0.065 | 0.608/0.291/0.394 |
| 10 | 1,906/2,235 | 0.984 | 0.784/0.364/0.497 | 0.475/0.444/0.459 | 0.192/0.042/0.069 | 0.619/0.287/0.392 |
| 15 | 1,858/2,139 | 0.984 | 0.808/0.372/0.510 | 0.500/0.474/0.486 | 0.217/0.045/0.075 | 0.623/0.284/0.390 |
| 20 | 1,812/2,050 | 0.984 | 0.800/0.364/0.500 | 0.522/0.471/0.495 | 0.227/0.050/0.081 | 0.603/0.267/0.370 |
| 30 | 1,755/1,928 | 0.984 | 0.865/0.405/0.552 | 0.636/0.512/0.568 | 0.267/0.050/0.084 | 0.650/0.279/0.391 |
| 40 | 1,717/1,833 | 0.984 | 0.750/0.364/0.490 | 0.690/0.588/0.635 | 0.333/0.053/0.091 | 0.633/0.282/0.391 |

### GSE120762 NT (mamba4_union xspecies)  (xsp-mamba4_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_gse120762_nt`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 8,058 (CDS 7,654, uORF 223, novel 8, dORF 9)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.967/0.952/0.960 | 0.341/0.637 | 0.212/0.875 | 0.000/0.000 | 0.410 |
| pred_preddepth (standalone theta=1) | 0.949/0.950/0.949 | 0.181/0.888 | 0.167/0.875 | 0.014/0.333 | 0.284 |
- standalone theta=1 over-call: total 9,436 calls vs 8,058 real (all 1.171x; novel 5.25x, dORF 23.667x)
| Poisson CDS-anchored (th=0.02, CDSrec 0.88) | 0.948/0.877/0.911 | 0.443/0.193 | 0.375/0.750 | 0.000/0.000 | 0.225 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 11,081/8,158 | 0.949 | 0.110/0.886/0.196 | 0.100/0.875/0.179 | 0.009/0.300/0.018 | 0.115/0.779/0.200 |
| 10 | 10,445/8,144 | 0.949 | 0.136/0.888/0.235 | 0.111/0.875/0.197 | 0.011/0.300/0.021 | 0.137/0.780/0.234 |
| 15 | 9,893/8,106 | 0.949 | 0.160/0.886/0.271 | 0.127/0.875/0.222 | 0.012/0.300/0.024 | 0.157/0.774/0.261 |
| 20 | 9,436/8,058 | 0.949 | 0.181/0.888/0.301 | 0.167/0.875/0.280 | 0.014/0.333/0.027 | 0.174/0.765/0.284 |
| 30 | 8,845/7,970 | 0.949 | 0.201/0.863/0.326 | 0.259/0.875/0.400 | 0.019/0.375/0.036 | 0.194/0.725/0.306 |
| 40 | 8,445/7,886 | 0.949 | 0.233/0.851/0.366 | 0.350/1.000/0.518 | 0.017/0.333/0.032 | 0.212/0.712/0.327 |

### GSE120762 LPS (mamba4_union xspecies)  (xsp-mamba4_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_gse120762_lps`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 9,196 (CDS 8,144, uORF 574, novel 57, dORF 61)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.971/0.972/0.971 | 0.541/0.702 | 0.348/0.702 | 0.126/0.213 | 0.515 |
| pred_preddepth (standalone theta=1) | 0.950/0.953/0.952 | 0.424/0.828 | 0.259/0.772 | 0.062/0.279 | 0.469 |
- standalone theta=1 over-call: total 10,137 calls vs 9,196 real (all 1.102x; novel 2.982x, dORF 4.459x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.949/0.920/0.934 | 0.636/0.293 | 0.439/0.632 | 0.148/0.066 | 0.337 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 11,845/9,749 | 0.952 | 0.342/0.829/0.485 | 0.224/0.776/0.348 | 0.046/0.271/0.079 | 0.307/0.703/0.427 |
| 10 | 11,227/9,606 | 0.952 | 0.375/0.831/0.517 | 0.240/0.764/0.365 | 0.053/0.275/0.090 | 0.333/0.696/0.451 |
| 15 | 10,629/9,390 | 0.952 | 0.400/0.825/0.538 | 0.235/0.783/0.361 | 0.058/0.277/0.096 | 0.347/0.684/0.460 |
| 20 | 10,137/9,196 | 0.952 | 0.424/0.828/0.561 | 0.259/0.772/0.388 | 0.062/0.279/0.102 | 0.360/0.672/0.469 |
| 30 | 9,489/8,886 | 0.952 | 0.437/0.825/0.571 | 0.325/0.800/0.462 | 0.075/0.288/0.119 | 0.367/0.652/0.470 |
| 40 | 9,032/8,643 | 0.952 | 0.402/0.776/0.530 | 0.375/0.805/0.512 | 0.073/0.268/0.115 | 0.341/0.588/0.431 |

### T-cell GSE155087 (mamba4_union xspecies)  (xsp-mamba4_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_gse155087_tcell`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 11,417 (CDS 9,159, uORF 1,107, novel 275, dORF 211)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.976/0.984/0.980 | 0.664/0.617 | 0.488/0.654 | 0.254/0.341 | 0.524 |
| pred_preddepth (standalone theta=1) | 0.957/0.966/0.961 | 0.646/0.677 | 0.435/0.713 | 0.225/0.332 | 0.533 |
- standalone theta=1 over-call: total 11,544 calls vs 11,417 real (all 1.011x; novel 1.64x, dORF 1.474x)
| Poisson CDS-anchored (th=0.02, CDSrec 0.87) | 0.952/0.868/0.908 | 0.809/0.065 | 0.570/0.193 | 0.294/0.024 | 0.119 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 13,464/12,792 | 0.961 | 0.579/0.671/0.622 | 0.381/0.699/0.493 | 0.168/0.332/0.223 | 0.489/0.569/0.526 |
| 10 | 12,758/12,372 | 0.961 | 0.622/0.682/0.651 | 0.397/0.699/0.506 | 0.188/0.329/0.240 | 0.518/0.566/0.541 |
| 15 | 12,114/11,887 | 0.961 | 0.640/0.682/0.660 | 0.421/0.712/0.529 | 0.209/0.329/0.256 | 0.530/0.557/0.543 |
| 20 | 11,544/11,417 | 0.961 | 0.646/0.677/0.661 | 0.435/0.713/0.540 | 0.225/0.332/0.268 | 0.528/0.538/0.533 |
| 30 | 10,780/10,726 | 0.961 | 0.643/0.662/0.652 | 0.470/0.740/0.575 | 0.261/0.345/0.297 | 0.522/0.512/0.517 |
| 40 | 10,268/10,192 | 0.961 | 0.624/0.662/0.642 | 0.488/0.767/0.597 | 0.263/0.373/0.308 | 0.499/0.495/0.497 |

### Wang liver (mamba4_union xspecies)  (xsp-mamba4_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_wang_liver`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 14,088 (CDS 11,963, uORF 774, novel 673, dORF 162)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.960/0.963/0.962 | 0.458/0.441 | 0.409/0.683 | 0.145/0.290 | 0.408 |
| pred_preddepth (standalone theta=1) | 0.936/0.942/0.939 | 0.380/0.596 | 0.329/0.728 | 0.092/0.420 | 0.384 |
- standalone theta=1 over-call: total 15,952 calls vs 14,088 real (all 1.132x; novel 2.214x, dORF 4.574x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.89) | 0.933/0.889/0.911 | 0.597/0.143 | 0.422/0.459 | 0.239/0.136 | 0.295 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 18,319/15,082 | 0.939 | 0.347/0.601/0.440 | 0.281/0.719/0.405 | 0.075/0.418/0.128 | 0.273/0.550/0.365 |
| 10 | 17,511/14,842 | 0.939 | 0.369/0.598/0.457 | 0.301/0.722/0.425 | 0.084/0.427/0.141 | 0.288/0.547/0.378 |
| 15 | 16,671/14,470 | 0.939 | 0.380/0.595/0.464 | 0.319/0.725/0.443 | 0.089/0.432/0.148 | 0.297/0.549/0.386 |
| 20 | 15,952/14,088 | 0.939 | 0.380/0.596/0.464 | 0.329/0.728/0.453 | 0.092/0.420/0.151 | 0.296/0.545/0.384 |
| 30 | 14,890/13,522 | 0.939 | 0.382/0.604/0.468 | 0.354/0.740/0.479 | 0.093/0.417/0.152 | 0.300/0.549/0.388 |
| 40 | 14,084/13,054 | 0.939 | 0.374/0.624/0.468 | 0.370/0.748/0.495 | 0.091/0.438/0.151 | 0.300/0.562/0.392 |

### Janich liver (attn_union xspecies)  (xsp-attn_union universe)

> **SUPERSEDED (2026-07-30).** This Janich block was computed on the broken 2,434-tx universe (salmon run on untrimmed ARTseq reads, 0.03% mapping). Quarantined under `data/_quarantine/janich_stale_universe_2026-07-30/`. See the re-run further down this file, on the rebuilt 10,431-tx universe.
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/attn_union/mouse_janich_liver`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 2,050 (CDS 1,656, uORF 143, novel 51, dORF 101)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.991/0.998/0.995 | 0.669/0.804 | 0.540/0.667 | 0.140/0.436 | 0.453 |
| pred_preddepth (standalone theta=1) | 0.991/0.998/0.995 | 0.721/0.650 | 0.460/0.569 | 0.195/0.228 | 0.455 |
- standalone theta=1 over-call: total 2,052 calls vs 2,050 real (all 1.001x; novel 1.235x, dORF 1.168x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.994/0.923/0.957 | 1.000/0.161 | 0.522/0.235 | 0.375/0.030 | 0.223 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 2,367/2,313 | 0.994 | 0.680/0.656/0.668 | 0.414/0.545/0.471 | 0.155/0.226/0.184 | 0.468/0.498/0.483 |
| 10 | 2,253/2,235 | 0.995 | 0.717/0.658/0.686 | 0.430/0.540/0.479 | 0.169/0.225/0.193 | 0.482/0.487/0.484 |
| 15 | 2,148/2,139 | 0.995 | 0.746/0.676/0.710 | 0.431/0.544/0.481 | 0.190/0.234/0.210 | 0.485/0.482/0.484 |
| 20 | 2,052/2,050 | 0.995 | 0.721/0.650/0.684 | 0.460/0.569/0.509 | 0.195/0.228/0.210 | 0.461/0.449/0.455 |
| 30 | 1,919/1,928 | 0.995 | 0.714/0.633/0.671 | 0.578/0.634/0.605 | 0.225/0.225/0.225 | 0.466/0.430/0.447 |
| 40 | 1,833/1,833 | 0.995 | 0.657/0.697/0.676 | 0.686/0.706/0.696 | 0.196/0.175/0.185 | 0.455/0.424/0.439 |

### GSE120762 NT (attn_union xspecies)  (xsp-attn_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/attn_union/mouse_gse120762_nt`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 8,058 (CDS 7,654, uORF 223, novel 8, dORF 9)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.968/0.938/0.953 | 0.351/0.596 | 0.179/0.875 | 0.043/0.222 | 0.410 |
| pred_preddepth (standalone theta=1) | 0.954/0.955/0.954 | 0.184/0.892 | 0.123/0.875 | 0.009/0.333 | 0.267 |
- standalone theta=1 over-call: total 9,594 calls vs 8,058 real (all 1.191x; novel 7.125x, dORF 39.111x)
| Poisson CDS-anchored (th=0.02, CDSrec 0.88) | 0.947/0.879/0.912 | 0.545/0.300 | 0.350/0.875 | 0.067/0.111 | 0.340 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 11,332/8,158 | 0.954 | 0.112/0.908/0.200 | 0.080/0.875/0.146 | 0.006/0.300/0.011 | 0.109/0.795/0.192 |
| 10 | 10,653/8,144 | 0.954 | 0.140/0.912/0.242 | 0.089/0.875/0.161 | 0.007/0.300/0.013 | 0.131/0.798/0.225 |
| 15 | 10,063/8,106 | 0.954 | 0.165/0.901/0.279 | 0.103/0.875/0.184 | 0.007/0.300/0.015 | 0.148/0.785/0.249 |
| 20 | 9,594/8,058 | 0.954 | 0.184/0.892/0.305 | 0.123/0.875/0.215 | 0.009/0.333/0.017 | 0.161/0.770/0.267 |
| 30 | 8,948/7,970 | 0.954 | 0.215/0.882/0.346 | 0.189/0.875/0.311 | 0.011/0.375/0.021 | 0.181/0.734/0.290 |
| 40 | 8,529/7,886 | 0.954 | 0.247/0.878/0.386 | 0.259/1.000/0.412 | 0.014/0.500/0.028 | 0.194/0.721/0.306 |

### T-cell GSE155087 (attn_union xspecies)  (xsp-attn_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/attn_union/mouse_gse155087_tcell`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 11,417 (CDS 9,159, uORF 1,107, novel 275, dORF 211)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.975/0.981/0.978 | 0.680/0.581 | 0.501/0.684 | 0.165/0.337 | 0.497 |
| pred_preddepth (standalone theta=1) | 0.956/0.964/0.960 | 0.636/0.630 | 0.454/0.727 | 0.162/0.355 | 0.501 |
- standalone theta=1 over-call: total 11,658 calls vs 11,417 real (all 1.021x; novel 1.604x, dORF 2.19x)
| Poisson CDS-anchored (th=0.02, CDSrec 0.87) | 0.953/0.872/0.911 | 0.723/0.054 | 0.580/0.171 | 0.304/0.033 | 0.108 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 13,558/12,792 | 0.960 | 0.564/0.631/0.596 | 0.393/0.696/0.502 | 0.125/0.344/0.183 | 0.457/0.543/0.496 |
| 10 | 12,841/12,372 | 0.960 | 0.607/0.636/0.621 | 0.415/0.705/0.523 | 0.140/0.346/0.200 | 0.481/0.540/0.509 |
| 15 | 12,215/11,887 | 0.960 | 0.633/0.633/0.633 | 0.436/0.721/0.543 | 0.151/0.346/0.210 | 0.488/0.533/0.509 |
| 20 | 11,658/11,417 | 0.960 | 0.636/0.630/0.633 | 0.454/0.727/0.559 | 0.162/0.355/0.223 | 0.484/0.519/0.501 |
| 30 | 10,852/10,726 | 0.960 | 0.636/0.599/0.617 | 0.484/0.740/0.585 | 0.191/0.379/0.254 | 0.477/0.492/0.485 |
| 40 | 10,345/10,192 | 0.960 | 0.607/0.604/0.606 | 0.498/0.767/0.604 | 0.187/0.389/0.253 | 0.454/0.486/0.469 |

### GSE120762 LPS (attn_union xspecies)  (xsp-attn_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/attn_union/mouse_gse120762_lps`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 9,196 (CDS 8,144, uORF 574, novel 57, dORF 61)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.970/0.966/0.968 | 0.544/0.681 | 0.331/0.754 | 0.109/0.279 | 0.503 |
| pred_preddepth (standalone theta=1) | 0.953/0.956/0.955 | 0.442/0.840 | 0.261/0.842 | 0.045/0.311 | 0.455 |
- standalone theta=1 over-call: total 10,317 calls vs 9,196 real (all 1.122x; novel 3.228x, dORF 6.885x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.92) | 0.950/0.916/0.933 | 0.656/0.300 | 0.370/0.474 | 0.038/0.033 | 0.344 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 12,129/9,749 | 0.955 | 0.349/0.851/0.495 | 0.238/0.868/0.374 | 0.036/0.314/0.065 | 0.296/0.731/0.422 |
| 10 | 11,450/9,606 | 0.955 | 0.394/0.860/0.541 | 0.255/0.861/0.394 | 0.041/0.319/0.073 | 0.325/0.730/0.450 |
| 15 | 10,822/9,390 | 0.955 | 0.421/0.847/0.563 | 0.243/0.850/0.378 | 0.042/0.308/0.074 | 0.334/0.710/0.454 |
| 20 | 10,317/9,196 | 0.955 | 0.442/0.840/0.579 | 0.261/0.842/0.398 | 0.045/0.311/0.079 | 0.339/0.692/0.455 |
| 30 | 9,608/8,886 | 0.955 | 0.464/0.830/0.595 | 0.324/0.880/0.473 | 0.051/0.308/0.087 | 0.343/0.664/0.452 |
| 40 | 9,143/8,643 | 0.955 | 0.436/0.816/0.569 | 0.363/0.902/0.517 | 0.059/0.342/0.101 | 0.317/0.618/0.419 |

### Wang liver (attn_union xspecies)  (xsp-attn_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/attn_union/mouse_wang_liver`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 14,088 (CDS 11,963, uORF 774, novel 673, dORF 162)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.963/0.958/0.960 | 0.445/0.431 | 0.381/0.670 | 0.128/0.340 | 0.391 |
| pred_preddepth (standalone theta=1) | 0.939/0.944/0.941 | 0.348/0.689 | 0.311/0.731 | 0.073/0.531 | 0.364 |
- standalone theta=1 over-call: total 16,843 calls vs 14,088 real (all 1.196x; novel 2.349x, dORF 7.315x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.90) | 0.938/0.901/0.919 | 0.572/0.133 | 0.413/0.452 | 0.255/0.222 | 0.298 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 19,686/15,082 | 0.941 | 0.319/0.671/0.432 | 0.260/0.720/0.382 | 0.061/0.531/0.109 | 0.242/0.593/0.343 |
| 10 | 18,690/14,842 | 0.941 | 0.341/0.678/0.454 | 0.282/0.718/0.405 | 0.068/0.542/0.121 | 0.257/0.595/0.359 |
| 15 | 17,734/14,470 | 0.941 | 0.348/0.681/0.461 | 0.298/0.725/0.422 | 0.070/0.534/0.123 | 0.263/0.598/0.365 |
| 20 | 16,843/14,088 | 0.941 | 0.348/0.689/0.462 | 0.311/0.731/0.437 | 0.073/0.531/0.128 | 0.262/0.593/0.364 |
| 30 | 15,485/13,522 | 0.941 | 0.349/0.693/0.464 | 0.337/0.739/0.462 | 0.077/0.543/0.135 | 0.267/0.590/0.367 |
| 40 | 14,523/13,054 | 0.941 | 0.330/0.714/0.451 | 0.357/0.759/0.486 | 0.075/0.573/0.133 | 0.265/0.604/0.369 |

### Janich liver (attn_union xspecies)  (xsp-attn_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/attn_union/mouse_janich_liver`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 7,737 (CDS 6,154, uORF 616, novel 338, dORF 226)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.984/0.987/0.985 | 0.682/0.737 | 0.601/0.645 | 0.159/0.425 | 0.518 |
| pred_preddepth (standalone theta=1) | 0.978/0.983/0.981 | 0.693/0.718 | 0.553/0.618 | 0.149/0.301 | 0.515 |
- standalone theta=1 over-call: total 8,047 calls vs 7,737 real (all 1.04x; novel 1.118x, dORF 2.013x)
| Poisson CDS-anchored (th=0.05, CDSrec 0.94) | 0.979/0.938/0.958 | 0.862/0.182 | 0.608/0.225 | 0.362/0.075 | 0.264 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 9,444/8,675 | 0.980 | 0.616/0.721/0.665 | 0.485/0.573/0.525 | 0.108/0.281/0.156 | 0.455/0.589/0.513 |
| 10 | 8,943/8,410 | 0.981 | 0.662/0.724/0.691 | 0.511/0.576/0.542 | 0.121/0.284/0.169 | 0.476/0.582/0.524 |
| 15 | 8,473/8,082 | 0.981 | 0.693/0.725/0.709 | 0.539/0.602/0.569 | 0.136/0.295/0.186 | 0.485/0.577/0.527 |
| 20 | 8,047/7,737 | 0.981 | 0.693/0.718/0.705 | 0.553/0.618/0.584 | 0.149/0.301/0.200 | 0.476/0.562/0.515 |
| 30 | 7,457/7,252 | 0.981 | 0.701/0.717/0.709 | 0.602/0.652/0.626 | 0.182/0.308/0.229 | 0.471/0.547/0.506 |
| 40 | 7,045/6,877 | 0.981 | 0.659/0.697/0.678 | 0.621/0.681/0.649 | 0.178/0.298/0.223 | 0.452/0.540/0.492 |

### Janich liver (mamba4_union xspecies)  (xsp-mamba4_union universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_janich_liver`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 7,737 (CDS 6,154, uORF 616, novel 338, dORF 226)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.983/0.988/0.986 | 0.712/0.745 | 0.646/0.654 | 0.233/0.274 | 0.568 |
| pred_preddepth (standalone theta=1) | 0.975/0.959/0.967 | 0.754/0.492 | 0.631/0.592 | 0.209/0.102 | 0.486 |
- standalone theta=1 over-call: total 7,079 calls vs 7,737 real (all 0.915x; novel 0.938x, dORF 0.487x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.92) | 0.976/0.919/0.946 | 0.857/0.078 | 0.692/0.325 | 0.240/0.026 | 0.205 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 7,732/8,675 | 0.967 | 0.707/0.452/0.551 | 0.554/0.551/0.553 | 0.184/0.103/0.132 | 0.592/0.394/0.473 |
| 10 | 7,513/8,410 | 0.967 | 0.746/0.466/0.573 | 0.589/0.545/0.567 | 0.188/0.101/0.132 | 0.614/0.397/0.482 |
| 15 | 7,296/8,082 | 0.967 | 0.757/0.481/0.588 | 0.613/0.566/0.589 | 0.200/0.102/0.136 | 0.625/0.402/0.489 |
| 20 | 7,079/7,737 | 0.967 | 0.754/0.492/0.595 | 0.631/0.592/0.611 | 0.209/0.102/0.137 | 0.619/0.400/0.486 |
| 30 | 6,780/7,252 | 0.967 | 0.743/0.500/0.598 | 0.680/0.630/0.654 | 0.256/0.112/0.156 | 0.618/0.408/0.491 |
| 40 | 6,563/6,877 | 0.967 | 0.707/0.537/0.610 | 0.678/0.672/0.675 | 0.250/0.107/0.150 | 0.609/0.428/0.503 |

### Janich liver (mamba4_union xspecies _mm10)  (xsp-mamba4_union_mm10 universe)
- run: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/heldout_xspecies/mamba4_union/mouse_janich_liver_mm10`
- **primary floor = ORFs >= 20 aa** (headline tables); sweep below spans full range.
- reference = observed Hepatocytes calls (>= 20 aa): 7,737 (CDS 6,154, uORF 616, novel 338, dORF 226)

| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |
|---|---|---|---|---|---|
| pred_obsdepth (shape @ real depth) | 0.984/0.990/0.987 | 0.709/0.722 | 0.629/0.651 | 0.285/0.270 | 0.573 |
| pred_preddepth (standalone theta=1) | 0.973/0.964/0.969 | 0.767/0.390 | 0.603/0.571 | 0.357/0.088 | 0.447 |
- standalone theta=1 over-call: total 6,954 calls vs 7,737 real (all 0.899x; novel 0.947x, dORF 0.248x)
| Poisson CDS-anchored (th=0.1, CDSrec 0.90) | 0.973/0.898/0.934 | 0.750/0.044 | 0.664/0.299 | 0.900/0.040 | 0.171 |

_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:
| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 | non-canon P/R/F1 |
|---|---|---|---|---|---|---|
| 5 | 7,438/8,675 | 0.969 | 0.726/0.348/0.470 | 0.545/0.519/0.532 | 0.311/0.087/0.137 | 0.614/0.327/0.426 |
| 10 | 7,285/8,410 | 0.969 | 0.762/0.362/0.491 | 0.569/0.520/0.543 | 0.319/0.089/0.140 | 0.634/0.334/0.437 |
| 15 | 7,132/8,082 | 0.969 | 0.771/0.383/0.511 | 0.584/0.540/0.561 | 0.328/0.090/0.141 | 0.640/0.343/0.447 |
| 20 | 6,954/7,737 | 0.969 | 0.767/0.390/0.517 | 0.603/0.571/0.587 | 0.357/0.088/0.142 | 0.636/0.344/0.447 |
| 30 | 6,714/7,252 | 0.969 | 0.750/0.402/0.523 | 0.645/0.605/0.624 | 0.439/0.096/0.157 | 0.636/0.358/0.458 |
| 40 | 6,536/6,877 | 0.969 | 0.726/0.423/0.534 | 0.644/0.647/0.645 | 0.481/0.099/0.165 | 0.627/0.383/0.476 |
