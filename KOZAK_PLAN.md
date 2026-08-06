# Plan: Kozak start-context ablation (remove / empirical / learned)

Finalized 2026-07-16. Questions: (1) how do metrics hold if the heuristic Kozak is removed, (2) does an
empirically-fit Kozak beat it, (3) can the model learn the context weights itself. Scope decisions (user):
**4 arms (V0/V1/V2/V3), one-hot backend only, train Fibroblast -> test Hepatocytes.**

## Background

The ORF-track `ext` (v2) start channel (channel 3, `build_orf_track.py`) is:
`start_ext[i] = START_W[codon] x (0.5 + 0.5*kz)`, `kz = 0.5*[purine at -3] + 0.5*[G at +4]`.
The `{0.5,0.5}` weights at exactly `-3`/`+4` are a hand-picked heuristic. Occupancy channels 0-2 and stop
channel 4 are ATG-based and identical across all arms; ONLY channel 3 changes, so every arm isolates the
start-context model.

## Shared setup (isolates channel 3)

- Train: Fibroblast (`data/packed/`). Test: Hepatocytes (`data/packed_Hepatocytes/`). Offsets identical
  (md5-verified), so a new orf_track drops in via `RIBO_ORF_TRACK`.
- `train_loto.py --train_tissues Fibroblast --holdout Hepatocytes --emb_backend onehot --use_orf_track
  --n_attn_layers 2 --cov_norm global_mean` + standard hyperparams. Single-tissue train, cross-tissue test.
- Anchor V0 is trained fresh under THIS protocol (existing orf_v2_attn models are train-on-8, not comparable).

## Arms (channel 3 only differs)

| arm | channel 3 (start propensity) | script |
|-----|------------------------------|--------|
| V0 heuristic | `START_W[codon] x (0.5+0.5*kz_heuristic)` | existing `orf_track_v2.npy` (reproduce as check) |
| V1 no-Kozak  | `START_W[codon]` (drop the Kozak factor) | `build_orf_track.py --mode ext --kozak none` |
| V2 empirical | `START_W[codon] x m_pwm` (PWM from annotated starts) | `--kozak pwm --pwm_file kozak_pwm.json` |
| V3 learnable | `START_W[codon]` in the track + a learnable conv over one-hot produces the gate | `model.py --learn_start_context` on the V1 track |

## Q2 empirical PWM (`build_kozak_pwm.py`)

- Source: annotated CDS starts (`tx2cds.tsv` has_start_codon=1, start at `utr5_len`) with flanks from
  `fibroblast_universe.fa`. Context positions = {-6,-5,-4,-3,-2,-1,+4} (exclude the ATG triplet -- constant,
  captured by START_W). PFM -> log-odds vs the universe background nt frequency (pseudocount 1).
- Candidate multiplier: `m_pwm = clip( 0.5 + 0.5*(score - p5)/(p95 - p5), 0.5, 1.0 )`, score = summed
  log-odds over the 7 positions; p5/p95 from the annotated-start score distribution (stored in the JSON).
  Matches the heuristic's [0.5,1] range so V2-vs-V0 is a clean swap.
- NO test-label leakage: PWM is fit from annotation SEQUENCE, not Hepatocytes Ribo-seq (same posture as
  the literature-based heuristic).

## Q3 learnable module (`model.py --learn_start_context`, one-hot only)

- Input = the V1 (no-Kozak) track (channel 3 = codon weight only). Inside the model, a
  `Conv1d(4 -> 1, kernel 11)` reads the one-hot slice `feats[:, :, 0:4]` with asymmetric pad (6 left, 4
  right) so output[i] depends on one-hot[i-6 .. i+4]; `gate = sigmoid(conv)`; the start channel is
  recomputed `w x (0.5 + 0.5*gate)` before `in_proj`. The 11x4 kernel = a LEARNED position x nucleotide
  Kozak matrix, extracted post-training and plotted vs the empirical PWM and the heuristic {-3 purine,+4 G}.
  - AS-BUILT correction (this pre-registration text was left intact as a record): implemented with
    `kernel_size=10`, pad (6,3), so output[i] reads one-hot[i-6..i+3] = Kozak {-6,-5,-4,-3,-2,-1,+1,+2,+3,+4},
    a 4x10 kernel (NOT k=11 / pad (6,4) / i-6..i+4). `model.py` and the Outcome section below are
    authoritative; methods.md 5N, results.md Task 20, and figures/kozak/FIGURE_DATA_INPUTS.md all use k=10.
- Behavioral Q3 (free): V1-vs-V0 -- the one-hot backbone already has -3/+4 in its receptive field, so if it
  can learn Kozak, V1 recovers V0.

## Eval per arm (foreground what Kozak touches)

Kozak modulates the START channel -> it should move non-AUG / uORF detection, not CDS. So rank metrics:
1. **CTG non-AUG drop-in** (precision/recall/F1) -- reuse `ribocode_dropin_ctg.sbatch` + `compare_dropin_ctg.py`. Sharpest test.
2. **uORF**: 5'UTR profile Pearson (`eval_extra.py`) + uORF localization AUROC (`eval_localization.py`).
3. ATG drop-in F1 (expect ~flat) + whole-pc profile Pearson (sanity, ~flat).
4. V3: learned-weights figure (heatmap/logo of the 11x4 kernel vs empirical PWM vs heuristic).

Each arm needs a `pred_profiles.npz` dump (`dump_pred_profiles.py`, GPU ~15 min) before the drop-ins.

## Execution order

1. `build_kozak_pwm.py` -> `data/kozak_pwm.json` (+ print PWM; first look at real Kozak).           [head node]
2. Extend `build_orf_track.py` with `--kozak {heuristic,none,pwm}` + `--pwm_file`; build
   `orf_track_v2_nokozak.npy`, `orf_track_v2_pwm.npy`; reproduce `orf_track_v2.npy` as a check.        [short]
3. `model.py` + `train_loto.py`: add `--learn_start_context`; CPU smoke both paths.                    [head node]
4. Train 4 arms (V0/V1/V2 on their tracks; V3 = V1 track + `--learn_start_context`).                   [gpu, ~2-3h each parallel]
5. Per arm: dump pred_profiles -> localization + ATG drop-in + CTG drop-in.                            [gpu + medium]
6. Aggregate; results.md Task 20 + methods.md + figure (learned/empirical/heuristic Kozak + metric bars).

## Files

New: `KOZAK_PLAN.md`, `scripts/build_kozak_pwm.py`, `scripts/train_kozak_arms.sbatch`,
`scripts/eval_kozak_arms.sbatch`, `scripts/plot_kozak_weights.py`, `data/kozak_pwm.json`,
`data/packed/orf_track_v2_{nokozak,pwm}.npy`. Edited: `build_orf_track.py`, `model.py`, `train_loto.py`,
`dataset.py` (learn_start_context plumbing). Results under `results/kozak/{v0,v1,v2,v3}_onehot_fib2hep/`.

## Pre-training finding (2026-07-16): is the Kozak context different for CTG/alt ORFs?

`kozak_context_alt_orfs.py` on RiboCode's real CTG-aware Hepatocytes calls (dropin_ctg/real_collapsed.txt),
context scored under the annotated-CDS PWM:

| group | n | -3 purine | +4 G | median annPWM score | % below p5 (floor) |
|-------|---|-----------|------|---------------------|--------------------|
| annotated CDS (ATG) | 32332 | 0.858 | 0.515 | ~1.55 | (5% by def) |
| RiboCode ATG ORFs (all) | 12935 | 0.801 | 0.503 | 1.35 | 7.6% |
| RiboCode CTG ORFs (all) | 1443 | 0.599 | 0.423 | 0.14 | 15.0% |
| **ATG-uORF only** | 1346 | 0.498 | 0.326 | -0.69 | 25.0% |
| **CTG-uORF only** | 1329 | 0.600 | 0.426 | 0.15 | 14.5% |

At face value CTG ORFs have a very different, weaker, G-rich context (G at -4/-3/-2/-1 vs the canonical
C-rich A(-3) Kozak). BUT that is mostly a LOCATION confound: CTG ORFs are ~92% uORFs, and controlling for
it (ATG-uORF vs CTG-uORF) the codon difference vanishes -- both uORF groups are G-rich and weak, and
CTG-uORF is if anything marginally MORE Kozak-like than ATG-uORF (-3 purine 0.60 vs 0.50). The real
dichotomy is **uORF vs CDS, not CTG vs ATG**: uORF starts (either codon) sit in a weak G-rich 5'UTR
context, the strong C-rich Kozak is a canonical-CDS property.

**Prediction this sets up:** the annotated-CDS PWM (V2) is a MISCALIBRATED prior for the non-canonical
starts the project targets -- it pushes ~15-25% of uORFs (both codons) to the 0.5 floor. So V2 may NOT beat
V1 (no Kozak) on uORF / CTG detection, because canonical Kozak is weak/absent at these sites. The 4 running
arms test this directly; if so, the lesson is that a CDS-fit context prior is the wrong tool for alt-ORF
detection (a uORF-fit PWM or a learned/location-aware gate would be needed). Caveat: RiboCode uORF calls
are a soft truth, but the ATG-uORF vs CTG-uORF contrast is internally consistent (same caller).

## Outcome (2026-07-16) -- LANDED (results.md Task 20)

All 4 arms trained + evaluated (Fib->Hep). Full table in `results/kozak/kozak_summary.md`.

**In-distribution the arms are equivalent:** best Fibroblast-val Pearson 0.6418-0.6458 (spread 0.004),
V1 (no-Kozak) actually the lowest. So the cross-tissue test deltas are generalization, not fit quality;
single seed each, soft CTG truth -- directional, not decisive.

- **Q1 (remove):** V1 no-Kozak is directionally BEST across ~7 held-out metrics (CTG F1 0.369 vs V0 0.340;
  uORF frame0 Pearson 0.503 vs 0.461; non-canon AUROC 0.848; pc 0.587; 5'UTR 0.536; ATG F1 0.923). The
  heuristic was a mild net negative.
- **Q2 (empirical PWM):** V2 CTG F1 0.343 == V0 0.340 -- a wash, below V1. Pre-registered prediction
  CONFIRMED: fitting a canonical-CDS Kozak better does not help the ~92%-uORF CTG starts.
- **Q3 (learned):** V3's learned 4x10 kernel rediscovers Kozak unsupervised -- -3 purine (A+G +0.759 vs
  C+T -0.808), correctly weak +4, and **r = 0.807 vs the empirical PWM**. But the explicit sigmoid gate is
  the WORST detector (most conservative: 720 CTG calls, recall 0.238, CTG F1 0.320). The backbone already
  learns start context implicitly, so an explicit gate is redundant/harmful.

**Recommendation: drop the heuristic (`build_orf_track.py --kozak none`, the V1 track).** Redundant with
what the sequence backbone learns on its own (Q3 r=0.807), simplifies the track, directionally improves
held-out non-canonical detection. Neither PWM nor learned gate earns its complexity. Follow-up before a
deployment rebuild: 3-seed confirmation of the V1 edge (deltas are modest, single-seed).

## Standing rules

ASCII only (no en/em dashes); no "our"; no AI-authorship; no commit/push without ask; scratch to
/data/tmp/emalekos; group ceph is 15 TB SHARED (watch EDQUOT); maintain methods.md/results.md/logs.
