# Model architecture

A dual-head dilated-CNN with a global-context mixer (`scripts/model.py`, `RiboSignalModel`). About 5M
parameters. Two mixers are released: **`mamba4`** (primary, GPU only) and **`attn`** (two transformer layers, CPU-compatible). Both use one-hot input and the no-Kozak ORF track.

```{figure} /img/arch_attn/arch_rinalmo.png
:width: 100%

The `attn` model end to end (the CPU-compatible release; `mamba4` differs only in the mixer, below). **Left:** the per-nucleotide input encoding, drawn as the literal matrices
the model receives (one-hot sequence, the 5-channel ORF-candidate track, RNA-seq coverage), which is where
the `(B, L, 10)` tensor comes from. **Centre:** the encoder, a 1x1 conv projection into 10 residual dilated
blocks (dilation 1 to 512, receptive field ~4 kb) followed by 2 pre-norm transformer blocks. **Right:** the
two heads -- profile (shape) and count (depth). Regenerate with
`cd figures/arch_attn && python make_arch_panel.py`.
```

## Inputs (per transcript, per nucleotide)

```{list-table}
:header-rows: 1
:widths: 24 14 62

* - Input
  - Dim
  - What it carries
* - Sequence
  - 4
  - one-hot A/C/G/T (default). Optionally swapped for a frozen RNA-FM embedding (RiNALMo 1280-d, Orthrus
    512-d, HydraRNA 1024-d) -- see {doc}`/benchmarks`.
* - RNA-seq coverage
  - 1
  - per-nt depth from unique-mapper STAR, globally mean-normalised. This is what makes the prediction
    tissue-specific.
* - ORF track
  - 5
  - a **reading-frame prior computed from sequence alone**: `orf_f0/f1/f2` (nt inside an ATG..in-frame-stop
    ORF, by frame), `start_ext` (graded start-codon propensity), `is_stop`. Defined identically for CDS,
    uORFs, dORFs, and lncRNA ORFs, so it does not leak the annotated CDS. See
    [the ORF-candidate track](orf-track) below.
```

(orf-track)=
## The ORF-candidate track, in detail

Built by `scripts/build_orf_track.py`. This is the input people most often misread, so it is worth being
precise: **the track encodes every candidate ORF at once, not one chosen ORF.**

### Two different kinds of information

The track tells the model two separate things about starts:

1. **"a candidate start codon sits here"** -- the `start_ext` channel, a single graded value at one nucleotide.
2. **"this nucleotide is inside an ORF body"** -- the `orf_f0/f1/f2` occupancy channels, a filled span running
   from a start to its in-frame stop.

Occupancy is bucketed by the **frame of the start codon** (`start_pos % 3`), so two ORFs that overlap in
different frames light up two channels simultaneously. `is_stop` marks the first nt of **every** stop codon in
**any** frame, not just the one closing a particular ORF. Within a single occupancy channel the scan steps by
three, so a span's start and stop are always in-frame with each other.

### Worked example

Running the production builder on a 24-nt window (`orf_track(seq, mode="ext", kozak="none")`):

```text
seq        A U G C A U G C U G A A G C C U A A C U A G C C
orf_f0     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
orf_f1     0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0
orf_f2     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
start_ext  1 . . . 1 . . .5 . . .15 . . . . . . . . . . . .
is_stop    . . . . . . . . 1 . . . . . . 1 . . . 1 . . . .
```

AUG at 0 opens an ORF in frame 0 (closing at the in-frame stop at 15); AUG at 4 opens a second ORF in frame 1
(closing at 19). **Both are occupied across positions 4 to 17.** The three stops land in three different
frames (8, 15, 19). The CUG at 7 and AAG at 10 get graded start weights but no span, for the reason below.

### Start-codon weights

`start_ext` is the product of a **per-codon weight** and a **context multiplier**:

```text
start_ext[i] = START_W[codon at i] x mult[i]
```

`START_W` covers AUG plus the nine near-cognates one substitution away:

```{list-table}
:header-rows: 1
:widths: 30 14 56

* - Codon
  - Weight
  - Note
* - AUG
  - 1.00
  - canonical
* - CUG
  - 0.50
  - strongest non-AUG start
* - GUG, ACG
  - 0.35
  -
* - UUG
  - 0.30
  -
* - AUA
  - 0.25
  -
* - AUU, AUC
  - 0.20
  -
* - AAG, AGG
  - 0.15
  -
```

```{admonition} These weights are a hand-picked prior, not fitted
:class: warning
The values are **ordinal, taken from the non-AUG initiation literature** (AUG far above CUG, then GUG/ACG,
then the rest). The rank ordering is well supported; the exact magnitudes are judgment calls. Nothing in the
pipeline tunes them against data.
```

### Why occupancy opens at AUG only

The occupancy channels use `start_mask = atg`. Non-AUG initiation therefore reaches the model **only** as a
start pulse, never as a filled ORF body:

```text
AUG-initiated ORF, 100 to 250        CUG-initiated ORF, 100 to 250
  start_ext  1.0 at 100                start_ext  0.5 at 100
  orf_fN     span 100 to 252           orf_fN     all zeros   <- no span
  is_stop    1 at 250                  is_stop    1 at 250
```

This is a deliberate density trade. AUG is 1 of 64 codons; the ten start codons together are 10 of 64. If
every one opened a span running to its next in-frame stop, the spans would overlap so heavily that occupancy
would sit near 1 almost everywhere in every frame, and the channel would stop carrying information. Keeping it
AUG-only preserves a sharp **~37% per frame** prior. A middle-ground mode (`--mode atgctg`) extends occupancy
to CUG alone; opening it to the full ten diluted an ATG lncRNA lift from +0.056 to +0.033.

The practical consequence: the model is **not blind to non-AUG ORFs**, but it gets less scaffolding for them,
so its non-canonical calls lean harder on measured periodicity than on the input prior.

### The context multiplier

```{list-table}
:header-rows: 1
:widths: 22 78

* - `--kozak`
  - Multiplier
* - `none` *(default)*
  - `mult = 1.0`. Codon weight only, no context.
* - `heuristic`
  - `mult = 0.5 + 0.5 x (0.5[purine at -3] + 0.5[G at +4])`, so `mult` is 0.50, 0.75, or 1.00.
* - `pwm`
  - Empirical PWM fitted by `build_kozak_pwm.py` from annotated CDS starts, scoring positions
    {-6..-1, +4} (the AUG triplet is excluded, since `START_W` already covers it) as log2-odds vs background,
    min-max calibrated onto [0.5, 1.0]. Uses annotation **sequence** only, never Ribo-seq signal.
```

## Body and heads

```{mermaid}
flowchart LR
    I["inputs<br/><small>seq + coverage + ORF track</small>"]
    B["dilated CNN body<br/><small>10 blocks, 256 ch</small>"]
    X["global mixer<br/><small>2-4 transformer layers, or Mamba</small>"]
    H1["profile head<br/><small>softmax over positions</small>"]
    H2["count head<br/><small>log total P-sites</small>"]
    I --> B --> X
    X --> H1
    X --> H2
    classDef s fill:#eaf2fb,stroke:#3a6ea5,color:#16314a;
    class I,B,X,H1,H2 s;
```

- **Profile head** -- a softmax over the transcript's positions: *where* along the transcript P-sites sit
  (scale-invariant shape).
- **Count head** -- the *total* predicted P-sites (depth). Shape x total = the predicted profile fed to the
  ORF caller.
- **Loss** = multinomial profile NLL + `count_weight` x count MSE.

## The mixer knob

The CNN body is fixed; the global-context **mixer** is the architecture variable we sweep. **Two settings
are released**, and the choice is a capability trade, not just an accuracy one:

```{list-table}
:header-rows: 1
:widths: 20 18 62

* - Mixer
  - Flag
  - Notes
* - **Bi-Mamba x4** (`mamba4`)
  - `--mixer mamba --n_attn_layers 4`
  - **primary model.** Best held-out test Pearson (0.680 vs 0.659 across 3 seeds, non-overlapping
    ranges) and wins 5/5 mouse datasets on cross-species CDS F1. 7,521,026 params. **GPU only** --
    `causal_conv1d` is CUDA-only.
* - **transformer x2** (`attn`)
  - *(default)*
  - **released for CPU-bound users** and the supplemental model in the paper. 5,071,106 params, runs
    anywhere. Measurably weaker but the only option without a GPU.
* - transformer x4 (`attn4`)
  - `--n_attn_layers 4`
  - looked best on a single seed; the gap vanishes across seeds. 6 layers *regresses*. Not released.
```

```{admonition} Which model should I use?
:class: tip
**`mamba4` if you have a GPU, `attn` otherwise.** The accuracy gap is real but modest (+0.021 test
Pearson, +0.021 mean cross-species non-canonical F1); the CPU/GPU distinction is absolute. Both are
maintained, both take the same inputs, and both are drop-in for the same downstream ORF-calling path.
```

Only the mixer region differs between the two candidate models, so the panels are directly comparable:

```{figure} /img/arch_attn/arch_rinalmo_mamba.png
:width: 100%

The `mamba4` variant (primary model). Identical to `attn` except the second dashed region: 4 **Bi-Mamba**
blocks (`LayerNorm` -> forward scan + reverse scan -> sum -> `Dropout` -> residual; `d_state` 16,
`d_conv` 4, `expand` 2) replace the 2 transformer blocks. Full-transcript context at O(L) instead of
O(L^2), at 7,521,026 parameters vs 5,071,106. Note it **cannot run on CPU** -- `causal_conv1d` is
CUDA-only.
```

```{admonition} No Kozak, ever
:class: important
The deployed model uses the **no-Kozak** ORF track (`--kozak none`). An earlier Kozak start-context heuristic
up-weighted near-cognate starts and fed spurious uORF/novel calls; ablation showed it is redundant or harmful,
and a learned kernel rediscovers Kozak unsupervised anyway. All current retrains set `--kozak none`.
```
