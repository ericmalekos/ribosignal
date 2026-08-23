# Panel A + B caption, rewritten in plain English (2026-08-17)

Same facts as the current caption, in language a passer-by can follow, plus the ORF / start /
stop track description and the non-canonical start weights that were missing.

Every number below was checked against an artifact. See "Provenance" at the bottom.

---

## Drop-in caption

**(A)** The model reads sequence through ten stacked convolution blocks. Each block looks further
out than the last, so by the end every position has seen 4,093 nt around it. Two attention layers
follow and let every position see every other one, so the model reads the transcript end to end.
5,071,106 weights in total.

The output splits in two. One head says *where* along the transcript the ribosomes sit, as
fractions that add to 1; the other says *how many* there are. Because the first head only ever
outputs fractions, sequencing deeper does not change it, and all of the depth sits in the second
head.

No Ribo-seq is used at inference. The inputs are the sequence itself and RNA-seq coverage, nothing
else. Ten channels enter the model: four for the base at each position, five for the ORF track, one
for coverage. The coverage profile drawn in the input block is schematic.

The ORF track marks every possible ORF at once, not a chosen one. Three of its channels are
occupancy: a position switches on if it lies anywhere from an ATG through the next in-frame stop,
and which of the three lights up records that ATG's reading frame. Overlapping ORFs are common, so
a position is often on in two frames at once. A fourth channel marks the first base of every
candidate start codon and grades it by how readily that codon is known to initiate: ATG 1.0,
CTG 0.5, GTG and ACG 0.35, TTG 0.3, ATA 0.25, ATT and ATC 0.2, AAG and AGG 0.15. These rankings
come from the non-AUG initiation literature and are set by hand, not learned. A fifth channel marks
the first base of every stop codon (TAA, TAG, TGA). Occupancy stays ATG-only on purpose: switching
it on at all ten start codons would cover most of the transcript and the channel would stop
carrying information. The weaker starts enter through the graded channel alone.

**(B)** Held-out hepatocytes, one exemplar per ORF class, picked rather than randomly sampled. The
annotated CDS row is cropped: its called ORF runs to 3,323 nt, but the window ends at 1,819, where
the observed signal stops. "r" is the Pearson correlation across the whole transcript, including
the part outside the window.

---

## Shorter ORF-track paragraph, if the full one does not fit

Keeps every fact, roughly half the length. Swap it for the fourth paragraph of (A).

> The ORF track marks every possible ORF at once. Three channels switch on from each ATG through
> its next in-frame stop, one channel per reading frame, so overlapping ORFs stay separable. A
> fourth marks every candidate start codon, graded by how readily it initiates (ATG 1.0, CTG 0.5,
> GTG and ACG 0.35, TTG 0.3, ATA 0.25, ATT and ATC 0.2, AAG and AGG 0.15, set from the literature,
> not learned). A fifth marks every stop codon. Only ATG opens an occupancy channel; the weaker
> starts enter through the graded channel alone, which is what stops occupancy from covering the
> whole transcript.

---

## One sentence to have in your pocket

A reviewer who knows the field will ask whether the start weights are modulated by Kozak context.
They are not, and that is deliberate: the deployed model is built with `--kozak none`, so the
graded channel carries the raw codon weight with no context multiplier. The ablation found the
heuristic gate redundant at best and harmful cross-tissue.

Add to the caption only if there is room:

> Start weights are raw codon weights; Kozak context is deliberately not applied.

---

## What the jargon became

| original wording | plain wording |
|---|---|
| residual dilated blocks | stacked convolution blocks, each looking further out than the last |
| 4,093 nt receptive field | every position has seen 4,093 nt around it |
| two transformer blocks read the whole transcript | two attention layers let every position see every other one |
| profile head | the head that says *where* the signal falls |
| scale-free | outputs fractions that add to 1, so sequencing deeper does not change it |
| count head | the head that says *how much* signal there is |
| ORF-candidate track encodes every candidate ATG ORF | the ORF track marks every possible ORF at once |
| whole-transcript Pearson | Pearson across the whole transcript, including the part outside the window |

---

## Two accuracy notes on the rewrite

**Occupancy includes the stop codon.** The three frame channels run from the ATG's first base
*through the last base of the in-frame stop*, not up to it. Verified on the panel's own 24-nt
display window `ATGCATGCTGAAGCCTAACTAGCC`: channel `orf_f0` is on across positions 0 to 17, and the
in-frame stop occupies 15 to 17. So "from an ATG through the next in-frame stop" is the correct
phrasing and "between an ATG and the next stop" is not.

**Two frames on at once is the normal case, not an edge case.** In that same window the ATG at
position 0 and the ATG at position 4 sit in different frames (0 % 3 = 0 and 4 % 3 = 1), so `orf_f0`
and `orf_f1` are both switched on across the overlap. `orf_f2` is empty because no ATG falls in
that frame. This is why the track needs three channels rather than one, and it is worth saying out
loud in the caption because a reader who assumes one ORF per transcript will misread the panel.

Occupancy channel values are binary 0 or 1, not counts.

---

## Provenance

Nothing here was transcribed from a previous caption or from memory.

| claim | traced to |
|---|---|
| 5,071,106 weights | summed from `best.pt` state_dict by `scripts/export_arch_spec.py` |
| 4,093 nt | `1 + 2*(k-1)*sum(dilations)` = `1 + 4*1023`, computed from checkpoint tensor shapes |
| ten blocks, two attention layers | `n_blocks` from the checkpoint, `n_attn_layers` from `args.json` |
| ten input channels, coverage last | `in_proj.weight` shape; `input_channels.order` in `arch_spec_attn.json` |
| the ten start weights | `START_W` imported from `scripts/build_orf_track.py` |
| five ORF-track channels and their meanings | `orf_track_spec.json`, `channels` |
| occupancy is binary, spans ATG through stop, two frames overlap | computed by calling `build_orf_track.orf_track()` on the display window |
| Kozak not applied | `kozak: "none"` in `orf_track_spec.json`; policy in `docs/PIPELINE_POLICY.md` |
| ORF runs to 3,323 nt, window ends at 1,819 | per-nt exemplar TSV; 1 P-site of 13,258 lies at or beyond 1,819 |
| "r" is `r_full` | `scripts/plot_profile_exemplars.py` |

Portable copies of both spec files are already in the poster package at
`figures/arch_attn/arch_spec_attn.json` and `figures/arch_attn/orf_track_spec.json`.
