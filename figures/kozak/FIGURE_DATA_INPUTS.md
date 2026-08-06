# FIGURE_DATA_INPUTS -- figures/kozak/

Traces every panel of each figure in this folder back to its data source. Update on every figure edit.

## learned_vs_empirical_vs_heuristic.png (+ .json)

Q3 payoff figure: the start-context weights the V3 model LEARNED vs the empirical Kozak PWM (V2) vs the
hand-picked heuristic (V0). Three stacked position x nucleotide heatmaps (Kozak positions
[-6,-5,-4,-3,-2,-1,+1,+2,+3,+4] on x, nucleotides A,C,G,T on y). Generator: `scripts/plot_kozak_weights.py`.

| panel | content | data source |
|-------|---------|-------------|
| top (V3 LEARNED) | the `start_ctx_conv` Conv1d(4->1, k=10) weight, reshaped (4 nt x 10 pos), raw conv scale | `results/kozak/v3_learned_onehot_fib2hep/best.pt` (state-dict key `*start_ctx_conv.weight`; bias from `*.bias`) |
| middle (V2 EMPIRICAL) | log2-odds PWM at the 7 context positions {-6..-1,+4}; codon columns (+1,+2,+3) masked grey (not defined -- captured by START_W) | `data/kozak_pwm.json` (`pwm` block), fit by `scripts/build_kozak_pwm.py` from 32,332 annotated CDS ATG starts |
| bottom (V0 HEURISTIC) | `kz = 0.5*[purine@-3] + 0.5*[G@+4]`: nonzero only at -3 (A,G = 0.5) and +4 (G = 0.5); codon columns masked | hardcoded in `plot_kozak_weights.py` (`heuristic_matrix()`); matches `build_orf_track.py` ext-mode `kz` |
| suptitle | learned-vs-empirical Pearson r over the 7 shared context positions (per-position mean-centered) | computed in-script from the top + middle panels; also saved to the sibling `.json` |

Sibling `learned_vs_empirical_vs_heuristic.json`: `{run, learned_available, learned_vs_empirical_pearson_r,
bias, kernel_pos, nt_order, learned_kernel_4x10}` -- consumed by `scripts/aggregate_kozak.py` for the Q3
row of `results/kozak/kozak_summary.md`.

Position mapping (conv geometry): `Conv1d` k=10 with pad (6,3) -> output[i] reads one-hot[i-6 .. i+3]; with
the start-codon A at index i, that is Kozak {-6,-5,-4,-3,-2,-1,+1,+2,+3,+4}. One-hot channel order A,C,G,T
is fixed by `dataset.py:onehot_encode`. See methods.md 5N.
