# X_gtf2i_isoforms -- GTF2I observed Ribo-seq across all 35 isoforms, held-out hepatocytes

**Status:** exploratory (the `X_` prefix marks it as not part of the numbered figure set).
Generated 2026-08-19 in response to a direct question; kept because it is the clearest worked
example of the posture-A multimapper convention in the project.

**What it shows:** GTF2I has 35 annotated isoforms in the union universe. Under posture A
(methods.md 2.3) a footprint is counted at weight 1 on *every* compatible isoform of its gene, so
one gene's ~13,000 footprints are reported 35 times over.

## Panel -> data

| element | source |
|---|---|
| per-nt observed P-sites, panel (a) | `data/packed_union_Hepatocytes/target_counts.npy`, sliced by `offsets.npy` / `lengths.npy`, rows keyed by `tx_order.txt` |
| isoform identity, length, RNA TPM | `data/union_universe.tsv` (`gene_name == "GTF2I"`, cols `tx_id`, `length`, `max_tpm`) |
| totals + spread, panel (b) | summed from the same `target_counts.npy` slices |

Hepatocytes is the **LOTO held-out tissue**. `packed_union_Hepatocytes` is on the 84,472-transcript
union universe, which carries all 35 isoforms; `packed_Hepatocytes` is the smaller 36,668-tx
universe and carries only **6** of them, so it is the wrong pack for this figure.

## The numbers

- 35 isoforms, **452,992** summed P-sites, **34.2x** inflation over any single isoform.
- Per-isoform totals span 11,161 to 13,258: a spread of only ~16% of the max.
- RNA TPM spans 1.11 to 50.31, a **45x** range.
- The isoform with the MOST observed P-sites (`ENST00000901263.1`, 13,258) has **TPM 1.25**. The
  highest-RNA isoform (`ENST00000620879.4`, TPM 50.31) has 13,138 -- slightly fewer.

## Why this matters

**Observed Ribo-seq carries almost no information about which isoform is expressed.** RNA spans
45x across these isoforms while observed footprints span 16%, and the rank order is close to
uncorrelated. This is the direct empirical argument against selecting a gene's representative
isoform by Ribo-seq: at the isoform level the Ribo-seq signal is nearly constant by construction,
because posture A copies the same reads onto every compatible isoform. `make_representative_tx.py`
selects on **RNA TPM** for exactly this reason.

It is also why the ~18x project-wide inflation figure is a *per-gene-family* average, not a
constant: GTF2I sits at 34.2x because it is unusually isoform-rich.

## Regenerating

```
$RIBOSEQ/conda_envs/riboseq/bin/python3 figures/X_gtf2i_isoforms/make_gtf2i_isoforms.py
```
Needs h5py-free numpy + matplotlib only. Writes `X_gtf2i_isoforms.{png,pdf}` and
`X_gtf2i_values.json` (`model: null` -- observed data only, no prediction).

---

## Added 2026-08-20: two further views

| figure | script | what it adds |
|---|---|---|
| `X_gtf2i_frames.{png,pdf}` | `make_gtf2i_frames.py` | one panel per isoform, P-sites coloured by frame **relative to that isoform's own CDS start** (`data/human_ribocode_annot_primary/transcripts_cds.txt`). `--top N --ncol N --tag S`. |
| `X_gtf2i_frames_top10.*` | same, `--top 10 --ncol 1` | top 10 by RNA TPM, one full-width panel each |
| `X_gtf2i_model_top10.*` | `make_gtf2i_model.py --top 10` | adds each isoform's exon structure below its profile, in the same transcript coordinates |

**Frame convention:** for 0-based position `i` and 1-based `cds_start`,
`frame = (i - (cds_start-1)) % 3`. Frame 0 is in-frame with that isoform's CDS. Because each
isoform has a different `cds_start`, identical shared footprints receive different frame labels
across isoforms -- that is the point, not a bug.

**Gene-model track:** the x-axis is TRANSCRIPT position, so exons tile end to end and no introns
appear. What the track carries is where the JUNCTIONS fall. Junctions are cumulative exon lengths
from the GENCODE v49 primary GTF ordered by `exon_number`; **verified that exon lengths sum exactly
to the packed transcript length for all 35 isoforms.**

**Findings these two added**, both recorded in `docs/isoform_selection_by_riboseq.md`:
frame-0 fraction spans only 85.3-85.8% across all 35 isoforms, so periodicity does not discriminate
them either; and the CDS starts differ by hundreds of nt but are near-multiples of 3, leaving the
isoforms in register with one another.

Values files: `X_gtf2i_frames{,_top10}_values.json`, `X_gtf2i_model_top10_values.json`, all
`model: null`.

---

## Superseded in part, 2026-08-20

These figures use the pack's RNA `coverage.npy` and `union_universe.tsv` `max_tpm` as the abundance
axis. **Both are now known to be unsound for isoform-level abundance**: pack coverage is
isoform-multimapped and length-biased (Spearman(length, cov/nt) ~ -0.8, and |rho| <= 0.17 against
salmon TPM in all 8 tissues), and `max_tpm` is a cross-tissue maximum that inverts within-tissue
rankings.

The sound per-tissue values now live at `data/tpm/chothani_alnmode/<tissue>_mean_tpm.tsv`, with the
GTF2I summary at `data/tpm/chothani_alnmode/GTF2I_per_tissue.json`. **The figures themselves remain
valid** -- their point is that Ribo-seq cannot discriminate isoforms, which the new data confirms
rather than contradicts. Only the RNA axis should be re-read against the new tables.
