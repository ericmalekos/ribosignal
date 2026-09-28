# Call ORFs

RiboCode needs its own prepared annotation once.

```bash
mkdir -p calls
prepare_transcripts -g ref/chr22.gtf -f ref/chr22.fa -o calls/annot
```

Three arms per architecture. `real` runs the caller on the observed P-sites and is the reference
the others are judged against; `pred_preddepth` calls from the predicted shape and the predicted
depth, so no observed P-site enters the call (the transcripts it runs on are the ones chosen on the
predict page); the Poisson arm trades recall for precision.

```bash
for ARCH in attn mamba4; do
  P=pred/$ARCH/pred_profiles.npz
  python $RIBO_SCRIPTS/ribocode_dropin.py --profiles "$P" --annot calls/annot \
         --variant real --out calls/$ARCH --min_aa 30 --pval 0.05
  python $RIBO_SCRIPTS/ribocode_dropin.py --profiles "$P" --annot calls/annot \
         --variant pred_preddepth --out calls/$ARCH --min_aa 30 --pval 0.05
  python $RIBO_SCRIPTS/ribocode_dropin.py --profiles "$P" --annot calls/annot \
         --variant pred_preddepth --pred_poisson --pred_scale 0.05 \
         --out calls/${ARCH}_poisson --min_aa 30 --pval 0.05
done
```

Ribo-TISH is a second caller on the same profiles, which separates a model failure from a caller
failure. Its GTF must be restricted to the predicted transcripts, keeping only genes whose every
transcript was predicted, because Ribo-TISH skips its BAM path only when the profile covers a
whole gene. `--emit-gtf` makes the wrapper write that GTF itself. Building it from
`pack_meta.tsv` instead is wrong: the pack holds every transcript, the prediction only the scored
ones, and Ribo-TISH then goes looking for a BAM that is not there.

```bash
for ARCH in attn mamba4; do
  python $RIBO_SCRIPTS/orfcallers/ribotish_dropin.py \
         --profiles pred/$ARCH/pred_profiles.npz --gtf ref/chr22.gtf --genome ref/chr22.fa \
         --emit-gtf calls/scored_$ARCH.gtf \
         --variant pred_preddepth --out calls/ribotish_$ARCH \
         --longest --minaalen 5 --fpth 0.05 --numproc 16
done
```

Pass `--longest`. Without it Ribo-TISH reports every in-frame downstream ATG as a separate
`Truncated` ORF, which was 86% of calls in testing.
