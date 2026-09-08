# Call ORFs

RiboCode needs its own prepared annotation once.

```bash
prepare_transcripts -g ref/chr22.gtf -f ref/chr22.fa -o calls/annot
```

Three arms per architecture. `real` runs the caller on the observed P-sites and is the reference
the others are judged against; `pred_preddepth` is the fully de novo call using no Ribo-seq at all;
the Poisson arm trades recall for precision.

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
failure. Its GTF must be restricted to exactly the scored transcripts, because it only skips its
BAM path when the profile covers every transcript of a gene.

```bash
cut -f1 pack/demo/pack_meta.tsv | tail -n +2 > calls/tx_ids.txt
awk -v FL=calls/tx_ids.txt 'BEGIN{while((getline l < FL)>0) k[l]=1}
     !/^#/ { if (match($0,/transcript_id "[^"]+"/)) {
       t=substr($0,RSTART+15,RLENGTH-16); if (t in k) print } }' ref/chr22.gtf > calls/scored.gtf

for ARCH in attn mamba4; do
  python $RIBO_SCRIPTS/orfcallers/ribotish_dropin.py \
         --profiles pred/$ARCH/pred_profiles.npz --gtf calls/scored.gtf --genome ref/chr22.fa \
         --variant pred_preddepth --out calls/ribotish_$ARCH \
         --longest --minaalen 5 --fpth 0.05 --numproc 16
done
```

Pass `--longest`. Without it Ribo-TISH reports every in-frame downstream ATG as a separate
`Truncated` ORF, which was 86% of calls in testing.
