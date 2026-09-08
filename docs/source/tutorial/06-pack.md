# Build the pack and the ORF track

A pack is the model's input format: one row per transcript, with coverage and P-sites on the same
axis as the sequence.

```bash
python $RIBO_SCRIPTS/build_pack.py --fasta ref/transcripts.fa \
       --coverage pack/coverage.hd5 --psites pack/psites.hd5 --out pack/demo
```

The ORF track marks candidate ORF positions and reading frames from sequence alone. `--kozak none`
is the setting both released checkpoints were trained with.

```bash
python $RIBO_SCRIPTS/build_orf_track.py --pack pack/demo --fasta ref/transcripts.fa \
       --mode ext --kozak none

export RIBO_ONEHOT_FASTA=$PWD/ref/transcripts.fa
export RIBO_ORF_TRACK=$PWD/pack/demo/orf_track_v2_nokozak.npy
```
