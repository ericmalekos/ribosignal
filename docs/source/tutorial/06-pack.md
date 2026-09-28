# Build the pack and the ORF track

A pack is the model's input format: one row per transcript, with coverage and P-sites on the same
axis as the sequence.

```bash
python $RIBO_SCRIPTS/build_pack.py --fasta ref/transcripts.fa \
       --coverage pack/coverage.hd5 --psites pack/psites.hd5 --out pack/demo
```

Every transcript in the FASTA goes in, whatever its biotype, provided the BAMs were aligned against
the same annotation: the pack joins on transcript id and checks lengths. The one exception is
length. Transcripts longer than 10,000 nt are dropped, because the checkpoints never saw one and
attention memory grows with the square of the length. On chr22 that keeps 11,575 of 11,615.

`--psites` is optional. For a sample with no Ribo-seq, leave it out: the pack's P-site target is
then all zeros, and `build_pack.py` also writes `expressed_tx.txt`, the transcripts with RNA-seq
coverage, which the predict page uses instead.

The ORF track marks candidate ORF positions and reading frames from sequence alone. `--kozak none`
is the setting both released checkpoints were trained with.

```bash
python $RIBO_SCRIPTS/build_orf_track.py --pack pack/demo --fasta ref/transcripts.fa \
       --mode ext --kozak none

export RIBO_ONEHOT_FASTA=$PWD/ref/transcripts.fa
export RIBO_ORF_TRACK=$PWD/pack/demo/orf_track_v2_nokozak.npy
```
