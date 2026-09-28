# RiboSignal

Predicts a per-nucleotide ribosome P-site profile for a transcript from its sequence and matched
RNA-seq coverage. It takes any transcript, coding or non-coding, up to 10,000 nt, the longest the
checkpoints were trained on. No ribosome-profiling experiment is needed at inference.

The tutorial below runs the whole pipeline on public data, on chromosome 22, and ends with a
number: predicted ORF calls scored against real Ribo-seq from the same donor. The Ribo-seq is there
to score the prediction. The Nextflow pipeline takes RNA-seq alone, and the pack and predict
pages show how to run the scripts without Ribo-seq.

```{toctree}
:maxdepth: 1
:caption: Tutorial

tutorial/01-environment
tutorial/02-annotation
tutorial/03-sequencing-data
tutorial/04-adapters
tutorial/05-alignment
tutorial/06-pack
tutorial/07-predict
tutorial/08-call-orfs
tutorial/09-score
tutorial/11-cluster
tutorial/12-nextflow
```
