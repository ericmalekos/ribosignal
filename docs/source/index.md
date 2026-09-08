# RiboSignal

Predicts a per-nucleotide ribosome P-site profile for a transcript from its mature mRNA sequence
and matched RNA-seq coverage. No ribosome-profiling experiment is needed at inference.

The tutorial below runs the whole pipeline on public data, on chromosome 22, and ends with a
number: predicted ORF calls scored against real Ribo-seq from the same donor.

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
tutorial/10-traps
tutorial/11-cluster
```
