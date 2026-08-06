# RiboSignal: predicting ribosome profiling signal from sequence and RNA-seq

## Introduction

Ribosome profiling (Ribo-seq) maps translating ribosomes across the transcriptome and is the
standard way to find open reading frames that genome annotation misses, including upstream ORFs in
5' untranslated regions and ORFs on long non-coding RNAs. Peptides from these non-canonical ORFs
have roles in stress response, immunity, and disease, and several are candidate immunotherapy
targets. The obstacle is practical. Ribo-seq is technically demanding, needs large amounts of fresh
material, and exists for only a small number of cell types. RNA-seq is routine and available for
nearly every tissue and condition. If the per-nucleotide ribosome signal could be predicted from
sequence and RNA-seq coverage, ORF discovery would become possible wherever expression data already
exist. We asked whether a neural network can learn that mapping well enough for a standard ORF
caller to run on predicted signal in place of measured signal.

## Materials and Methods

The model is a dilated convolutional encoder followed by a sequence mixer, either transformer layers
(5.1 million parameters, CPU-runnable) or Bidirectional Mamba blocks (7.5 million, GPU only). It
takes three per-nucleotide inputs: one-hot sequence, RNA-seq coverage, and a track marking candidate
ORF starts and stops. It outputs a per-nucleotide P-site distribution and a predicted total count.
We trained on seven human cell types from a published Ribo-seq compendium and held out hepatocytes
entirely. Predicted profiles were then passed to RiboCode, an established ORF caller, in place of
experimental counts. A calibration step samples the predicted profile as Poisson counts at a scaled
depth, tuned so that recall of annotated coding sequences reaches a set target. We applied the same
model to mouse macrophages and to four human mass spectrometry datasets, and tested the predicted
ORFs against matched proteomes and MHC class I immunopeptidomes.

## Results and Discussion

On held out hepatocytes the calibrated model recovers annotated coding sequences at F1 0.90
(precision 0.911, recall 0.898), upstream ORFs at precision 0.766, and non-canonical ORFs at
precision 0.591. Calibration is the decisive control. Without it, precision on non-canonical ORFs
drops to 0.381 while recall rises, so the dial trades yield for confidence at fixed canonical
recall. Applied to mouse macrophages, a tissue and species the model never saw, predicted ORFs gave
a 1,533 sequence proteomic search database that identified 37 novel peptides while costing only 36
canonical spectrum matches. Enumerating every near-cognate ORF instead produces 2.5 million
sequences and 166 peptides, but destroys 29,633 canonical spectrum matches through false discovery
rate inflation. The pattern holds on a deep human proteome: a 2,404 sequence predicted database
found 11 novel peptides at zero canonical cost, where the 3.9 million sequence near-cognate
enumeration found 38 and destroyed 11,278. Across four human datasets the two sequence mixers were
indistinguishable for this task, 49 against 48 novel peptides, so the smaller CPU-runnable model
loses nothing in discovery.

## Conclusions

Predicted ribosome profiles are accurate enough to stand in for a ribosome profiling experiment
when the goal is finding open reading frames. Any tissue with RNA-seq can now be searched for
non-canonical translation, and because the predicted ORF sets are small and evidence based, they
work as mass spectrometry search databases without sacrificing canonical protein identifications.
The method needs no ribosome profiling, no GPU, and no foundation model embeddings. Ongoing work
extends it to additional immunopeptidomes and to N-terminal extensions initiating at non-AUG start
codons.
