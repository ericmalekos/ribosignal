# Score against the real Ribo-seq

The run ends with a number rather than an assertion that it completed. `score_demo.py` reports
per-nucleotide profile correlation and ORF-call precision, recall and F1, predicted against
observed.

```bash
python $RIBO_SCRIPTS/demo/score_demo.py --pred-dir pred --calls-dir calls \
       --out RESULTS.json | tee RESULTS.txt
```

Read three things. **Profile correlation** is how well the predicted shape matches the observed
one. **Annotated-CDS F1** should be high; this is the class the annotation makes trustworthy.
**Non-canonical F1** is the hard case, and the uncalibrated arm over-calls it. Comparing the plain
`pred_preddepth` arm against the Poisson arm shows what the dial buys.
