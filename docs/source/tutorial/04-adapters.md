# Measure the adapters before trimming

Adapters are measured, never assumed. `SRR15513269` is already adapter-trimmed, so passing `-a`
would trim sequence that is not there, and `--discard-untrimmed` would destroy the library.

```bash
python $RIBO_SCRIPTS/demo/measure_adapter.py fastq/SRR15513269_1.fastq.gz | tee logs/adapter_rna.txt
python $RIBO_SCRIPTS/demo/measure_adapter.py fastq/SRR15513208.fastq.gz   | tee logs/adapter_ribo.txt
```

The RNA read should report `already trimmed`. Read the Ribo verdict out of the file rather than
hardcoding it.

```bash
grep -q "already trimmed" logs/adapter_rna.txt || echo "RNA verdict changed: read the file"
RIBO_ADAPTER=$(awk '/Use `cutadapt -a /{print $4}' logs/adapter_ribo.txt | head -1)
echo "ribo adapter: ${RIBO_ADAPTER:-<none detected>}"
```
