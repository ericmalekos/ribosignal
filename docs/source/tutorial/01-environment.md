# Install and check the environment

Everything runs inside one container image, which carries the aligners, the callers, the Python
stack and both released checkpoints. Nothing is downloaded at run time.

```bash
docker pull ghcr.io/ericmalekos/riboseq-model:latest
mkdir -p ~/ribosignal_demo && cd ~/ribosignal_demo
docker run --rm -it -v "$PWD:$PWD" -w "$PWD" ghcr.io/ericmalekos/riboseq-model:latest bash
```

Check the tools are present. The build already asserts this, so a failure here means a bad pull.

```bash
STAR --version && samtools --version | head -1 && cutadapt --version
gffread --version && ribotish --version && ribotricer --version
python -c "import torch, numpy, pysam, h5py, RiboCode; print('ok', torch.__version__)"
ls $RIBO_WEIGHTS
```

Two notes. **Salmon is in the image but this tutorial does not use it**; it is there for building
an expressed-transcript universe, which a chromosome subset makes unnecessary. **RiboTaper is not
in the image**: it pins Python 2.7 and is GPL-3, so install it separately if you want it.

Run on Linux. On macOS the bioconda STAR builds either abort on `--quantMode` or silently return
zero reads.
