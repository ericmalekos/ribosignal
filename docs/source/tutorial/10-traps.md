# Known traps

Every one of these was hit for real during development.

**The adapter is measured, never defaulted.** `SRR15513269` is already trimmed: Illumina TruSeq
appears in 4 reads out of 400,000. Passing `-a` would trim what is not there, and
`--discard-untrimmed` would have destroyed the library.

**RNA-seq and Ribo-seq want opposite trimming.** `--maximum-length` and `--discard-untrimmed` are
correct for footprints and wrong for RNA. On one library the wrong call kept 16.4% of reads where
the right one kept 99.0%.

**STAR can return zero reads and exit 0.** On macOS, bioconda's STAR 2.7.10b accepts
`--readFilesCommand` and silently produces a valid empty BAM. Always assert a non-zero input-read
count from `Log.final.out`.

**`--genomeSAindexNbases` must be scaled for a small genome.** The mammalian default of 14 is
wrong for one chromosome, and STAR builds a poor index rather than erroring.

**`samtools quickcheck`, not a header check.** A BAM truncated mid-write keeps a perfectly valid
header, so grepping the header accepts a broken file.

**Ribo-TISH needs a restricted GTF.** It only skips its BAM path when the supplied profile covers
every transcript of a gene. Passing the full annotation makes it try to read a BAM.

**Coordinate conventions differ by caller and by strand.** RiboCode, Ribo-TISH and RiboTaper
disagree on the ORF start offset per strand, and disagree entirely on whether the stop codon is
inside the ORF. Compare call sets on the start coordinate, and calibrate the offsets rather than
assuming them.

**Callers disagree with each other on real data.** On observed hepatocyte Ribo-seq with no model
involved, three callers agree on 914 canonical ORFs out of 1,160 to 1,852, and on 5 dORFs out of
14 to 25. Do not read a cross-caller count difference as biology.
