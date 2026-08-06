#!/usr/bin/env bash
# Shared config for held-out Ribo-seq -> transcriptome P-site processing.
# Source with DATASET set. Replicates the Chothani posture-A pipeline
# (expression_context_human phase28 align + phase30 filter + RiboCode) so the
# held-out P-sites are comparable to the training data.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
GENOMES=/private/groups/carpenterlab/emalekos/genomes
RIBOSEQ_ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
RIBOCODE_BIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
HELDOUT_DIR="$NEW/scripts/heldout"

case "${DATASET:?set DATASET}" in
  human_ruizorera)
    SPECIES=human
    FASTQ_DIR="$NEW/data/external/ruizorera2024/fastq/human"
    RIBO_SRRS="ERR12549926 ERR12549927 ERR12549928 ERR12549929 ERR12549930"
    RNA_SRRS="ERR12549891 ERR12549892 ERR12549893 ERR12549894 ERR12549895"  # 5 CM-matched RNAseq (PE)
    ADAPTER_MODE=trimn        # 35 nt N-padded pre-clipped (== Chothani); strip N + len>=20 only
    STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
    ANNOT="$ECH/data/ribocode_annot"
    NCRNA_TX="$ECH/data/phase30/ncrna_filter_tx.txt"
    TX_TO_GENE="$ECH/data/phase30/tx_to_gene.tsv"
    ;;
  mouse_wang_liver)
    SPECIES=mouse
    FASTQ_DIR="$NEW/data/external/wang2021_mouse/fastq/mouse_liver"
    RIBO_SRRS="SRR5262890 SRR5262891"
    RNA_SRRS="SRR5262874 SRR5262875"  # 2 P42 liver mRNAseq (SE); re-fetch via download (deleted earlier)
    ADAPTER_MODE=truseq       # raw TruSeq: -a AGATCGGAAGAGC, footprint 20-40 nt
    STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
    ANNOT="$ECH/data/ribocode_annot_mouse"
    NCRNA_TX="$ECH/data/phase31/mouse_ncrna_tx.txt"
    TX_TO_GENE="$NEW/data/heldout_refs/mouse_tx_to_gene.tsv"
    ;;
  mouse_gse120762_nt)
    SPECIES=mouse
    FASTQ_DIR="$NEW/proteogenomics/data/gse120762/raw"
    RIBO_SRRS="SRR7956050 SRR7956052"                 # RibosomeProfiling CHX_NT rep1,2 (SE, TruSeq)
    RNA_SRRS="SRR7956038 SRR7956039 SRR7956040"       # BMDM WT/NT RNA-seq rep1-3 (PE)
    ADAPTER_MODE=truseq       # raw TruSeq: -a AGATCGGAAGAGC, footprint 20-40 nt (97.7% adapter-bearing)
    STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
    ANNOT="$ECH/data/ribocode_annot_mouse"
    NCRNA_TX="$ECH/data/phase31/mouse_ncrna_tx.txt"
    TX_TO_GENE="$NEW/data/heldout_refs/mouse_tx_to_gene.tsv"
    ;;
  mouse_gse120762_lps)
    SPECIES=mouse
    FASTQ_DIR="$NEW/proteogenomics/data/gse120762/raw"
    RIBO_SRRS="SRR7956051 SRR7956053"                 # RibosomeProfiling CHX_LPS rep1,2 (SE, TruSeq)
    RNA_SRRS="SRR7956041 SRR7956042"                  # BMDM LPS RNA-seq rep1,2 (PE; rep3 SRR7956043 dropped, md5-flaky, 2 reps pool fine)
    ADAPTER_MODE=truseq
    STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
    ANNOT="$ECH/data/ribocode_annot_mouse"
    NCRNA_TX="$ECH/data/phase31/mouse_ncrna_tx.txt"
    TX_TO_GENE="$NEW/data/heldout_refs/mouse_tx_to_gene.tsv"
    ;;
  mouse_gse155087_tcell)
    SPECIES=mouse
    FASTQ_DIR="$NEW/data/heldout_raw/gse155087_tcell/fastq"
    RIBO_SRRS="SRR12318326 SRR12318327 SRR12318328 SRR12318329 SRR12318330"   # WT/control Ribo Ctrl 1-5 (SE, TruSeq)
    RNA_SRRS="SRR12318320 SRR12318321 SRR12318322"                            # WT/control mRNA Ctrl 1-3 (SE)
    ADAPTER_MODE=truseq       # raw TruSeq: -a AGATCGGAAGAGC, footprint 20-40 nt (97.8% adapter-bearing, confirmed)
    STAR_INDEX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
    ANNOT="$ECH/data/ribocode_annot_mouse"
    NCRNA_TX="$ECH/data/phase31/mouse_ncrna_tx.txt"
    TX_TO_GENE="$NEW/data/heldout_refs/mouse_tx_to_gene.tsv"
    ;;
  *) echo "unknown DATASET: $DATASET" >&2; exit 1 ;;
esac

BAM_DIR="$NEW/data/heldout_bam/$DATASET"
PSITE_DIR="$NEW/data/heldout_psites/$DATASET"
LOG_DIR="$NEW/logs/heldout"
