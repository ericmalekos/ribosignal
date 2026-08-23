# References 4 and 5: the actual citations (2026-08-17)

Answering the `[author] et al. [year]` placeholders I flagged in the caption audit. I had NOT sent
these before. They also did not exist anywhere in the project: every dataset was tracked by
accession only, with no source publication recorded. That gap is now closed here.

**Both placeholders bundle three separate studies each.** Neither can be a single reference entry.
Six citations are needed, not two.

Every one was resolved from the accession, not from memory. Provenance for each is stated below.

---

## The six citations

### Reference 4: "Third mouse liver; human THP-1 and CAR-T"

**Third mouse liver (GSE243134)**
> Nadimpalli HP, ... Gatfield D. (2024). Diurnal control of iron responsive element containing
> mRNAs through iron regulatory proteins IRP1 and IRP2 is mediated by feeding rhythms.
> *Genome Biology* 25(1):128. doi:10.1186/s13059-024-03270-2. PMID 38773499.

**Human THP-1 (GSE208041)**
> Ansari SA, ... Uhlenhaut NH. (2022). Integrative analysis of macrophage ribo-Seq and RNA-Seq data
> define glucocorticoid receptor regulated inflammatory response genes into distinct regulatory
> classes. *Computational and Structural Biotechnology Journal* 20:5622-5638.
> doi:10.1016/j.csbj.2022.09.042. PMID 36284713.

**Human CAR-T (GSE304796)**
> Shi Y, ... Sadelain M. (2026). pTalpha enhances mRNA translation and potentiates CAR T cells for
> solid tumor eradication. *Cell* 189(2):401-417.e21. doi:10.1016/j.cell.2025.11.005.

### Reference 5: "Mouse T cell, BMDM, macrophage proteogenomics"

**Mouse CD4+ T cell (GSE155087)**
> Matheson LS, ... Turner M. (2022). Multiomics analysis couples mRNA turnover and translational
> control of glutamine metabolism to the differentiation of the activated CD4+ T cell.
> *Scientific Reports* 12(1):19657. doi:10.1038/s41598-022-24132-6. PMID 36385275.

**Mouse BMDM (GSE120762)**
> Jackson R, ... Flavell RA. (2018). The translation of non-canonical open reading frames controls
> mucosal immunity. *Nature* 564(7736):434-438. doi:10.1038/s41586-018-0794-7. PMID 30542152.

**Macrophage proteogenomics, the MS search data (iProX IPX0001245001 = PRIDE PXD021583)**
> Qie J, ... Ding C. (2022). Integrated proteomic and transcriptomic landscape of macrophages in
> mouse tissues. *Nature Communications* 13:7389. doi:10.1038/s41467-022-35095-7.

---

## Drop-in text

**Recommended, because it needs no renumbering.** It matches the style reference 3 already uses
(`Wang et al. 2021; Janich et al. 2015.`), so reference 6 stays reference 6 and no in-text citation
number moves.

```
4. Nadimpalli et al. 2024; Ansari et al. 2022; Shi et al. 2026. Third mouse liver; human THP-1 and CAR-T.
5. Matheson et al. 2022; Jackson et al. 2018; Qie et al. 2022. Mouse T cell; BMDM; macrophage proteogenomics.
```

**If you would rather give each study its own number,** this is the full list, but note that
RiboCode moves from 6 to 10 and every in-text `[6]` must move with it:

```
4. Nadimpalli et al. 2024. Third mouse liver Ribo-seq.
5. Ansari et al. 2022. Human THP-1 Ribo-seq.
6. Shi et al. 2026. Human CAR-T Ribo-seq.
7. Matheson et al. 2022. Mouse CD4+ T cell Ribo-seq.
8. Jackson et al. 2018. Mouse BMDM Ribo-seq.
9. Qie et al. 2022. Mouse tissue macrophage proteome; the MS search data.
10. Xiao et al. 2018. RiboCode, the ORF caller used throughout.
```

---

## Four things worth knowing before you paste

### 1. The CAR-T paper is 2026, not 2025

Submitted to GEO in August 2025, published in *Cell* in January 2026. If the placeholder year was
going to be guessed from the GEO submission date it would have been wrong.

### 2. "THP-1" and "macrophage" will read as the same system, and they are not

Ansari 2022's own title calls its data **macrophage** Ribo-seq, because THP-1 is a monocytic line
differentiated toward a macrophage phenotype. So reference 4 says "THP-1" and reference 5 says
"macrophage", and a reader scanning the list has no way to tell they are different systems.

They are: reference 4's is a **human cell line**, reference 5's is **mouse primary tissue
macrophages across 12 populations**. Consider "human THP-1 monocytic line" and "mouse tissue
macrophages" to disambiguate.

### 3. Two of the three mouse livers come from the same lab, nine years apart

This one I would want to know before standing next to panel 5.

| mouse liver | study | senior author |
|---|---|---|
| Janich (reference 3) | Genome Research 2015 | **Gatfield D** (Lausanne) |
| GSE243134 (reference 4) | Genome Biology 2024 | **Gatfield D** (Lausanne) |
| Wang (reference 3) | Nucleic Acids Research 2021 | Xie Z (Sun Yat-Sen) |

"Cross-study comparison" is still accurate, since they are three separate studies and three separate
GEO series. But two of the three share a senior author and very likely share protocol lineage, so
they are not three independent labs. Only Wang is fully independent of the other two.

This matters for takeaway 3 specifically: the experiment the model beats on recall is **Wang**,
which is both the shallowest and the only lab-independent one. That is arguably a stronger result
than the caption currently claims, not a weaker one, but do not say "three independent experiments"
if asked.

### 4. Reference 3 checks out

Not asked, but I verified it while I was in there, since it is the same class of claim.

- **Wang et al. 2021** = Wang H, ... Xie Z. *Tissue- and stage-specific landscape of the mouse
  translatome.* Nucleic Acids Res 49(11):6165-6180. doi:10.1093/nar/gkab482. PMID 34107020.
  Confirmed by walking SRR5262890 to GSM2493749 to GSE94982 to the GEO-declared citation.
- **Janich et al. 2015** = Janich P, ... Gatfield D. *Ribosome profiling reveals the rhythmic liver
  translatome and circadian clock regulation by upstream open reading frames.* Genome Research
  25(12):1848-1859. doi:10.1101/gr.195404.115. PMID 26486724.

GSE39561 (THP-1) needs no citation, since that dataset is dropped.

---

## Provenance, and where I am less certain

Four of the six are **GEO-declared**: the GEO series page itself names the PMID, so the
dataset-to-paper link is asserted by the submitters and is not my inference.

| accession | how the paper was resolved | confidence |
|---|---|---|
| GSE208041 | GEO Citation field, PMID 36284713 | certain |
| GSE155087 | GEO Citation field, PMID 36385275 | certain |
| GSE120762 | GEO Citation field, PMID 30542152 | certain |
| GSE94982 (ref 3, Wang) | GEO Citation field, PMID 34107020 | certain |
| **GSE243134** | **no GEO citation.** Matched on a verbatim-identical title plus contributor overlap (Nadimpalli HP, Gatfield D on both) | high, but inferred |
| **GSE304796** | **no GEO citation.** Matched on a verbatim-identical title plus author overlap (Shi Y first, Sadelain lab) | high, but inferred |
| IPX0001245001 / PXD021583 | recorded with its DOI in `proteogenomics/DATA_PROVENANCE.md`; DOI resolved to confirm authors and title | certain |

The two inferred ones are strong matches, the GEO series titles are word-for-word the paper titles,
but they are my inference rather than something the submitters declared. Flagging so you can decide
whether that needs a second pair of eyes before it goes to print.
