## Preprocessing
This section was done by Mingshen Shen and Bomin Lee jointly. All preprocessing work was conducted on the King’s College London (KCL) HPC cluster.

### Cell types
- Nurr (dopaminergic neurons)
- Olig (oligodendrocytes)  
- NeuN (neurons)  
- Neg (microglia-enriched)

### Rat reference genome
mRatBN7.2

### Required tools
- FastQC (v0.11.5) — quality control  
- Trim Galore (v0.6.10) — trimming  
- MultiQC (v1.15) — generating a summary report for all samples  
- Bowtie2 (v2.5.4) — alignment  
- Picard (v2.18.29) — removing duplicates  
- Samtools (v1.22) — merging replicates  
- MACS2 (v2.2.7.1) — peak calling  
- Bedtools (v2.31.0) — calculating FRiP scores

## 1. Cell Type-Specific Regulatory Landscapes
This subproject contains comprehensive cell type-specific enrichment analyses in a rat rotenone model of Parkinson’s disease (PD). The analyses integrate genomic, epigenomic, and statistical genetics approaches to identify which brain cell types show the strongest genetic contribution to PD pathogenesis.

### Aims
- To integrate cell type–specific H3K27ac enhancer landscapes from rat midbrain.
- To combine these enhancer annotations with Parkinson’s disease GWAS summary statistics.
- To quantify cell type–specific PD heritability enrichment using LD Score Regression (LDSC).

### Key Findings
- **Dopaminergic enhancers** show the strongest enrichment of PD heritability.  
- **Oligodendrocyte enhancers** display moderate to strong PD heritability enrichment.  
- Enhancer annotations for other two cell types exhibit lower or non-significant enrichment.  
- These results support a **multi–cell-type genetic architecture** underlying Parkinson’s disease risk.

### Required packages or tools
**R Environment**
- R 4.4.2

**Core R Packages**
- ChIPseeker 1.42.1
- clusterProfiler 4.14.6
- org.Rn.eg.db 3.20.0
- TxDb.Rnorvegicus.UCSC.rn7.refGene 3.15.0
- enrichplot 1.26.6
- ggplot2 3.5.2

**External Tools**
- Bedtools v2.31.0 — for peak processing  
- HOMER v4.x — motif analysis  
- LDSC — Python3 fork: abdenlab/ldsc  


## 2. Deconvolution analysis using CHAS
This subproject deconvolutes the same bulk H3K27ac data using CHAS, a deconvolution tool. So that we can determine which DARs were specific to individual cell types. Furthermore, this approach investigated cell type–specific gene regulatory changes and dysregulated pathways across the substantia nigra (SN) using GO analysis. 

### Aims
- To deconvolute a bulk H3K27ac ChIP-seq datasets using CHAS.
- To identify cell type–specific differentially acetylated regions (DARs) across the substantia nigra (SN).
- To investigate dysregulated gene regulatory pathways and biological processes through GO enrichment analysis.

### Key Findings
- **Oligodendrocytes** show a consistent pattern of **hypoacetylation** in the substantia nigra (SN).
- Hypoacetylated regions in oligodendrocytes are strongly associated with **myelin-related functions**, suggesting widespread disruption of oligodendroglial regulation.
- The possible reduction in microglia may be linked to **oxidative stress–induced apoptotic pathways**.

### Required packages or tools
**R Environment**
- R 4.4.3

**Core R Packages**
- dplyr 1.1.4
- edgeR 4.4.2
- pheatmap 1.0.13
- clusterProfiler 4.14.6
- org.Rn.eg.db 3.20.0
- ggplot2 3.5.2
- ggrepel 0.9.6

**External Tools**
- Subread v2.0.6 — for constructing matrix for CHAS input  

## Data Access
The raw and processed data used in this project are stored on the King’s College London (KCL) High Performance Computing (HPC) system. Due to institutional data governance policies, these datasets are not publicly available.

Researchers who require access to the data for collaborative purposes may contact the project contributors listed below.

## Contributors
- **Bomin Lee** — Responsible for preprocessing and the analysis of Cell Type-Specific Regulatory Landscapes.
- **Mingshen Shen** — Responsible for preprocessing and Deconvolution analysis using CHAS.
- **Maria Tsalenchuk** — Provided the raw ChIP-seq datasets used in this project.
- **Paulina Urbanaviciute**, **Sarah J. Marzi** — Supervision and project guidance.

## Reference

This project uses published enhancer annotations, midbrain cell-type datasets, and statistical genetics methods as described in:

- **Nott A et al.** Cell type–specific enhancer–promoter interactome maps of the human brain. *Science*. 2019; 366(6469):1134–1139.

- **Agarwal D et al.** A single-cell atlas of the human midbrain identifies cell-specific genetic risk for Parkinson's disease. *Nature Neuroscience*. 2020; 23:939–951.

This project makes use of the CHAS deconvolution framework as described in:

- **Murphy KB et al.** *CHAS infers cell type–specific signatures in bulk brain histone acetylation studies of neurological and psychiatric disorders.* Cell Reports Methods. 2025; 5.

The rotenone-exposed rat H3K27ac data used for this project were originally reported in:

- **Tsalenchuk M et al.** *Unique nigral and cortical pathways implicated by epigenomic and transcriptional analyses in a rotenone Parkinson’s model.* npj Parkinson's Disease. 2025; 11: 217.


