# Rat Midbrain ChIP-seq Analysis - Local Setup Guide

This guide explains how to run the ChIP-seq analysis pipeline in a local environment.

## 📋 Table of Contents
1. [Environment Setup](#environment-setup)
2. [Data Preparation](#data-preparation)
3. [Pipeline Execution](#pipeline-execution)
4. [Script Descriptions](#script-descriptions)
5. [Troubleshooting](#troubleshooting)

---

## 🔧 Environment Setup

### 1. Create Conda Environment

```bash
# Navigate to project directory
cd projects/rat-midbrain

# Create Conda environment (run only once)
conda env create -f environment.yml

# Activate environment
conda activate ratlas_env
```

### 2. Download Reference Genome (if needed)

#### Download Rat genome (rn7):
```bash
# Download genome FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz
gunzip rn7.fa.gz

# Build Bowtie2 index (if alignment is needed)
bowtie2-build rn7.fa rn7_index/rn7
```

#### Edit config.yaml:
```yaml
genome:
  fasta: "/path/to/rn7.fa"
  bowtie2_index: "/path/to/rn7_index/rn7"
```

---

## 📁 Data Preparation

### 1. Place FASTQ Files

Place FASTQ files in the `0_data/` folder:

```bash
projects/rat-midbrain/0_data/
├── IGF131357_R1.fastq.gz
├── IGF131357_R2.fastq.gz
├── IGF131358_R1.fastq.gz
├── IGF131358_R2.fastq.gz
├── ...
└── IGF131377_R2.fastq.gz
```

### 2. Sample Information

**Sample composition by cell type:**
- **NeuN** (Neurons): IGF131357, IGF131358, IGF131359
  - Input: IGF131373
- **Nurr** (Dopaminergic neurons): IGF131366, IGF131367
  - Input: IGF131376
- **Olig** (Oligodendrocytes): IGF131360, IGF131361, IGF131362
  - Input: IGF131374
- **Neg** (Negative control): IGF131369, IGF131370, IGF131371
  - Input: IGF131377

---

## 🚀 Pipeline Execution

### Method 1: Run Individual Scripts (Recommended)

Run each step sequentially:

```bash
# Activate environment
conda activate ratlas_env

# Step 1: Preprocessing (requires FASTQ files)
bash scripts/01_preprocessing_qc_local.sh

# Step 2: Alignment (requires reference genome)
bash scripts/02_alignment.sh

# Step 3: Peak Calling (requires BAM files)
bash scripts/03_peak_calling_local.sh

# Step 4: QC Metrics (FRiP calculation)
bash scripts/04_qc_metrics_local.sh

# Step 5: Peak Processing
bash scripts/05_peak_processing.sh

# Step 6: Annotation (requires R)
Rscript scripts/06_annotation.R

# Step 7: Enrichment Analysis (requires R)
Rscript scripts/07_enrichment.R

# Step 8: Motif Analysis
bash scripts/08_motif_analysis.sh
```

### Method 2: Run Python Main Pipeline

```bash
# Run full pipeline
python main.py

# Start from specific step
python main.py --start-from peak_calling

# Stop at specific step
python main.py --stop-at qc_metrics

# Skip specific steps
python main.py --skip preprocessing alignment

# List available steps
python main.py --list-steps
```

---

## 📝 Script Descriptions

### Local Version Scripts (No SLURM Required)

#### `01_preprocessing_qc_local.sh`
- **Function**: FASTQ QC and adapter trimming
- **Input**: `0_data/*_R1.fastq.gz`, `*_R2.fastq.gz`
- **Output**: `01_preprocessing_qc/trimmed_fastq/`
- **Differences from original**:
  - SLURM array job → for loop
  - Removed `module load` (uses conda)
  - Uses local paths

#### `03_peak_calling_local.sh`
- **Function**: Peak calling by cell type (merged replicates)
- **Input**: `02_alignment/filtered_bam/*.bam`
- **Output**: `03_peak_calling/macs2_output/`
- **Differences from original**:
  - Only paths modified to local
  - MACS2 parameters unchanged

#### `04_qc_metrics_local.sh`
- **Function**: FRiP score calculation (fragment-level)
- **Input**: BAM files + Peak files
- **Output**: `04_qc_metrics/frip_scores.csv`
- **Differences from original**:
  - Only paths modified to local
  - Logic completely identical

---

## ⚠️ Troubleshooting

### 1. "command not found" Error

```bash
# Check if Conda environment is activated
conda activate ratlas_env

# Verify specific tool installation
which trim_galore
which macs2
which bedtools
```

### 2. Out of Memory

```bash
# If memory error occurs in MACS2, run samples separately
# Or adjust --buffer-size option
```

### 3. Bowtie2 Index Error

```bash
# Set BOWTIE2_INDEX environment variable
export BOWTIE2_INDEX=/path/to/rn7_index/rn7

# Or modify directly in script
# Edit BOWTIE2_INDEX variable in scripts/02_alignment.sh
```

### 4. R Package Installation Issues

```R
# Manual installation in R console
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "ChIPseeker",
    "clusterProfiler",
    "TxDb.Rnorvegicus.UCSC.rn7.refGene",
    "org.Rn.eg.db"
))
```

### 5. HOMER Genome Installation

```bash
# Install HOMER genome (required for motif analysis)
perl $(which configureHomer.pl) -install rn7

# Or check latest version
homer2 listGenomes
```

---

## 📊 Results

### Main Output Files:

```
01_preprocessing_qc/
├── multiqc_report.html          ← Overall QC summary
└── trimmed_fastq/               ← Trimmed FASTQ files

03_peak_calling/
└── macs2_output/
    ├── NeuN/NeuN_all_peaks.narrowPeak
    ├── Nurr/Nurr_all_peaks.narrowPeak
    ├── Olig/Olig_all_peaks.narrowPeak
    └── Neg/Neg_all_peaks.narrowPeak

04_qc_metrics/
├── frip_scores.csv              ← FRiP score results
└── qc_summary.txt               ← QC summary

06_annotation/
├── annotation_summary.csv       ← Peak annotation
└── plots/                       ← Visualizations

07_enrichment/
├── results/                     ← GO/KEGG results
└── figures/                     ← Enrichment plots

08_motif_analysis/
└── motif_results/               ← HOMER motif results
    └── */homerResults.html      ← Results HTML
```

---

## 🎯 Quick Start (If You Already Have BAM Files)

If you already have BAM files from completed alignment:

```bash
# 1. Place BAM files in correct location
cp /path/to/bams/*.bam projects/rat-midbrain/02_alignment/filtered_bam/

# 2. Start from peak calling
conda activate ratlas_env
bash scripts/03_peak_calling_local.sh
bash scripts/04_qc_metrics_local.sh
bash scripts/05_peak_processing.sh
Rscript scripts/06_annotation.R
Rscript scripts/07_enrichment.R
bash scripts/08_motif_analysis.sh
```

---

## 📧 Contact

If you encounter problems or have questions:
- Create GitHub Issues
- Email: bomin.lee@kcl.ac.uk

---

## 📚 References

- [MACS2 Documentation](https://github.com/macs3-project/MACS)
- [ChIPseeker Documentation](https://bioconductor.org/packages/ChIPseeker)
- [HOMER Documentation](http://homer.ucsd.edu/homer/)
- [Trim Galore Documentation](https://github.com/FelixKrueger/TrimGalore)
