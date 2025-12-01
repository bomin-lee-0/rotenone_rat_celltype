# LDSC Cell-Type Specific Enhancer Enrichment Analysis Workflow

## 📋 Overview

This document summarizes the complete workflow for analyzing heritability enrichment in cell type-specific enhancer regions from Parkinson's disease GWAS data using LDSC (Linkage Disequilibrium Score Regression).

## 🎯 Research Objectives

- **Main Objective**: Quantify the extent to which enhancer regions of each brain cell type contribute to Parkinson's disease heritability
- **Hypothesis**: Parkinson's disease-related variants will be enriched in enhancer regions of specific brain cell types
- **Expected Results**: Confirm enrichment values and statistical significance for each cell type

## 📊 Datasets

### Input Data
1. **GWAS Summary Statistics**: Parkinson's disease GWAS data (`parkinson_gwas.sumstats.gz`)
2. **Cell Type-Specific Enhancer BED Files**: 8 datasets
   - Neg_cleaned, Neg_unique
   - NeuN_cleaned, NeuN_unique
   - Nurr_cleaned, Nurr_unique
   - Olig_cleaned, Olig_unique
3. **Reference Panel**: BaselineLD v2.2 (97 functional annotations)
4. **1000G EUR Phase 3**: LD calculation and weight files

### Data Structure
```
/scratch/prj/eng_waste_to_protein/repositories/bomin/
├── 0_data/
│   ├── raw/cleaned_data/           # Original BED files
│   └── reference/ldsc_reference/   # BaselineLD, 1000G files
├── ldsc_results/
│   ├── annotations/               # Generated annotation files
│   ├── simple_ld_scores/         # LD score files
│   └── sumstats/                 # Preprocessed GWAS data
└── ldsc_results_final/           # Final result files
```

## 🔬 Methodology

### 1. LDSC Theoretical Background

LDSC uses the following linear regression model:

```
E[χ²ⱼ] = 1 + N × Σₖ τₖ × lₖⱼ
```

**Variable Descriptions:**
- `χ²ⱼ`: GWAS chi-square statistic of the j-th SNP
- `N`: GWAS sample size
- `τₖ`: Per-SNP heritability coefficient of the k-th annotation
- `lₖⱼ`: LD score of the j-th SNP for the k-th annotation

### 2. Enrichment Calculation Formula

```
Enrichment = (Proportion of h²g) / (Proportion of SNPs)
```

**Interpretation:**
- Ratio of the heritability proportion explained by the annotation to the SNP proportion
- Greater than 1 indicates a more important region than average (enriched)
- Less than 1 indicates a less important region than average (depleted)

### 3. Statistical Significance Testing

**Z-score Calculation:**
```
Z = coefficient / coefficient_SE
```

**p-value Calculation:**
```
p-value = 2 × (1 - Φ(|Z|))
```
Where Φ is the cumulative distribution function of the standard normal distribution

## 🛠 Step-by-Step Workflow

### Step 1: Environment Setup and Data Preparation

```bash
# Set working directory
cd /scratch/prj/eng_waste_to_protein/repositories/bomin

# Verify LDSC Python3 version
cd 1_preprocessing/ldsc-python3
python ldsc.py --help
```

### Step 2: BaselineLD File Validation

```python
def verify_baseline_files():
    """Validate BaselineLD file format and count"""
    ref_dir = Path('0_data/reference/ldsc_reference')

    # Check file existence
    assert len(list(ref_dir.glob('baselineLD.*.l2.ldscore.gz'))) == 22
    assert len(list(ref_dir.glob('baselineLD.*.l2.M'))) == 22

    # Validate format with first chromosome
    with gzip.open(ref_dir / 'baselineLD.1.l2.ldscore.gz', 'rt') as f:
        header = f.readline().strip().split('\t')
        n_annot = len(header) - 3  # Exclude CHR, SNP, BP

    # Verify 97 BaselineLD annotations
    assert n_annot == 97
    return n_annot
```

### Step 3: Cell Type-Specific Annotation Generation

**Using Pre-generated Annotation Files** (coordinate issue resolved):

```python
def use_existing_annotations(dataset_name):
    """Use pre-generated annotation files"""
    annot_dir = Path('ldsc_results/annotations')

    # Verify 22 chromosome annotation files
    for chr_num in range(1, 23):
        annot_file = annot_dir / f'{dataset_name}.{chr_num}.annot.gz'
        assert annot_file.exists()

    return True
```

**Annotation File Structure:**
```
CHR  BP    SNP         BaselineLD_annotations...  CellType_enhancer
1    11008 rs575272151 [97 columns]              0
1    15274 rs62635286  [97 columns]              1
```

### Step 4: LD Score Calculation

**Using Pre-calculated LD Score Files** (already computed):

```bash
# Verify LD score files
ls ldsc_results/simple_ld_scores/Neg_cleaned.*.l2.ldscore.gz
ls ldsc_results/simple_ld_scores/Neg_cleaned.*.l2.M
ls ldsc_results/simple_ld_scores/Neg_cleaned.*.l2.M_5_50
```

### Step 5: LDSC Regression Execution

```python
def run_ldsc_regression(dataset_name):
    """Run LDSC regression"""

    # Set file paths
    gwas_file = 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz'
    ref_ld_prefix = f'ldsc_results/simple_ld_scores/{dataset_name}.'
    w_ld_prefix = '0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
    frq_prefix = '0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.'
    output_prefix = f'ldsc_results_final/{dataset_name}_h2'

    # Execute LDSC command
    cmd = [
        'python', 'ldsc.py', '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]

    result = subprocess.run(cmd, cwd='1_preprocessing/ldsc-python3')
    return output_prefix + '.log' if result.returncode == 0 else None
```

### Step 6: Result Extraction and Analysis

```python
def extract_enrichment_results(log_file, dataset_name):
    """Extract enrichment information from LDSC results"""

    with open(log_file, 'r') as f:
        content = f.read()

    results = {}

    # Extract total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        results['total_h2'] = float(h2_match.group(1))
        results['total_h2_se'] = float(h2_match.group(2))

    # Extract enrichment values
    lines = content.split('\n')
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_values = line.replace('Enrichment:', '').strip().split()
            if len(enrichment_values) >= 2:
                results['base_enrichment'] = float(enrichment_values[0])
                results['enhancer_enrichment'] = float(enrichment_values[1])

    # Extract proportion information
    for line in lines:
        if line.startswith('Proportion of h2g:'):
            prop_values = line.replace('Proportion of h2g:', '').strip().split()
            if len(prop_values) >= 2:
                results['h2g_proportion'] = float(prop_values[1])

        elif line.startswith('Proportion of SNPs:'):
            snp_values = line.replace('Proportion of SNPs:', '').strip().split()
            if len(snp_values) >= 2:
                results['snp_proportion'] = float(snp_values[1])

    # Extract coefficient and SE (for p-value calculation)
    for line in lines:
        if line.startswith('Coefficients:'):
            coeff_values = line.replace('Coefficients:', '').strip().split()
            if len(coeff_values) >= 2:
                results['enhancer_coeff'] = float(coeff_values[1])

        elif line.startswith('Coefficient SE:'):
            se_values = line.replace('Coefficient SE:', '').strip().split()
            if len(se_values) >= 2:
                results['enhancer_coeff_se'] = float(se_values[1])

    # Calculate p-value
    if 'enhancer_coeff' in results and 'enhancer_coeff_se' in results:
        z_score = results['enhancer_coeff'] / results['enhancer_coeff_se']
        from scipy.stats import norm
        results['p_value'] = 2 * (1 - norm.cdf(abs(z_score)))

    return results
```

## 📈 Results Interpretation

### Actual Analysis Results (Neg_cleaned vs Neg_unique)

| Item | Neg_cleaned | Neg_unique | Interpretation |
|------|-------------|------------|----------------|
| **Enrichment** | 53.3x | 54.1x | Both show very strong enrichment |
| **p-value** | 5.79e-06 | 2.05e-03 | Both statistically significant |
| **Heritability Contribution** | 14.86% | 4.96% | Neg_cleaned is more comprehensive |
| **SNP Ratio** | 0.28% | 0.092% | Neg_unique is more selective |

### Biological Interpretation

1. **Strong Enrichment (>50x)**:
   - PD-related variants are highly concentrated in brain enhancer regions
   - Much stronger than typical enhancer enrichment (3-10x)

2. **Neg_cleaned vs Neg_unique**:
   - **Neg_cleaned**: More comprehensive and stable signal
   - **Neg_unique**: More selective and specific signal

3. **Statistical Reliability**:
   - Both datasets significant at p < 0.01
   - Neg_cleaned provides stronger statistical evidence

## 🔧 Key Technical Solutions

### 1. Coordinate System Issue Resolution
- **Problem**: Coordinate mismatch between BED files and BaselineLD SNPs
- **Solution**: Use pre-generated annotation files (coordinate mapping completed)

### 2. File Format Compatibility
- **Problem**: LDSC annotation file format requirements
- **Solution**: Integrate BaselineLD (97) + enhancer (1) = 98 annotations

### 3. GWAS Data Compatibility
- **Problem**: SNP ID format mismatch (chr:pos vs rsID)
- **Solution**: Use GWAS file with rsID format

### 4. Statistical Validation
- **Problem**: Validate high enrichment values (>50x)
- **Solution**: Confirm relative stability through comparison with BaselineLD-only analysis

## 📝 Quality Control Checklist

### Data Validation
- [ ] All 22 chromosomes present for BaselineLD files
- [ ] Annotation file format and column count verified (98)
- [ ] LD score file generation confirmed
- [ ] GWAS file SNP ID format verified

### Analysis Validation
- [ ] LDSC regression completed successfully
- [ ] Total heritability in reasonable range (0.01-0.02)
- [ ] Enrichment values extracted successfully
- [ ] p-value calculation accuracy confirmed

### Results Validation
- [ ] Enrichment > 1 (expected direction)
- [ ] Statistical significance confirmed (p < 0.05)
- [ ] Comparable across other cell types
- [ ] Biological interpretation validity

## 🔄 Future Extensions

### Additional Analyses
1. **Remaining 6 cell type analyses**: NeuN, Nurr, Olig each cleaned/unique
2. **Cell type comparison**: Which cell type shows the strongest enrichment?
3. **Conditional analysis**: Analysis considering multiple cell types simultaneously
4. **Comparison with other diseases**: Compare with Alzheimer's, Huntington's, etc.

### Methodology Improvements
1. **More refined annotations**: Improve cell type specificity
2. **Multi-trait analysis**: Simultaneous analysis of multiple PD-related phenotypes
3. **Functional validation**: Suggest candidate regions for experimental validation

## 📚 References and Tools

### Main Tools
- **LDSC**: Bulik-Sullivan et al. (2015) Nature Genetics
- **BaselineLD v2.2**: Gazal et al. (2017) Nature Genetics
- **1000 Genomes Phase 3**: European ancestry reference panel

### Analysis Environment
- **Python**: 3.9+
- **Required Packages**: pandas, numpy, scipy, pathlib
- **System**: Linux HPC environment

---

**Document Date**: 2025-07-30
**Analysis Completed**: Neg_cleaned, Neg_unique
**Next Steps**: Proceed with remaining 6 cell type analyses
