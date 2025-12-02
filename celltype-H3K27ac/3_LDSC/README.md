# LDSC Partitioned Heritability Analysis for Parkinson's Disease

## Overview
Cell type-specific genetic enrichment analysis of Parkinson's disease GWAS signals
across four brain cell types using LDSC (Linkage Disequilibrium Score Regression).

### Data Scale
- **GWAS:** 17.4M SNPs (37,688 cases / 1.4M controls)
- **Source:** GCST009325 - Nalls et al. (2019)

## Cell Types Analyzed
| Cell Type | Marker | Biological Role |
|-----------|--------|-----------------|
| Olig | Oligodendrocytes | Myelin formation |
| Nurr | Dopaminergic Neurons | Dopamine production |
| NeuN | General Neurons | Neural signaling |
| Neg | Microglia | Immune response |

## Folder Structure
```
├── 0.Data/          # GWAS, enhancer BED files, references
├── 1.Scripts/       # LDSC, visualization, utility scripts
├── 2.Results/       # Output files and plots
└── 3.Documentation/ # Workflow guides
```

## Quick Start
```bash
# Run complete analysis
python main.py --all

# Or step by step
python main.py --step coordinate  # rn7 → hg38 → hg19 liftover
python main.py --step ldsc        # LDSC partitioned heritability
python main.py --step visualize   # Manhattan plots
```

## Estimated Runtime
- First run: ~1-2 hours
- Subsequent runs: ~10 minutes (cached)

## Output
- Enrichment values per cell type
- P-values and standard errors
- Manhattan plots
- Summary report (CSV/MD)

## References
- Nalls et al. (2019) Lancet Neurology - PD GWAS
- Finucane et al. (2015) Nature Genetics - LDSC method
