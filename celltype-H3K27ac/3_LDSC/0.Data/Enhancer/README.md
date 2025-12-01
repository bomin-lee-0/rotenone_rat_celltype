# Enhancer Data Directory

This folder stores brain cell type-specific enhancer BED files.

## File Description

### Cell Type-Specific Enhancer Files

- **Olig_cleaned.bed / Olig_unique.bed**: Oligodendrocyte enhancer
- **Nurr_cleaned.bed / Nurr_unique.bed**: Dopaminergic neuron enhancer
- **NeuN_cleaned.bed / NeuN_unique.bed**: General neuron enhancer
- **Neg_cleaned.bed / Neg_unique.bed**: Microglia enhancer

### Processing Methods

- **cleaned**: Preprocessed enhancer regions
- **unique**: Extracted unique enhancer regions only

## Coordinate System

- Original: rn7 (rat genome)
- Conversion required: rn7 → hg38 → hg19 (human genome)
- Conversion tool: `1.Scripts/Utils/setup_liftover.py`
