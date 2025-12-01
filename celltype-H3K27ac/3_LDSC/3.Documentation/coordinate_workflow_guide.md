# Coordinate Conversion Workflow Guide

## ⚠️ Important: Coordinate System Mismatch Issue

The reason for `Enrichment ratio: 0.000` in the current analysis is due to **coordinate system mismatch**:

- **GWAS data (GCST009325)**: hg19 (GRCh37) coordinate system
- **Enhancer BED files**: Based on rn7 coordinate system

## 🔄 Solution: Coordinate Conversion

### Step 1: Download LiftOver Tool and Chain Files

```bash
# Set up coordinate conversion environment (automatic download attempt)
python setup_liftover.py
```

If automatic download fails, download manually:

#### Files Required for Manual Download:

1. **liftOver executable**
   ```bash
   cd 0_data/reference/liftover_data
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
   chmod +x liftOver
   ```

2. **rn7 → hg38 chain file**
   ```bash
   wget http://hgdownload.soe.ucsc.edu/goldenPath/rn7/liftOver/rn7ToHg38.over.chain.gz
   ```

3. **hg38 → hg19 chain file**
   ```bash
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
   ```

### Step 2: Run Coordinate Conversion

```bash
# Convert all BED files from rn7 → hg38 → hg19
python setup_liftover.py
```

This command performs the following conversions:
- `0_data/raw/cleaned_data/*.bed` → `0_data/processed/hg19_coordinates/*_hg19.bed`
- `0_data/raw/unique_data/*.bed` → `0_data/processed/hg19_coordinates/*_hg19.bed`

### Step 3: Verify Conversion

After conversion is complete, the following directory structure will be created:

```
0_data/
├── raw/
│   ├── cleaned_data/*.bed         # Original (rn7)
│   └── unique_data/*.bed          # Original (rn7)
└── processed/
    ├── hg38_coordinates/          # Intermediate conversion (hg38)
    └── hg19_coordinates/          # Final conversion (hg19) ⭐
        ├── cleaned_data_Olig_cleaned_hg19.bed
        ├── cleaned_data_Nurr_cleaned_hg19.bed
        ├── unique_data_Olig_unique_hg19.bed
        └── ...
```

### Step 4: Re-run Batch Analysis

After coordinate conversion is complete:

```bash
# Clear cache and re-analyze with converted coordinates
python run_complete_batch_pipeline.py --force-reload
```

## 🔍 Expected Results After Conversion

Expected changes after coordinate conversion:

- **Before**: `Enrichment ratio: 0.000` (coordinate mismatch)
- **After**: `Enrichment ratio: 1.2-3.5` (normal enrichment)

## ⚡ Quick Execution Guide

```bash
# 1. Coordinate conversion (only once)
python setup_liftover.py

# 2. Re-run batch analysis
python run_complete_batch_pipeline.py --force-reload
```

## 📊 Conversion Quality Check

After conversion, verify the following:

1. **Conversion rate**: Should typically be 90% or higher
2. **Number of regions**: Should be similar to original
3. **Enrichment ratio**: Reasonable value, not zero

## 🚫 Troubleshooting

### Problem 1: liftOver Download Failed
- Manually download and save to `0_data/reference/liftover_data/`
- Check execution permission: `chmod +x liftOver`

### Problem 2: Chain File Missing
- Download directly from UCSC site
- Check file integrity

### Problem 3: Low Conversion Rate (<80%)
- Check original BED file format
- Check chromosome naming convention (chr1 vs 1)

---

**Important**: This conversion process is crucial for the accuracy of the entire analysis.
Accurate enrichment analysis is not possible without conversion.
