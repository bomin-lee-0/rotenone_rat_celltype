#!/usr/bin/env python3
"""
Manhattan plot for cell type-specific enhancer intersect SNPs
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import gzip
from pathlib import Path
import os

# Set project root directory
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_ROOT / "0.Data"
RESULTS_DIR = PROJECT_ROOT / "2.Results"

def load_celltype_annotations():
    """Load cell type-specific annotation data"""

    print("🧬 Loading cell type-specific annotation data...")

    celltypes = {
        'Microglia': 'Neg_cleaned',
        'Neuron': 'NeuN_cleaned',
        'Oligodendrocyte': 'Olig_cleaned',
        'Dopaminergic': 'Nurr_cleaned'
    }

    all_annotations = {}

    for celltype, file_prefix in celltypes.items():
        print(f"  Loading {celltype}...")

        # First test with chromosome 1 only
        # If LDSC annotation file doesn't exist, load region info from BED file
        annot_file = DATA_DIR / "Results" / "annotations" / f"{file_prefix}.1.annot.gz"
        bed_file = DATA_DIR / "processed" / "hg19_coordinates" / f"cleaned_data_{file_prefix}_hg19.bed"

        try:
            if annot_file.exists():
                # Use LDSC annotation file
                df = pd.read_csv(annot_file, sep='\t', compression='gzip')

                # Last column is enhancer annotation
                enhancer_col = df.columns[-1]  # Last column

                # Select only intersect SNPs
                intersect_snps = df[df[enhancer_col] == 1]['SNP'].tolist()

                print(f"    Chromosome 1: {len(intersect_snps)} SNPs intersect with {celltype} enhancer")

                all_annotations[celltype] = {
                    'intersect_snps': set(intersect_snps),
                    'total_snps': len(df),
                    'intersect_count': len(intersect_snps)
                }
            elif bed_file.exists():
                # Load enhancer regions from BED file
                bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name'])

                print(f"    Loaded {len(bed_df)} enhancer regions from BED file")
                print(f"    Note: SNP-enhancer overlap will be calculated after GWAS data loading")

                all_annotations[celltype] = {
                    'bed_regions': bed_df,
                    'intersect_snps': set(),  # Calculate later
                    'total_snps': 0,
                    'intersect_count': 0
                }
            else:
                print(f"    Warning: Cannot find annotation file and BED file")
                all_annotations[celltype] = {
                    'intersect_snps': set(),
                    'total_snps': 0,
                    'intersect_count': 0
                }

        except Exception as e:
            print(f"    Error: {e}")
            all_annotations[celltype] = {
                'intersect_snps': set(),
                'total_snps': 0,
                'intersect_count': 0
            }

    return all_annotations

def compute_snp_overlaps(gwas_df, annotations):
    """Calculate overlap between BED regions and GWAS SNPs"""

    print("\n🔍 Calculating SNP-enhancer overlap...")

    for celltype, data in annotations.items():
        if 'bed_regions' in data and data['bed_regions'] is not None:
            print(f"  Calculating {celltype} overlap...")

            bed_df = data['bed_regions']
            intersect_snps = []

            # Check enhancer region overlap for each SNP
            for _, snp in gwas_df.iterrows():
                snp_chr = snp['CHR']
                snp_pos = snp['BP']

                # Check if SNP is contained in any enhancer region
                overlaps = bed_df[
                    (bed_df['chr'].astype(str) == str(snp_chr)) &
                    (bed_df['start'] <= snp_pos) &
                    (bed_df['end'] >= snp_pos)
                ]

                if len(overlaps) > 0:
                    intersect_snps.append(snp['SNP'])

            data['intersect_snps'] = set(intersect_snps)
            data['intersect_count'] = len(intersect_snps)

            print(f"    {len(intersect_snps)} SNPs overlap with {celltype} enhancer")

    return annotations

def load_gwas_with_positions():
    """Load GWAS data with position information"""

    print("📊 Loading GWAS data...")

    # GWAS data
    gwas_file = DATA_DIR / "GWAS" / "GCST009325.h.tsv.gz"

    # Sample loading (memory saving)
    sample_size = 100000  # Test: reduced to 100k

    print(f"  Loading GWAS data with sampling (n={sample_size:,})...")

    df = pd.read_csv(gwas_file, sep='\t', compression='gzip', nrows=sample_size)

    # Check and standardize column names
    print(f"  Columns: {', '.join(df.columns[:10])}...")

    # Find P-value column
    p_col = None
    for col in ['p_value', 'P', 'pvalue', 'PVAL']:
        if col in df.columns:
            p_col = col
            break

    if p_col:
        df['P'] = df[p_col]
    elif 'Z' in df.columns:
        # Calculate P-value from Z-score
        df['P'] = 2 * (1 - norm.cdf(np.abs(df['Z'])))
    else:
        print("  Warning: Cannot find P-value column!")
        df['P'] = 0.5

    df['-log10P'] = -np.log10(np.maximum(df['P'], 1e-50))

    # Add position information
    chr_col = next((c for c in ['chromosome', 'CHR', 'chr'] if c in df.columns), None)
    bp_col = next((c for c in ['base_pair_location', 'BP', 'bp', 'POS'] if c in df.columns), None)
    snp_col = next((c for c in ['variant_id', 'rsid', 'SNP', 'ID'] if c in df.columns), None)

    if chr_col:
        df['CHR'] = df[chr_col]
    else:
        df['CHR'] = 1

    if bp_col:
        df['BP'] = df[bp_col]
    else:
        df['BP'] = range(len(df))

    if snp_col:
        df['SNP'] = df[snp_col]
    elif 'rsid' in df.columns:
        df['SNP'] = df['rsid']
    else:
        df['SNP'] = [f"SNP_{i}" for i in range(len(df))]

    print(f"  Loading complete: {len(df):,} SNPs")
    print(f"  P-value range: {df['P'].min():.2e} - {df['P'].max():.2e}")

    return df

def create_celltype_manhattan_plots(gwas_df, annotations):
    """Create cell type-specific Manhattan plots"""

    print("🗽 Creating cell type-specific Manhattan plots...")

    # 4 cell type subplots
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    axes = axes.flatten()

    celltypes = ['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic']
    colors_intersect = ['red', 'darkgreen', 'purple', 'orange']
    colors_background = ['lightcoral', 'lightgreen', 'plum', 'moccasin']

    # Significance thresholds
    genome_wide = -np.log10(5e-8)
    suggestive = -np.log10(1e-5)

    for i, celltype in enumerate(celltypes):
        ax = axes[i]

        print(f"  Creating {celltype} plot...")

        # Intersect SNPs for this cell type
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_count = len(intersect_snps)

        # Mark intersect status in GWAS data
        gwas_df['is_intersect'] = gwas_df['SNP'].isin(intersect_snps)

        # Non-intersect SNPs (background)
        non_intersect = gwas_df[~gwas_df['is_intersect']]
        intersect_data = gwas_df[gwas_df['is_intersect']]

        print(f"    Non-intersect SNPs: {len(non_intersect):,}")
        print(f"    Intersect SNPs: {len(intersect_data):,}")

        # Background SNPs (small, transparent)
        if len(non_intersect) > 0:
            # Sample if too many
            if len(non_intersect) > 50000:
                non_intersect_sample = non_intersect.sample(n=50000, random_state=42)
            else:
                non_intersect_sample = non_intersect

            ax.scatter(non_intersect_sample['BP'], non_intersect_sample['-log10P'],
                      c=colors_background[i], alpha=0.3, s=1, label='Background SNPs')

        # Intersect SNPs (large, opaque)
        if len(intersect_data) > 0:
            ax.scatter(intersect_data['BP'], intersect_data['-log10P'],
                      c=colors_intersect[i], alpha=0.8, s=20,
                      label=f'{celltype} Enhancer SNPs (n={len(intersect_data)})')

        # Significance thresholds
        ax.axhline(y=genome_wide, color='red', linestyle='--', alpha=0.7,
                  linewidth=1.5, label='p=5×10⁻⁸')
        ax.axhline(y=suggestive, color='blue', linestyle='--', alpha=0.5,
                  linewidth=1, label='p=1×10⁻⁵')

        # Label top intersect SNPs
        if len(intersect_data) > 0:
            top_intersect = intersect_data.nlargest(3, '-log10P')
            for _, snp in top_intersect.iterrows():
                if snp['-log10P'] > suggestive:
                    ax.annotate(f"{snp['SNP']}",
                               xy=(snp['BP'], snp['-log10P']),
                               xytext=(5, 5), textcoords='offset points',
                               fontsize=8, alpha=0.8,
                               bbox=dict(boxstyle='round,pad=0.2',
                                       facecolor=colors_intersect[i], alpha=0.7))

        # Axis settings
        ax.set_xlabel('SNP Position (Chr 1)', fontsize=12, fontweight='bold')
        ax.set_ylabel('-log₁₀(P-value)', fontsize=12, fontweight='bold')
        ax.set_title(f'{celltype} Enhancer Manhattan Plot\n'
                    f'({intersect_count} intersect SNPs)',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=10)

        # Set Y-axis range
        max_y = max(gwas_df['-log10P'].max(), genome_wide + 2)
        ax.set_ylim(0, max_y)

    plt.tight_layout()

    # Save
    plots_dir = RESULTS_DIR / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(plots_dir / 'celltype_manhattan_plots.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(plots_dir / 'celltype_manhattan_plots.pdf',
                bbox_inches='tight', facecolor='white')

    print(f"  Saved: {plots_dir / 'celltype_manhattan_plots.png'}")

    plt.show()

    # Statistics summary
    print(f"\n📈 Cell type-specific Manhattan plot statistics:")
    print("="*60)

    for celltype in celltypes:
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_data = gwas_df[gwas_df['SNP'].isin(intersect_snps)]

        if len(intersect_data) > 0:
            significant_intersect = intersect_data[intersect_data['-log10P'] > suggestive]
            genome_wide_intersect = intersect_data[intersect_data['-log10P'] > genome_wide]

            print(f"\n🧬 {celltype}:")
            print(f"  Total intersect SNPs: {len(intersect_data)}")
            print(f"  Suggestive (p<1e-5): {len(significant_intersect)}")
            print(f"  Genome-wide (p<5e-8): {len(genome_wide_intersect)}")

            if len(intersect_data) > 0:
                print(f"  Max -log10(P): {intersect_data['-log10P'].max():.2f}")

                # Top 3 SNPs
                top_snps = intersect_data.nlargest(3, '-log10P')
                print(f"  Top SNPs:")
                for _, snp in top_snps.iterrows():
                    print(f"    {snp['SNP']}: p={snp['P']:.2e}")
        else:
            print(f"\n🧬 {celltype}: No intersect SNPs found")

def create_comparison_manhattan(gwas_df, annotations):
    """4 cell type comparison Manhattan plot"""

    print("\n🔄 Creating cell type comparison Manhattan plot...")

    # Compare 4 cell types in a single plot
    fig, ax = plt.subplots(figsize=(16, 10))

    celltypes = ['Microglia', 'Neuron', 'Oligodendrocyte', 'Dopaminergic']
    colors = ['red', 'green', 'blue', 'orange']

    # Plot points for each cell type
    for i, celltype in enumerate(celltypes):
        intersect_snps = annotations[celltype]['intersect_snps']
        intersect_data = gwas_df[gwas_df['SNP'].isin(intersect_snps)]

        if len(intersect_data) > 0:
            ax.scatter(intersect_data['BP'], intersect_data['-log10P'],
                      c=colors[i], alpha=0.7, s=15, label=f'{celltype} (n={len(intersect_data)})')

    # Significance thresholds
    genome_wide = -np.log10(5e-8)
    suggestive = -np.log10(1e-5)

    ax.axhline(y=genome_wide, color='red', linestyle='--', alpha=0.8,
              linewidth=2, label='Genome-wide (p=5×10⁻⁸)')
    ax.axhline(y=suggestive, color='blue', linestyle='--', alpha=0.6,
              linewidth=1.5, label='Suggestive (p=1×10⁻⁵)')

    # Axis settings
    ax.set_xlabel('SNP Position (Chr 1)', fontsize=14, fontweight='bold')
    ax.set_ylabel('-log₁₀(P-value)', fontsize=14, fontweight='bold')
    ax.set_title('Cell Type-Specific Enhancer Manhattan Plot Comparison',
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Save
    plots_dir = RESULTS_DIR / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(plots_dir / 'celltype_comparison_manhattan.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(plots_dir / 'celltype_comparison_manhattan.pdf',
                bbox_inches='tight', facecolor='white')

    print(f"  Saved: {plots_dir / 'celltype_comparison_manhattan.png'}")

    plt.show()

if __name__ == "__main__":
    print("🗽 Creating cell type-specific enhancer Manhattan plots!")
    print("="*60)

    # 1. Load cell type-specific annotations
    annotations = load_celltype_annotations()

    # 2. Load GWAS data
    gwas_data = load_gwas_with_positions()

    # 3. Calculate SNP-enhancer overlap when using BED files
    annotations = compute_snp_overlaps(gwas_data, annotations)

    # 4. Create cell type-specific Manhattan plots
    create_celltype_manhattan_plots(gwas_data, annotations)

    # 5. Create comparison Manhattan plot
    create_comparison_manhattan(gwas_data, annotations)

    print(f"\n✅ Cell type-specific Manhattan plots created!")
    print(f"   Saved files:")
    plots_dir = RESULTS_DIR / "Plots"
    print(f"   - {plots_dir / 'celltype_manhattan_plots.png'} (4-panel)")
    print(f"   - {plots_dir / 'celltype_comparison_manhattan.png'} (overlay)")
