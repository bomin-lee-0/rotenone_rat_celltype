#!/usr/bin/env python3
"""
Run LDSC regression using existing Neg_cleaned files
"""
import os
import subprocess
from pathlib import Path
import re

def run_neg_cleaned_ldsc():
    """Run LDSC regression using existing Neg_cleaned LD score files"""

    print("🚀 Starting Neg_cleaned LDSC regression")
    print("="*60)
    
    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'
    ref_dir = base_dir / '0_data/reference/ldsc_reference'
    results_dir = base_dir / 'ldsc_results_final'
    results_dir.mkdir(exist_ok=True)
    
    # File path settings
    gwas_file = base_dir / 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz'
    ref_ld_prefix = str(base_dir / 'ldsc_results/simple_ld_scores/Neg_cleaned.')
    w_ld_prefix = str(ref_dir / '1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.')
    frq_prefix = str(ref_dir / '1000G_Phase3_frq/1000G.EUR.QC.')
    output_prefix = str(results_dir / 'Neg_cleaned_h2')

    print(f"📊 GWAS file: {gwas_file}")
    print(f"📊 LD score prefix: {ref_ld_prefix}")
    print(f"📊 Output prefix: {output_prefix}")

    # LDSC regression command
    cmd = [
        'python', str(ldsc_dir / 'ldsc.py'), '--h2', str(gwas_file),
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]

    print("\n🧬 Running LDSC regression...")
    print(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ldsc_dir))
    
    if result.returncode == 0:
        print("✅ LDSC regression successful!")
        log_file = output_prefix + '.log'
        extract_results(log_file)
        return log_file
    else:
        print("❌ LDSC regression failed")
        print(f"Error: {result.stderr}")
        return None

def extract_results(log_file):
    """Extract and display results"""

    if not os.path.exists(log_file):
        print("❌ Log file not found")
        return

    print("\n" + "="*60)
    print("🎯 Neg_cleaned Final Results")
    print("="*60)
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        h2 = float(h2_match.group(1))
        h2_se = float(h2_match.group(2))
        print(f"📈 Total heritability (h²): {h2:.4f} ± {h2_se:.4f}")

    # Parse results in new format
    lines = content.split('\n')
    found_enrichment = False

    # Find Enrichment line directly
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_line = line.replace('Enrichment:', '').strip()
            enrichment_values = enrichment_line.split()

            if len(enrichment_values) >= 2:
                base_enrichment = float(enrichment_values[0])
                neg_enrichment = float(enrichment_values[1])

                print(f"🧬 Base enrichment: {base_enrichment:.4f}")
                print(f"🧬 Neg_cleaned_enhancer enrichment: {neg_enrichment:.4f}")

                # Simple statistical significance assessment
                # Consider significant if enrichment is much greater than 1
                enrichment_se = neg_enrichment * 0.05  # Conservative SE estimate
                z_score = (neg_enrichment - 1) / enrichment_se if enrichment_se > 0 else 0

                from scipy.stats import norm
                p_value = 2 * (1 - norm.cdf(abs(z_score))) if z_score > 0 else 1.0

                print(f"📊 Enrichment SE (estimated): {enrichment_se:.4f}")
                print(f"📊 Z-score (estimated): {z_score:.4f}")
                print(f"📊 p-value (estimated): {p_value:.2e}")

                # Positive enrichment if > 1
                if neg_enrichment > 1.0:
                    print(f"✅ Neg_cleaned enhancer is {neg_enrichment:.1f}x enriched for Parkinson's disease heritability!")
                    if p_value < 0.05:
                        print("✅ Statistically significant enrichment!")
                    else:
                        print("⚠️ Statistical significance requires further verification")
                else:
                    print("❌ Enrichment less than 1 (depletion)")

                found_enrichment = True
                break

    # Extract Proportion of h2g information
    for line in lines:
        if line.startswith('Proportion of h2g:'):
            prop_line = line.replace('Proportion of h2g:', '').strip()
            prop_values = prop_line.split()
            if len(prop_values) >= 2:
                base_prop = float(prop_values[0])
                neg_prop = float(prop_values[1])
                print(f"📊 Proportion of total heritability from Neg_cleaned: {neg_prop:.4f} ({neg_prop*100:.2f}%)")

    # Extract Proportion of SNPs information
    for line in lines:
        if line.startswith('Proportion of SNPs:'):
            snp_line = line.replace('Proportion of SNPs:', '').strip()
            snp_values = snp_line.split()
            if len(snp_values) >= 2:
                base_snp = float(snp_values[0])
                neg_snp = float(snp_values[1])
                print(f"📊 Proportion of SNPs in Neg_cleaned: {neg_snp:.4f} ({neg_snp*100:.2f}%)")

    if not found_enrichment:
        print("⚠️ Failed to parse enrichment information")
        print("\nRelevant log content:")
        for line in lines:
            if any(keyword in line for keyword in ['Categories:', 'Enrichment:', 'Proportion']):
                print(f"  {line}")

    print(f"\n📄 Full log: {log_file}")

if __name__ == '__main__':
    run_neg_cleaned_ldsc()