#!/usr/bin/env python3
"""
Final Neg_cleaned LDSC regression execution and results analysis
"""
import subprocess
import os
import re

def fix_all_m_files():
    """Fix M files for all chromosomes correctly"""

    print('🔧 Fixing all chromosome M files...')

    for chr_num in range(1, 23):
        baseline_m_file = f'/scratch/prj/eng_waste_to_protein/repositories/bomin/0_data/reference/ldsc_reference/baselineLD.{chr_num}.l2.M'
        celltype_m_file = f'/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/simple_ld_scores/Neg_cleaned.{chr_num}.l2.M'
        output_m_file = f'/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/test_combined_ld_scores/Neg_cleaned.{chr_num}.l2.M'

        if not os.path.exists(baseline_m_file) or not os.path.exists(celltype_m_file):
            print(f'  ⚠️ Chr{chr_num}: File not found')
            continue

        # Read BaselineLD M (97 values)
        with open(baseline_m_file, 'r') as f:
            baseline_values = f.read().strip().split('\t')

        # Extract enhancer value from Neg_cleaned M
        with open(celltype_m_file, 'r') as f:
            lines = f.read().strip().split('\n')
            if len(lines) >= 2:
                enhancer_value = lines[1]  # Second line
            else:
                enhancer_value = lines[0].split('\t')[1] if '\t' in lines[0] else lines[0]

        # Combine: 97 + 1 = 98
        combined_values = baseline_values + [enhancer_value]

        # Save as single line
        with open(output_m_file, 'w') as f:
            f.write('\t'.join(combined_values) + '\n')

        print(f'{chr_num}', end=' ', flush=True)

    print('\n✅ All M files fixed')

def run_ldsc_regression():
    """Run LDSC regression"""

    print('🧬 Running LDSC regression...')

    ldsc_path = '/scratch/prj/eng_waste_to_protein/repositories/bomin/1_preprocessing/ldsc-python3/ldsc.py'
    gwas_file = '/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/sumstats/parkinson_gwas_fixed.sumstats.gz'
    ref_ld_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/test_combined_ld_scores/Neg_cleaned.'
    w_ld_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
    frq_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.'
    output_prefix = '/scratch/prj/eng_waste_to_protein/repositories/bomin/ldsc_results/test_results/Neg_cleaned_FINAL'

    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    cmd = [
        'python', ldsc_path, '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--overlap-annot',
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]

    result = subprocess.run(cmd, capture_output=True, text=True,
                          cwd='/scratch/prj/eng_waste_to_protein/repositories/bomin/1_preprocessing/ldsc-python3')

    if result.returncode == 0:
        print('✅ LDSC regression successful!')
        return output_prefix + '.log'
    else:
        print('❌ LDSC regression failed')
        print(f'Error: {result.stderr[:1000]}')
        return None

def extract_results(log_file):
    """Extract enrichment and p-value from log file"""

    if not log_file or not os.path.exists(log_file):
        print('❌ Log file not found.')
        return

    print('📊 Analyzing results...')

    with open(log_file, 'r') as f:
        log_content = f.read()

    # Extract results
    print('\n' + '='*60)
    print("🧬 Parkinson's Disease GWAS - Neg_cleaned Enhancer Analysis Results")
    print('='*60)

    # Total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', log_content)
    if h2_match:
        h2_value = h2_match.group(1)
        h2_se = h2_match.group(2)
        print(f'📈 Total heritability (h²): {h2_value} ± {h2_se}')

    # Find Neg_cleaned enhancer enrichment
    lines = log_content.split('\n')
    for i, line in enumerate(lines):
        if 'Neg_cleaned_enhancer' in line and ('Enrichment:' in line or 'Coefficient:' in line):
            print(f'🎯 {line.strip()}')
        elif 'Categories:' in line:
            # Print results by category
            print(f'\n📋 Results by category:')
            for j in range(1, min(10, len(lines) - i)):
                next_line = lines[i + j].strip()
                if next_line and not next_line.startswith('Analysis') and 'L2' in next_line:
                    if 'Neg_cleaned_enhancer' in next_line:
                        print(f'🎯 **{next_line}**')
                    elif j <= 5:  # Print first few only
                        print(f'   {next_line}')

    print(f'\n📄 Full log: {log_file}')
    return True

def main():
    """Main execution"""

    print('🚀 Starting final Neg_cleaned LDSC analysis')

    # 1. Fix M files
    fix_all_m_files()

    # 2. Run LDSC regression
    log_file = run_ldsc_regression()

    # 3. Extract results
    if log_file:
        extract_results(log_file)
    else:
        print('❌ Analysis failed')

if __name__ == '__main__':
    main()
