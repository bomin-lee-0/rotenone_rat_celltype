#!/usr/bin/env python3
"""
Run single dataset LDSC with correct method
"""
import subprocess
import re
from pathlib import Path
from scipy.stats import norm

def generate_correct_ld_scores(dataset_name):
    """Generate correct LD scores"""
    print(f"🧬 Generating correct LD scores for {dataset_name}...")
    
    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'
    annot_dir = base_dir / 'ldsc_results/correct_annotations'
    output_dir = base_dir / 'ldsc_results/correct_ld_scores'
    output_dir.mkdir(exist_ok=True)
    
    completed = 0

    # Test with a few chromosomes first
    test_chromosomes = [1, 2, 3]

    for chr_num in test_chromosomes:
        output_prefix = str(output_dir / f'{dataset_name}.{chr_num}')

        # Skip if already completed
        if Path(f"{output_prefix}.l2.ldscore.gz").exists():
            print(f"  Chr {chr_num}: Already completed")
            completed += 1
            continue

        print(f"  Processing Chr {chr_num}...")

        # Input file paths
        annot_file = annot_dir / f'{dataset_name}.{chr_num}.annot.gz'
        bfile_prefix = str(base_dir / f'0_data/reference/ldsc_reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr_num}')

        if not annot_file.exists():
            print(f"    ❌ Annotation file not found")
            continue

        if not Path(f"{bfile_prefix}.bed").exists():
            print(f"    ❌ Plink file not found")
            continue

        # LDSC command
        cmd = [
            'python', 'ldsc.py', '--l2',
            '--bfile', bfile_prefix,
            '--ld-wind-cm', '1',
            '--annot', str(annot_file),
            '--out', output_prefix,
            '--print-snps', str(base_dir / '0_data/reference/ldsc_reference/listHM3.txt')
        ]
        
        try:
            result = subprocess.run(cmd, cwd=str(ldsc_dir),
                                  capture_output=True, text=True, timeout=600)

            if result.returncode == 0:
                completed += 1
                print(f"    ✅ Completed")
            else:
                print(f"    ❌ Failed: {result.stderr[:100]}")

        except subprocess.TimeoutExpired:
            print(f"    ⏰ Timeout")
        except Exception as e:
            print(f"    ❌ Exception: {e}")

    print(f"Completed chromosomes: {completed}/{len(test_chromosomes)}")
    return completed == len(test_chromosomes)

def run_correct_ldsc_regression(dataset_name):
    """Run correct LDSC regression"""
    print(f"\n📊 Running correct LDSC regression for {dataset_name}...")
    
    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'
    
    # File path settings
    gwas_file = str(base_dir / 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz')
    ref_ld_prefix = str(base_dir / f'ldsc_results/correct_ld_scores/{dataset_name}.')
    w_ld_prefix = str(base_dir / '0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.')
    frq_prefix = str(base_dir / '0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.')
    output_prefix = str(base_dir / f'ldsc_results_correct/{dataset_name}_h2')

    # Create output directory
    Path(output_prefix).parent.mkdir(exist_ok=True)

    # Check LD score files (3 chromosomes for testing)
    ld_files = list(Path(base_dir / 'ldsc_results/correct_ld_scores').glob(f'{dataset_name}.*.l2.ldscore.gz'))
    print(f"Available LD score files: {len(ld_files)}")

    if len(ld_files) < 3:
        print(f"❌ Insufficient LD score files")
        return None

    # LDSC command
    cmd = [
        'python', 'ldsc.py', '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]
    
    try:
        result = subprocess.run(cmd, cwd=str(ldsc_dir),
                              capture_output=True, text=True, timeout=900)

        if result.returncode == 0:
            print(f"✅ LDSC regression successful")
            return output_prefix + '.log'
        else:
            print(f"❌ LDSC regression failed")
            print(f"Error: {result.stderr}")
            return None

    except Exception as e:
        print(f"❌ Exception occurred: {e}")
        return None

def extract_correct_results(log_file, dataset_name):
    """Extract results for 98th annotation (our enhancer)"""
    if not log_file or not Path(log_file).exists():
        return None
        
    with open(log_file, 'r') as f:
        content = f.read()
    
    results = {'dataset': dataset_name}
    lines = content.split('\n')
    
    # Total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        results['total_h2'] = float(h2_match.group(1))
        results['total_h2_se'] = float(h2_match.group(2))
    
    # Check annotation count
    categories_line = None
    for line in lines:
        if line.startswith('Categories:'):
            categories_line = line
            categories = line.replace('Categories:', '').strip().split()
            results['n_annotations'] = len(categories)
            print(f"  Total annotation count: {len(categories)}")
            break

    # 98th (last) annotation = our enhancer
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_values = line.replace('Enrichment:', '').strip().split()
            if len(enrichment_values) >= 98:
                results['enhancer_enrichment'] = float(enrichment_values[-1])  # Last value
                print(f"  Enhancer enrichment: {enrichment_values[-1]}")
    
    for line in lines:
        if line.startswith('Proportion of h2g:'):
            prop_values = line.replace('Proportion of h2g:', '').strip().split()
            if len(prop_values) >= 98:
                results['h2g_proportion'] = float(prop_values[-1])
                
        elif line.startswith('Proportion of SNPs:'):
            snp_values = line.replace('Proportion of SNPs:', '').strip().split()
            if len(snp_values) >= 98:
                results['snp_proportion'] = float(snp_values[-1])
    
    # Coefficient and SE (last value)
    for line in lines:
        if line.startswith('Coefficients:'):
            coeff_values = line.replace('Coefficients:', '').strip().split()
            if len(coeff_values) >= 98:
                results['enhancer_coeff'] = float(coeff_values[-1])

        elif line.startswith('Coefficient SE:'):
            se_values = line.replace('Coefficient SE:', '').strip().split()
            if len(se_values) >= 98:
                results['enhancer_coeff_se'] = float(se_values[-1])

    # Calculate p-value
    if 'enhancer_coeff' in results and 'enhancer_coeff_se' in results:
        z_score = results['enhancer_coeff'] / results['enhancer_coeff_se']
        results['z_score'] = z_score
        results['p_value'] = 2 * (1 - norm.cdf(abs(z_score)))
    
    return results

def main():
    """Main execution"""
    dataset = 'Neg_cleaned'

    print("🔬 Testing correct LDSC method")
    print("="*50)

    # Step 1: Generate LD scores
    ld_success = generate_correct_ld_scores(dataset)

    if not ld_success:
        print(f"❌ LD score generation failed")
        return

    # Step 2: LDSC regression
    log_file = run_correct_ldsc_regression(dataset)

    if not log_file:
        print(f"❌ LDSC regression failed")
        return

    # Step 3: Extract results
    results = extract_correct_results(log_file, dataset)

    if results:
        print(f"\n🎯 {dataset} corrected results:")
        print("-" * 40)
        if 'n_annotations' in results:
            print(f"Number of annotations used: {results['n_annotations']}")
        if 'enhancer_enrichment' in results:
            print(f"Enrichment: {results['enhancer_enrichment']:.1f}x")
        if 'p_value' in results:
            print(f"p-value: {results['p_value']:.2e}")
        if 'h2g_proportion' in results:
            print(f"Heritability contribution: {results['h2g_proportion']*100:.2f}%")
        if 'snp_proportion' in results:
            print(f"SNP proportion: {results['snp_proportion']*100:.3f}%")

        print(f"\n✅ Success! Now ready to extend to all 22 chromosomes")
    else:
        print(f"❌ Result extraction failed")

if __name__ == '__main__':
    main()