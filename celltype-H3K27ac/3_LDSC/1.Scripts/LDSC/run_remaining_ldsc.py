#!/usr/bin/env python3
"""
Run LDSC regression for remaining 5 datasets
"""
import subprocess
import re
from pathlib import Path
from scipy.stats import norm

def run_ldsc_for_dataset(dataset_name):
    """Run LDSC regression for a specific dataset"""
    print(f"\n🧬 Starting {dataset_name} LDSC regression...")

    base_dir = Path('/scratch/prj/eng_waste_to_protein/repositories/bomin')
    ldsc_dir = base_dir / '1_preprocessing/ldsc-python3'

    # Set file paths
    gwas_file = str(base_dir / 'ldsc_results/sumstats/parkinson_gwas.sumstats.gz')
    ref_ld_prefix = str(base_dir / f'ldsc_results/simple_ld_scores/{dataset_name}.')
    w_ld_prefix = str(base_dir / '0_data/reference/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.')
    frq_prefix = str(base_dir / '0_data/reference/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC.')
    output_prefix = str(base_dir / f'ldsc_results_final/{dataset_name}_h2')

    # Check input file existence
    if not Path(gwas_file).exists():
        print(f"❌ GWAS file not found: {gwas_file}")
        return None

    # Check LD score file count (22 chromosomes)
    ld_files = list(Path(base_dir / 'ldsc_results/simple_ld_scores').glob(f'{dataset_name}.*.l2.ldscore.gz'))
    if len(ld_files) != 22:
        print(f"❌ Insufficient LD score files: {len(ld_files)}/22")
        return None

    # Construct LDSC command
    cmd = [
        'python', 'ldsc.py', '--h2', gwas_file,
        '--ref-ld-chr', ref_ld_prefix,
        '--w-ld-chr', w_ld_prefix,
        '--frqfile-chr', frq_prefix,
        '--out', output_prefix
    ]

    try:
        # Run LDSC
        result = subprocess.run(cmd, cwd=str(ldsc_dir),
                              capture_output=True, text=True, timeout=1800)

        if result.returncode == 0:
            print(f"✅ {dataset_name}: LDSC regression successful")
            return output_prefix + '.log'
        else:
            print(f"❌ {dataset_name}: LDSC regression failed")
            print(f"Error: {result.stderr}")
            return None

    except subprocess.TimeoutExpired:
        print(f"⏰ {dataset_name}: Timeout (30 minutes)")
        return None
    except Exception as e:
        print(f"❌ {dataset_name}: Exception occurred - {e}")
        return None

def extract_results(log_file, dataset_name):
    """Extract LDSC results"""
    if not log_file or not Path(log_file).exists():
        return None

    with open(log_file, 'r') as f:
        content = f.read()

    results = {'dataset': dataset_name}

    # Total heritability
    h2_match = re.search(r'Total Observed scale h2: ([\d\.]+) \(([\d\.]+)\)', content)
    if h2_match:
        results['total_h2'] = float(h2_match.group(1))
        results['total_h2_se'] = float(h2_match.group(2))

    lines = content.split('\n')

    # Extract enrichment values
    for line in lines:
        if line.strip().startswith('Enrichment:'):
            enrichment_values = line.replace('Enrichment:', '').strip().split()
            if len(enrichment_values) >= 2:
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

    # Extract coefficient and SE
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
        results['z_score'] = z_score
        results['p_value'] = 2 * (1 - norm.cdf(abs(z_score)))

    return results

def main():
    """Main execution function"""
    print("🔬 LDSC analysis for remaining 5 datasets")
    print("="*60)

    datasets = ['NeuN_cleaned', 'NeuN_unique', 'Nurr_cleaned', 'Nurr_unique', 'Olig_cleaned']
    all_results = []

    for dataset in datasets:
        log_file = run_ldsc_for_dataset(dataset)
        if log_file:
            results = extract_results(log_file, dataset)
            if results:
                all_results.append(results)

                # Print results
                print(f"\n📊 {dataset} Results:")
                if 'enhancer_enrichment' in results:
                    print(f"  Enrichment: {results['enhancer_enrichment']:.1f}x")
                if 'p_value' in results:
                    print(f"  p-value: {results['p_value']:.2e}")
                if 'h2g_proportion' in results:
                    print(f"  Heritability contribution: {results['h2g_proportion']*100:.2f}%")
                if 'snp_proportion' in results:
                    print(f"  SNP proportion: {results['snp_proportion']*100:.3f}%")

    # Overall results summary
    print(f"\n🎯 Overall Results Summary ({len(all_results)}/5 completed)")
    print("-" * 80)
    print("Dataset\t\tEnrichment\tp-value\t\th2g%\tSNP%")
    print("-" * 80)

    for r in all_results:
        enrich = r.get('enhancer_enrichment', 0)
        pval = r.get('p_value', 1)
        h2g = r.get('h2g_proportion', 0) * 100
        snp = r.get('snp_proportion', 0) * 100

        print(f"{r['dataset']:<15}\t{enrich:.1f}x\t{pval:.2e}\t{h2g:.2f}%\t{snp:.3f}%")

    return all_results

if __name__ == '__main__':
    results = main()
    print(f"\n✅ Completed: {len(results)} datasets")
