#!/usr/bin/env python3
"""
Common Utility Classes
======================
Common data loading and processing utilities for Parkinson's disease GWAS analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
from typing import Tuple, Dict, Any


class DataManager:
    """Data loading and preprocessing management class"""

    def __init__(self, base_dir: str = "."):
        self.base_dir = Path(base_dir)
        self.data_dir = self.base_dir / "0_data"
        self.cache_dir = self.base_dir / "shared_cache"
        self.cache_dir.mkdir(exist_ok=True)

        # File paths
        self.gwas_file = self.data_dir / "raw" / "GCST009325.h.tsv.gz"
        self.enhancer_file = self.data_dir / "processed" / "Olig_cleaned_hg19_final_sorted.bed"

        # Cached data
        self._gwas_df = None
        self._enhancers_df = None
        self._classified_gwas_df = None

    def load_gwas_data(self, force_reload: bool = False) -> pd.DataFrame:
        """Load and preprocess GWAS data"""
        cache_file = self.cache_dir / "gwas_processed.pkl"

        if not force_reload and cache_file.exists() and self._gwas_df is None:
            print("Loading cached GWAS data...")
            with open(cache_file, 'rb') as f:
                self._gwas_df = pickle.load(f)
            return self._gwas_df

        if self._gwas_df is not None and not force_reload:
            return self._gwas_df

        print("Loading and preprocessing GWAS data...")

        # Load original data
        gwas_df = pd.read_csv(self.gwas_file, sep='\t', compression='gzip')
        print(f"   Original SNPs: {len(gwas_df):,}")

        # Clean data
        gwas_df = gwas_df.dropna(subset=['p_value', 'chromosome', 'base_pair_location'])
        gwas_df = gwas_df[(gwas_df['p_value'] > 0) & (gwas_df['p_value'] <= 1)]
        gwas_df = gwas_df[gwas_df['chromosome'].isin(range(1, 23))]

        # Remove duplicates (keep most significant at same position)
        gwas_df = gwas_df.sort_values('p_value').drop_duplicates(
            subset=['chromosome', 'base_pair_location'], keep='first')

        # Calculate -log10(p)
        gwas_df['neg_log10_p'] = -np.log10(gwas_df['p_value'] + 1e-300)

        print(f"   Cleaned SNPs: {len(gwas_df):,}")

        # Save cache
        with open(cache_file, 'wb') as f:
            pickle.dump(gwas_df, f)

        self._gwas_df = gwas_df
        return gwas_df

    def load_enhancer_data(self, force_reload: bool = False) -> pd.DataFrame:
        """Load and preprocess enhancer data (including coordinate conversion)"""
        cache_file = self.cache_dir / "enhancers_processed.pkl"

        if not force_reload and cache_file.exists() and self._enhancers_df is None:
            print("Loading cached enhancer data...")
            with open(cache_file, 'rb') as f:
                self._enhancers_df = pickle.load(f)
            return self._enhancers_df

        if self._enhancers_df is not None and not force_reload:
            return self._enhancers_df

        print("Loading and preprocessing enhancer data...")

        # Check for converted coordinate file first
        converted_file = self._find_converted_enhancer_file()

        if converted_file and converted_file.exists():
            print(f"   Using converted coordinate file: {converted_file}")
            enhancers_df = pd.read_csv(converted_file, sep='\t', header=None,
                                      names=['CHR', 'START', 'END', 'NAME'])
        else:
            print("   Warning: Converted coordinate file not found. Using original file (may be inaccurate)")
            print("   For accurate analysis, run setup_liftover.py first.")
            enhancers_df = pd.read_csv(self.enhancer_file, sep='\t', header=None,
                                      names=['CHR', 'START', 'END', 'NAME'])

        # Clean data
        enhancers_df['CHR'] = enhancers_df['CHR'].str.replace('chr', '')
        numeric_mask = enhancers_df['CHR'].str.isnumeric()
        enhancers_df = enhancers_df[numeric_mask].copy()
        enhancers_df['CHR'] = enhancers_df['CHR'].astype(int)
        enhancers_df = enhancers_df[enhancers_df['CHR'].isin(range(1, 23))]

        print(f"   Enhancer regions: {len(enhancers_df):,}")

        # Save cache
        with open(cache_file, 'wb') as f:
            pickle.dump(enhancers_df, f)

        self._enhancers_df = enhancers_df
        return enhancers_df

    def _find_converted_enhancer_file(self) -> Path:
        """Find converted enhancer file"""
        # Generate converted file path
        base_dir = Path(".")
        hg19_dir = base_dir / "0_data" / "processed" / "hg19_coordinates"

        if not hg19_dir.exists():
            return None

        # Estimate converted filename from original file path
        # e.g., 0_data/raw/cleaned_data/Olig_cleaned.bed -> cleaned_data_Olig_cleaned_hg19.bed
        original_path = Path(self.enhancer_file)
        parent_name = original_path.parent.name  # cleaned_data or unique_data
        file_stem = original_path.stem  # Olig_cleaned

        converted_filename = f"{parent_name}_{file_stem}_hg19.bed"
        converted_path = hg19_dir / converted_filename

        return converted_path if converted_path.exists() else None

    def classify_snps(self, force_reload: bool = False) -> pd.DataFrame:
        """Classify SNPs as inside/outside enhancer regions"""
        cache_file = self.cache_dir / "classified_gwas.pkl"

        if not force_reload and cache_file.exists() and self._classified_gwas_df is None:
            print("Loading cached classification data...")
            with open(cache_file, 'rb') as f:
                self._classified_gwas_df = pickle.load(f)
            return self._classified_gwas_df

        if self._classified_gwas_df is not None and not force_reload:
            return self._classified_gwas_df

        print("Classifying SNPs...")

        # Load data
        gwas_df = self.load_gwas_data().copy()
        enhancers_df = self.load_enhancer_data()

        # Initialize classification
        gwas_df['in_enhancer'] = False

        # Find intersection with enhancer regions
        for _, enhancer in enhancers_df.iterrows():
            mask = ((gwas_df['chromosome'] == enhancer['CHR']) &
                   (gwas_df['base_pair_location'] >= enhancer['START']) &
                   (gwas_df['base_pair_location'] <= enhancer['END']))
            gwas_df.loc[mask, 'in_enhancer'] = True

        enhancer_count = gwas_df['in_enhancer'].sum()
        total_count = len(gwas_df)

        print(f"   SNPs in enhancer: {enhancer_count:,} ({enhancer_count/total_count:.4f})")
        print(f"   Background SNPs: {total_count-enhancer_count:,}")

        # Save cache
        with open(cache_file, 'wb') as f:
            pickle.dump(gwas_df, f)

        self._classified_gwas_df = gwas_df
        return gwas_df

    def get_all_data(self, force_reload: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Return all data (GWAS, Enhancer, Classified GWAS)"""
        gwas_df = self.load_gwas_data(force_reload)
        enhancers_df = self.load_enhancer_data(force_reload)
        classified_gwas_df = self.classify_snps(force_reload)

        return gwas_df, enhancers_df, classified_gwas_df

    def clear_cache(self):
        """Delete cache files"""
        import shutil
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
            self.cache_dir.mkdir(exist_ok=True)
        print("Cache cleared")


class StatisticalAnalyzer:
    """Statistical analysis utility class"""

    @staticmethod
    def calculate_enrichment_stats(classified_gwas_df: pd.DataFrame,
                                 significance_threshold: float = 5e-8) -> Dict[str, Any]:
        """Calculate enrichment statistics (including LDSC style)"""
        from scipy.stats import mannwhitneyu, fisher_exact
        import numpy as np

        enhancer_snps = classified_gwas_df[classified_gwas_df['in_enhancer']]
        non_enhancer_snps = classified_gwas_df[~classified_gwas_df['in_enhancer']]

        # Compare p-value distributions (Mann-Whitney U test)
        enhancer_neg_log_p = enhancer_snps['neg_log10_p']
        non_enhancer_neg_log_p = non_enhancer_snps['neg_log10_p']

        mw_stat, mw_pval = mannwhitneyu(enhancer_neg_log_p, non_enhancer_neg_log_p,
                                       alternative='greater')

        # Compare significant SNP ratios (Fisher's exact test)
        enhancer_sig = (enhancer_snps['p_value'] < significance_threshold).sum()
        enhancer_total = len(enhancer_snps)
        non_enhancer_sig = (non_enhancer_snps['p_value'] < significance_threshold).sum()
        non_enhancer_total = len(non_enhancer_snps)

        contingency = [[enhancer_sig, enhancer_total - enhancer_sig],
                      [non_enhancer_sig, non_enhancer_total - non_enhancer_sig]]
        odds_ratio, fisher_pval = fisher_exact(contingency, alternative='greater')

        # Enrichment ratio (traditional method)
        enhancer_sig_rate = enhancer_sig / enhancer_total if enhancer_total > 0 else 0
        non_enhancer_sig_rate = non_enhancer_sig / non_enhancer_total if non_enhancer_total > 0 else 0
        enrichment_ratio = enhancer_sig_rate / non_enhancer_sig_rate if non_enhancer_sig_rate > 0 else np.inf

        # LDSC style enrichment calculation
        ldsc_enrichment = StatisticalAnalyzer.calculate_ldsc_enrichment(classified_gwas_df)

        return {
            'total_snps': len(classified_gwas_df),
            'enhancer_snps': enhancer_total,
            'enhancer_sig_snps': enhancer_sig,
            'background_sig_snps': non_enhancer_sig,
            'enhancer_sig_rate': enhancer_sig_rate,
            'background_sig_rate': non_enhancer_sig_rate,
            'enrichment_ratio': enrichment_ratio,
            'ldsc_enrichment': ldsc_enrichment['enrichment'],
            'ldsc_enrichment_se': ldsc_enrichment['enrichment_se'],
            'ldsc_enrichment_p': ldsc_enrichment['enrichment_p'],
            'mann_whitney_pval': mw_pval,
            'fisher_pval': fisher_pval,
            'odds_ratio': odds_ratio,
            'mann_whitney_stat': mw_stat,
            'significance_threshold': significance_threshold
        }

    @staticmethod
    def calculate_ldsc_enrichment(classified_gwas_df: pd.DataFrame) -> Dict[str, float]:
        """Calculate LDSC style enrichment"""

        enhancer_snps = classified_gwas_df[classified_gwas_df['in_enhancer']]
        all_snps = classified_gwas_df

        # Calculate effective annotation count
        M_annot = len(enhancer_snps)  # Number of SNPs in annotation
        M_total = len(all_snps)       # Total number of SNPs

        # Calculate chi-square statistics (square of z-score)
        enhancer_chisq = enhancer_snps['neg_log10_p'].apply(lambda x: (x / np.log10(np.e))**2 if x > 0 else 0)
        all_chisq = all_snps['neg_log10_p'].apply(lambda x: (x / np.log10(np.e))**2 if x > 0 else 0)

        # LDSC regression approximation calculation
        sum_chisq_annot = enhancer_chisq.sum()
        sum_chisq_total = all_chisq.sum()

        # Calculate enrichment (mean chi-square per annotation / overall mean chi-square)
        mean_chisq_annot = sum_chisq_annot / M_annot if M_annot > 0 else 0
        mean_chisq_total = sum_chisq_total / M_total if M_total > 0 else 0

        ldsc_enrichment = mean_chisq_annot / mean_chisq_total if mean_chisq_total > 0 else 0

        # Standard error approximation calculation
        var_chisq_annot = enhancer_chisq.var() if M_annot > 1 else 0
        var_chisq_total = all_chisq.var() if M_total > 1 else 0

        # SE calculation (delta method approximation)
        se_numerator = np.sqrt(var_chisq_annot / M_annot) if M_annot > 0 else 0
        se_denominator = np.sqrt(var_chisq_total / M_total) if M_total > 0 else 0

        if mean_chisq_total > 0 and M_annot > 0:
            enrichment_se = ldsc_enrichment * np.sqrt(
                (se_numerator / mean_chisq_annot)**2 + (se_denominator / mean_chisq_total)**2
            )
        else:
            enrichment_se = 0

        # Calculate Z-score and p-value
        if enrichment_se > 0:
            z_score = (ldsc_enrichment - 1) / enrichment_se
            from scipy.stats import norm
            enrichment_p = 2 * (1 - norm.cdf(abs(z_score)))  # two-tailed test
        else:
            enrichment_p = 1.0

        return {
            'enrichment': ldsc_enrichment,
            'enrichment_se': enrichment_se,
            'enrichment_p': enrichment_p,
            'mean_chisq_annot': mean_chisq_annot,
            'mean_chisq_total': mean_chisq_total,
            'n_snps_annot': M_annot,
            'n_snps_total': M_total
        }


class ManhattanPlotData:
    """Manhattan plot data preparation class"""

    @staticmethod
    def prepare_plot_data(classified_gwas_df: pd.DataFrame,
                         max_points: int = 200000) -> Dict[str, Any]:
        """Prepare data for Manhattan plot"""
        print("Preparing Manhattan plot data...")

        # Calculate cumulative position by chromosome
        plot_data = []
        cumulative_pos = 0
        chr_centers = {}
        chr_colors = {}

        # Chromosome colors (alternating)
        colors = ['#4472C4', '#70AD47']  # Blue, green alternating

        for chrom in range(1, 23):
            chr_data = classified_gwas_df[classified_gwas_df['chromosome'] == chrom].copy()
            if len(chr_data) == 0:
                continue

            # Sort by position
            chr_data = chr_data.sort_values('base_pair_location')

            # Calculate plot position
            max_pos = chr_data['base_pair_location'].max()
            chr_data['plot_pos'] = chr_data['base_pair_location'] + cumulative_pos

            # Store chromosome center position (for labels)
            chr_centers[chrom] = cumulative_pos + max_pos / 2
            chr_colors[chrom] = colors[(chrom - 1) % 2]

            cumulative_pos += max_pos + 5e6  # 5Mb gap
            plot_data.append(chr_data)

        # Combine all data
        plot_df = pd.concat(plot_data, ignore_index=True)

        # Sample if too many points (for visualization performance)
        if len(plot_df) > max_points:
            # Keep all significant SNPs, sample the rest
            significant = plot_df[plot_df['neg_log10_p'] > -np.log10(1e-4)]
            non_significant = plot_df[plot_df['neg_log10_p'] <= -np.log10(1e-4)]

            remaining_points = max_points - len(significant)
            if len(non_significant) > remaining_points and remaining_points > 0:
                sampled_non_sig = non_significant.sample(remaining_points, random_state=42)
            else:
                sampled_non_sig = non_significant

            plot_df = pd.concat([significant, sampled_non_sig], ignore_index=True)
            print(f"   Sampled for visualization: {len(plot_df):,} SNPs")

        return {
            'plot_df': plot_df,
            'chr_centers': chr_centers,
            'chr_colors': chr_colors
        }


class ResultsManager:
    """Results saving and loading management class"""

    def __init__(self, base_dir: str = "."):
        self.base_dir = Path(base_dir)
        self.results_dir = self.base_dir / "2_analysis" / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)

    def save_enrichment_results(self, results: Dict[str, Any]) -> Path:
        """Save enrichment analysis results"""
        results_df = pd.DataFrame([results])
        output_file = self.results_dir / "enrichment_results.csv"
        results_df.to_csv(output_file, index=False)

        print(f"   Results saved: {output_file}")
        return output_file

    def load_enrichment_results(self) -> Dict[str, Any]:
        """Load enrichment analysis results"""
        results_file = self.results_dir / "enrichment_results.csv"
        if not results_file.exists():
            raise FileNotFoundError("Results file not found. Run analysis first.")

        results_df = pd.read_csv(results_file)
        return results_df.iloc[0].to_dict()

    def save_detailed_snp_data(self, classified_gwas_df: pd.DataFrame):
        """Save detailed SNP data (for visualization)"""
        # Save genome-wide significant SNPs
        significant_snps = classified_gwas_df[classified_gwas_df['p_value'] < 5e-8]
        sig_file = self.results_dir / "genome_wide_significant_snps.csv"
        significant_snps.to_csv(sig_file, index=False)

        # Save significant SNPs in enhancer regions
        enhancer_sig_snps = significant_snps[significant_snps['in_enhancer']]
        enh_sig_file = self.results_dir / "enhancer_significant_snps.csv"
        enhancer_sig_snps.to_csv(enh_sig_file, index=False)

        print(f"   Detailed SNP data saved")
        print(f"      - Significant SNPs: {len(significant_snps)}")
        print(f"      - Significant SNPs in enhancer: {len(enhancer_sig_snps)}")
