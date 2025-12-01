#!/usr/bin/env python3
"""
Chromosome Coordinate Conversion Utility
========================================
rn7 → hg38 → hg19 coordinate system conversion
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import subprocess
import tempfile
import os
from typing import Optional, Dict, Any, Tuple
import warnings
warnings.filterwarnings('ignore')

logger = logging.getLogger(__name__)


class CoordinateConverter:
    """Chromosome coordinate conversion class"""

    def __init__(self):
        self.base_dir = Path(".")
        self.reference_dir = self.base_dir / "0_data" / "reference"
        self.liftover_dir = self.reference_dir / "liftover_data"

        # LiftOver chain file paths
        self.chain_files = {
            'rn7_to_hg38': self.liftover_dir / "rn7ToHg38.over.chain.gz",
            'hg38_to_hg19': self.liftover_dir / "hg38ToHg19.over.chain.gz"
        }

        # LiftOver executable path
        self.liftover_cmd = self.liftover_dir / "liftOver"

        logger.info("Coordinate converter initialized")

    def check_liftover_availability(self) -> Dict[str, bool]:
        """Check LiftOver tool and chain file availability"""
        status = {
            'liftover_executable': self.liftover_cmd.exists() and os.access(self.liftover_cmd, os.X_OK),
            'rn7_to_hg38_chain': self.chain_files['rn7_to_hg38'].exists(),
            'hg38_to_hg19_chain': self.chain_files['hg38_to_hg19'].exists()
        }

        logger.info(f"LiftOver availability check:")
        for component, available in status.items():
            logger.info(f"  {component}: {'✅' if available else '❌'}")

        return status

    def detect_coordinate_system(self, bed_file: Path) -> str:
        """Estimate coordinate system of BED file"""
        logger.info(f"Estimating coordinate system: {bed_file}")

        # Read sample regions
        sample_df = pd.read_csv(bed_file, sep='\t', header=None, nrows=100,
                               names=['chr', 'start', 'end', 'name'])

        # Check coordinate range for chromosome 1
        chr1_data = sample_df[sample_df['chr'] == 'chr1']
        if len(chr1_data) == 0:
            return "unknown"

        max_pos = chr1_data['end'].max()
        min_pos = chr1_data['start'].min()

        logger.info(f"  Chromosome 1 coordinate range: {min_pos:,} - {max_pos:,}")

        # Heuristic-based coordinate system estimation
        if max_pos > 240000000:  # hg19/hg38 (~249Mb)
            if max_pos > 248000000:
                return "hg19_or_hg38"
            else:
                return "hg19_or_hg38"
        elif max_pos > 260000000:  # rn7 (~285Mb)
            return "rn7"
        else:
            return "unknown"

    def convert_bed_coordinates(self, input_bed: Path, output_bed: Path,
                               from_assembly: str, to_assembly: str) -> bool:
        """Convert BED file coordinates"""
        logger.info(f"Coordinate conversion: {from_assembly} → {to_assembly}")

        # Chain file mapping
        chain_mapping = {
            ('rn7', 'hg38'): 'rn7_to_hg38',
            ('hg38', 'hg19'): 'hg38_to_hg19'
        }

        chain_key = (from_assembly, to_assembly)
        if chain_key not in chain_mapping:
            logger.error(f"Unsupported conversion: {from_assembly} → {to_assembly}")
            return False

        chain_file = self.chain_files[chain_mapping[chain_key]]
        if not chain_file.exists():
            logger.error(f"Chain file not found: {chain_file}")
            return False

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_input:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_unmapped:
                try:
                    # Run LiftOver
                    cmd = [
                        str(self.liftover_cmd),
                        str(input_bed),
                        str(chain_file),
                        str(output_bed),
                        tmp_unmapped.name
                    ]

                    result = subprocess.run(cmd, capture_output=True, text=True)

                    if result.returncode == 0:
                        logger.info(f"Coordinate conversion successful: {output_bed}")

                        # Conversion statistics
                        if output_bed.exists():
                            original_count = sum(1 for _ in open(input_bed))
                            converted_count = sum(1 for _ in open(output_bed))
                            success_rate = converted_count / original_count * 100

                            logger.info(f"  Conversion success rate: {converted_count}/{original_count} ({success_rate:.1f}%)")

                        return True
                    else:
                        logger.error(f"LiftOver execution failed: {result.stderr}")
                        return False

                finally:
                    # Clean up temporary files
                    try:
                        os.unlink(tmp_input.name)
                        os.unlink(tmp_unmapped.name)
                    except:
                        pass

    def setup_liftover_environment(self):
        """Set up LiftOver environment"""
        logger.info("Setting up LiftOver environment...")

        # Create directory
        self.liftover_dir.mkdir(parents=True, exist_ok=True)

        # Download URLs for required files
        download_urls = {
            'liftOver': 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver',
            'rn7ToHg38.over.chain.gz': 'http://hgdownload.soe.ucsc.edu/goldenPath/rn7/liftOver/rn7ToHg38.over.chain.gz',
            'hg38ToHg19.over.chain.gz': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
        }

        logger.info("Automatic download is not supported.")
        logger.info("Please download the following files manually:")
        for filename, url in download_urls.items():
            target_path = self.liftover_dir / filename
            logger.info(f"  {filename}: {url}")
            logger.info(f"    → {target_path}")

        return False


class EnhancedDataManager:
    """Data manager with coordinate conversion support"""

    def __init__(self, enhancer_file: Path):
        self.enhancer_file = enhancer_file
        self.converter = CoordinateConverter()
        self.cache_dir = Path("coordinate_conversion_cache")
        self.cache_dir.mkdir(exist_ok=True)

        logger.info(f"Enhanced data manager initialized: {enhancer_file}")

    def get_converted_enhancer_data(self, target_assembly: str = "hg19") -> Optional[pd.DataFrame]:
        """Return coordinate-converted enhancer data"""

        # Check cache file
        cache_file = self.cache_dir / f"{self.enhancer_file.stem}_{target_assembly}.pkl"
        if cache_file.exists():
            logger.info(f"Loading cached converted data: {cache_file}")
            return pd.read_pickle(cache_file)

        # Estimate original coordinate system
        original_assembly = self.converter.detect_coordinate_system(self.enhancer_file)
        logger.info(f"Estimated original coordinate system: {original_assembly}")

        if original_assembly == target_assembly:
            logger.info("Coordinate conversion not required.")
            enhancer_df = pd.read_csv(self.enhancer_file, sep='\t', header=None,
                                    names=['CHR', 'START', 'END', 'NAME'])
            # Clean data
            enhancer_df['CHR'] = enhancer_df['CHR'].str.replace('chr', '')
            numeric_mask = enhancer_df['CHR'].str.isnumeric()
            enhancer_df = enhancer_df[numeric_mask].copy()
            enhancer_df['CHR'] = enhancer_df['CHR'].astype(int)
            enhancer_df = enhancer_df[enhancer_df['CHR'].isin(range(1, 23))]

            # Save cache
            enhancer_df.to_pickle(cache_file)
            return enhancer_df

        # Coordinate conversion required
        logger.warning("Coordinate conversion required but LiftOver tool is not configured.")
        logger.info("Using original data as-is for now.")
        logger.info("For accurate analysis, please perform coordinate conversion.")

        # Temporarily return original data
        enhancer_df = pd.read_csv(self.enhancer_file, sep='\t', header=None,
                                names=['CHR', 'START', 'END', 'NAME'])
        # Clean data
        enhancer_df['CHR'] = enhancer_df['CHR'].str.replace('chr', '')
        numeric_mask = enhancer_df['CHR'].str.isnumeric()
        enhancer_df = enhancer_df[numeric_mask].copy()
        enhancer_df['CHR'] = enhancer_df['CHR'].astype(int)
        enhancer_df = enhancer_df[enhancer_df['CHR'].isin(range(1, 23))]

        return enhancer_df


def main():
    """Test and validation"""
    converter = CoordinateConverter()

    # Environment check
    status = converter.check_liftover_availability()

    if not all(status.values()):
        logger.warning("LiftOver environment is not complete.")
        converter.setup_liftover_environment()

    # Estimate coordinate system for sample files
    sample_files = [
        Path("0_data/raw/cleaned_data/Olig_cleaned.bed"),
        Path("0_data/raw/unique_data/Olig_unique.bed")
    ]

    for sample_file in sample_files:
        if sample_file.exists():
            assembly = converter.detect_coordinate_system(sample_file)
            logger.info(f"{sample_file}: {assembly}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
