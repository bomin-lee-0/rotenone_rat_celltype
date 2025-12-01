#!/usr/bin/env python3
"""
LiftOver Tool Setup Script
==========================
Environment setup for rn7 → hg38 → hg19 coordinate conversion
"""

import os
import subprocess
from pathlib import Path
import urllib.request
import gzip
import shutil
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def setup_liftover_environment():
    """Set up LiftOver environment"""
    base_dir = Path(".")
    liftover_dir = base_dir / "0_data" / "reference" / "liftover_data"
    liftover_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Setting up LiftOver environment: {liftover_dir}")

    # Required files
    files_to_download = {
        'liftOver': {
            'url': 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver',
            'executable': True
        },
        'rn7ToHg38.over.chain.gz': {
            'url': 'http://hgdownload.soe.ucsc.edu/goldenPath/rn7/liftOver/rn7ToHg38.over.chain.gz',
            'executable': False
        },
        'hg38ToHg19.over.chain.gz': {
            'url': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',
            'executable': False
        }
    }

    # Download and setup
    for filename, info in files_to_download.items():
        file_path = liftover_dir / filename

        if file_path.exists():
            logger.info(f"✅ {filename} already exists")
            continue

        logger.info(f"📥 Downloading {filename}...")
        try:
            urllib.request.urlretrieve(info['url'], file_path)

            if info['executable']:
                os.chmod(file_path, 0o755)
                logger.info(f"✅ {filename} downloaded and execution permission set")
            else:
                logger.info(f"✅ {filename} download complete")

        except Exception as e:
            logger.error(f"❌ {filename} download failed: {e}")
            logger.info(f"   Manual download required: {info['url']}")

    # Verify installation
    liftover_cmd = liftover_dir / "liftOver"
    if liftover_cmd.exists():
        try:
            result = subprocess.run([str(liftover_cmd)], capture_output=True, text=True)
            logger.info("✅ liftOver tool installation verified")
        except Exception as e:
            logger.error(f"❌ liftOver execution failed: {e}")

    return liftover_dir


def convert_bed_file(input_bed: Path, output_bed: Path,
                    chain_file: Path, liftover_cmd: Path):
    """Convert BED file coordinates"""
    logger.info(f"Coordinate conversion: {input_bed} → {output_bed}")

    # File to store unmapped regions
    unmapped_file = output_bed.with_suffix('.unmapped')

    # Run liftOver
    cmd = [
        str(liftover_cmd),
        str(input_bed),
        str(chain_file),
        str(output_bed),
        str(unmapped_file)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Conversion result statistics
        if input_bed.exists() and output_bed.exists():
            original_count = sum(1 for _ in open(input_bed))
            converted_count = sum(1 for _ in open(output_bed))
            success_rate = converted_count / original_count * 100

            logger.info(f"Conversion complete: {converted_count}/{original_count} ({success_rate:.1f}%)")
            return True
        else:
            logger.error("Conversion result file was not generated.")
            return False

    except subprocess.CalledProcessError as e:
        logger.error(f"liftOver execution failed: {e}")
        logger.error(f"stderr: {e.stderr}")
        return False


def convert_all_enhancer_files():
    """Convert coordinates for all enhancer files"""
    logger.info("=" * 60)
    logger.info("🔄 Starting coordinate conversion for all enhancer files")
    logger.info("=" * 60)

    # Set up liftOver environment
    liftover_dir = setup_liftover_environment()
    liftover_cmd = liftover_dir / "liftOver"

    if not liftover_cmd.exists():
        logger.error("liftOver tool is not installed.")
        return False

    # Chain files
    rn7_to_hg38_chain = liftover_dir / "rn7ToHg38.over.chain.gz"
    hg38_to_hg19_chain = liftover_dir / "hg38ToHg19.over.chain.gz"

    if not rn7_to_hg38_chain.exists() or not hg38_to_hg19_chain.exists():
        logger.error("Required chain files are missing.")
        return False

    # Files to convert
    base_dir = Path(".")
    raw_dir = base_dir / "0_data" / "raw"

    # Create directory for converted files
    converted_dir = base_dir / "0_data" / "processed"
    converted_dir.mkdir(exist_ok=True)

    hg38_dir = converted_dir / "hg38_coordinates"
    hg19_dir = converted_dir / "hg19_coordinates"
    hg38_dir.mkdir(exist_ok=True)
    hg19_dir.mkdir(exist_ok=True)

    # List of files to convert
    bed_files = []
    for subdir in ['cleaned_data', 'unique_data']:
        subdir_path = raw_dir / subdir
        if subdir_path.exists():
            bed_files.extend(subdir_path.glob("*.bed"))

    logger.info(f"Number of files to convert: {len(bed_files)}")

    success_count = 0

    for bed_file in bed_files:
        file_stem = f"{bed_file.parent.name}_{bed_file.stem}"

        logger.info(f"\n🔄 Converting {bed_file.name}...")

        # Step 1: rn7 → hg38
        hg38_file = hg38_dir / f"{file_stem}_hg38.bed"
        logger.info("  Step 1: rn7 → hg38")

        if convert_bed_file(bed_file, hg38_file, rn7_to_hg38_chain, liftover_cmd):

            # Step 2: hg38 → hg19
            hg19_file = hg19_dir / f"{file_stem}_hg19.bed"
            logger.info("  Step 2: hg38 → hg19")

            if convert_bed_file(hg38_file, hg19_file, hg38_to_hg19_chain, liftover_cmd):
                logger.info(f"✅ {bed_file.name} conversion complete: {hg19_file}")
                success_count += 1
            else:
                logger.error(f"❌ {bed_file.name} step 2 conversion failed")
        else:
            logger.error(f"❌ {bed_file.name} step 1 conversion failed")

    logger.info(f"\n🎉 Coordinate conversion complete: {success_count}/{len(bed_files)} files")
    logger.info(f"Converted files location: {hg19_dir}")

    return success_count == len(bed_files)


if __name__ == "__main__":
    convert_all_enhancer_files()
