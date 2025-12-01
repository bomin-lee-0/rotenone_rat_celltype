#!/usr/bin/env python3
"""
Parkinson's Disease GWAS - Cell Type-Specific Enhancer Enrichment Analysis
Main execution script
"""

import os
import sys
import argparse
from pathlib import Path

# Add scripts directories to path
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "LDSC"))
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "Visualization"))
sys.path.append(str(Path(__file__).parent / "1.Scripts" / "Utils"))

def main():
    parser = argparse.ArgumentParser(
        description="Parkinson's Disease GWAS Cell Type-Specific Enhancer Enrichment Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  # 1. Coordinate conversion
  python main.py --step coordinate

  # 2. LDSC analysis
  python main.py --step ldsc

  # 3. Visualization
  python main.py --step visualize

  # Run full pipeline
  python main.py --all
        """
    )

    parser.add_argument(
        '--step',
        choices=['coordinate', 'ldsc', 'visualize'],
        help='Select step to run'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run full pipeline'
    )

    args = parser.parse_args()

    if args.all:
        print("=" * 60)
        print("Parkinson's Disease GWAS Cell Type-Specific Enhancer Enrichment Analysis")
        print("Running full pipeline")
        print("=" * 60)

        # Step 1: Coordinate conversion
        print("\n[1/3] Coordinate conversion...")
        try:
            from setup_liftover import convert_all_enhancer_files
            convert_all_enhancer_files()
            print("✅ Coordinate conversion completed")
        except Exception as e:
            print(f"⚠️  Coordinate conversion skipped (already completed or error occurred): {e}")

        # Step 2: LDSC analysis (optional - may fail if LDSC not configured)
        print("\n[2/3] LDSC analysis...")
        try:
            from ldsc_analysis_system import main as ldsc_main
            ldsc_main()
            print("✅ LDSC analysis completed")
        except Exception as e:
            print(f"⚠️  LDSC analysis skipped (LDSC setup required or error): {e}")
            print("   (Visualization will continue)")

        # Step 3: Visualization
        print("\n[3/3] Visualization...")
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent / "1.Scripts" / "Visualization"))
            import celltype_manhattan_plot
            # Run main code from the module
            celltype_manhattan_plot.annotations = celltype_manhattan_plot.load_celltype_annotations()
            celltype_manhattan_plot.gwas_data = celltype_manhattan_plot.load_gwas_with_positions()
            celltype_manhattan_plot.annotations = celltype_manhattan_plot.compute_snp_overlaps(
                celltype_manhattan_plot.gwas_data, celltype_manhattan_plot.annotations
            )
            celltype_manhattan_plot.create_celltype_manhattan_plots(
                celltype_manhattan_plot.gwas_data, celltype_manhattan_plot.annotations
            )
            celltype_manhattan_plot.create_comparison_manhattan(
                celltype_manhattan_plot.gwas_data, celltype_manhattan_plot.annotations
            )
            print("✅ Visualization completed")
        except Exception as e:
            print(f"❌ Visualization failed: {e}")
            import traceback
            traceback.print_exc()

        print("\n✅ Full pipeline completed!")

    elif args.step == 'coordinate':
        print("Running coordinate conversion...")
        try:
            from setup_liftover import convert_all_enhancer_files
            convert_all_enhancer_files()
            print("✅ Coordinate conversion completed")
        except Exception as e:
            print(f"❌ Coordinate conversion failed: {e}")
            import traceback
            traceback.print_exc()

    elif args.step == 'ldsc':
        print("Running LDSC analysis...")
        try:
            from ldsc_analysis_system import main as ldsc_main
            ldsc_main()
            print("✅ LDSC analysis completed")
        except Exception as e:
            print(f"❌ LDSC analysis failed: {e}")
            print("Please check if LDSC environment is configured.")
            import traceback
            traceback.print_exc()

    elif args.step == 'visualize':
        print("Running visualization...")
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent / "1.Scripts" / "Visualization"))
            import celltype_manhattan_plot
            # Run main code
            annotations = celltype_manhattan_plot.load_celltype_annotations()
            gwas_data = celltype_manhattan_plot.load_gwas_with_positions()
            annotations = celltype_manhattan_plot.compute_snp_overlaps(gwas_data, annotations)
            celltype_manhattan_plot.create_celltype_manhattan_plots(gwas_data, annotations)
            celltype_manhattan_plot.create_comparison_manhattan(gwas_data, annotations)
            print("✅ Visualization completed")
        except Exception as e:
            print(f"❌ Visualization failed: {e}")
            import traceback
            traceback.print_exc()

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
