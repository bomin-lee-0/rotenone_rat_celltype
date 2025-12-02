#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/build_index.out

cd /scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/bowtie2_index

module load bowtie2

bowtie2-build GCF_015227675.2_mRatBN7.2_genomic.fna mRatBN7.2
