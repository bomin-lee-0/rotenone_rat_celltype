#!/bin/bash

#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --output=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/samtools_transfer.out
#SBATCH --array=1-34
#SBATCH --time=04:00:00
#SBATCH --job-name=samtools_transfer

module load samtools

DATA_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen
BAM_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/BAM
SAM_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/SAM

cd $DATA_DIR

file=$(head -n $SLURM_ARRAY_TASK_ID filtered_sample_withoutR.txt | tail -1)

samtools view -b -o ${BAM_DIR}/${file}.bam ${SAM_DIR}/${file}.sam
