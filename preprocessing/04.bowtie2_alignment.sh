#!/bin/bash

#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --output=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/bowtie2_q_newref.out
#SBATCH --array=1-34
#SBATCH --time=06:00:00
#SBATCH --job-name=bowtie2_align

module load bowtie2


DATA_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/trimmed_exclude_output
SAM_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/SAM
GENOME=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/bowtie2_index/mRatBN7.2
FILELIST=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/filtered_sample_withoutR.txt

file=$(head -n $SLURM_ARRAY_TASK_ID $FILELIST | tail -1)

echo "****** Running alignment for: $file ******"

bowtie2 -q -p 8 -x $GENOME \
  -1 ${DATA_DIR}/${file}_R1_001_val_1.fq.gz \
  -2 ${DATA_DIR}/${file}_R2_001_val_2.fq.gz \
  -S ${SAM_DIR}/${file}.sam \
  --very-sensitive
