#!/bin/bash

#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --output=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/dedup.out
#SBATCH --array=1-17
#SBATCH --time=24:00:00
#SBATCH --job-name=dedup

module load picard

DATA_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen
BAM_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/BAM
TMP_DIR=/scratch/prj/bcn_marzi_lab/ratlas/Mingshen/alignment/BAM/tmp

cd $DATA_DIR

file=$(head -n $SLURM_ARRAY_TASK_ID sample_17.txt | tail -1)

picard MarkDuplicates \
  --INPUT ${BAM_DIR}/sorted/$file.sorted.bam \
  --OUTPUT ${BAM_DIR}/sorted/$file.dedup.bam \
  --TMP_DIR $TMP_DIR \
  --METRICS_FILE ${BAM_DIR}/sorted/$file.dedup.metrics \
  --REMOVE_DUPLICATES true \
  --VALIDATION_STRINGENCY LENIENT

