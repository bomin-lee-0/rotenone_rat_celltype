#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --time=00:20:00
#SBATCH --job-name=fastqc_batch
#SBATCH --array=1-80

module load fastqc

cd /scratch/prj/bcn_pd_pesticides/ratlas/Mingshen

sample=$(head -n $SLURM_ARRAY_TASK_ID filtered_samplename.txt | tail -n 1)

INPUT_FILE="/scratch/prj/bcn_pd_pesticides/ratlas/raw_reads/${sample}.fastq.gz"

fastqc -d . -o /scratch/prj/bcn_pd_pesticides/ratlas/Mingshen/fastqc_multiple "$INPUT_FILE"
