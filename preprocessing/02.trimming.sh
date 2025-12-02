#!/bin/bash

#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/ratlas/Mingshen/trimmed_output/trim.%a.out
#SBATCH --array=1-40
#SBATCH --time=01:00:00

module load trimgalore/0.6.6-gcc-13.2.0-python-3.11.6

cd /scratch/prj/bcn_pd_pesticides/ratlas/Mingshen

sample=$(head -n $SLURM_ARRAY_TASK_ID sample_R1_only.txt | tail -n 1)

INPUT_R1="/scratch/prj/bcn_pd_pesticides/ratlas/raw_reads/${sample}.fastq.gz"
INPUT_R2="/scratch/prj/bcn_pd_pesticides/ratlas/raw_reads/${sample/_R1_/_R2_}.fastq.gz"

trim_galore --paired --fastqc \
  --fastqc_args "--outdir trimmed_output/fastqc" \
  --three_prime_clip_R1 5 \
  --three_prime_clip_R2 5 \
  -o trimmed_output \
  "$INPUT_R1" "$INPUT_R2"
