#!/bin/bash
#SBATCH --job-name=trim_{identifier}
#SBATCH --output=trim_%j.log
#SBATCH --mail-user=your.email@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

module load trimmomatic/0.39


trimmomatic PE -threads 4 \
{Forward_Read_R1}.fastq.gz {Reverse_Read_R2}.fastq.gz \
trimmed_{Forward_Read_R1}.fastq.gz trimmed_unparied_{Forward_Read_R1}.fastq.gz \
trimmed_{Reverse_Read_R2}.fastq.gz trimmed_unpaired_{Reverse_Read_R2}.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36