#!/bin/bash
#SBATCH --job-name=trimmo
#SBATCH --output=%x_%j.log
#SBATCH --mail-user={your.name}@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load trimmomatic/0.39

for sample in $(ls *fastq.gz | cut -d "_" -f 1,2,3 | sort | uniq)
do
    fq1=$(ls ${sample}_R1*)
    fq2=$(ls ${sample}_R2*)
    trimmomatic PE -threads 4 \
    ${fq1} ${fq2} \
    ${sample}_R1_clean.fastq.gz ${sample}_R1_unpaired.fastq.gz \
    ${sample}_R2_clean.fastq.gz ${sample}_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:36
done
