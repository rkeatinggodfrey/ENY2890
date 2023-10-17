#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --output=hisat2_%j.out  # running + error log output
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 72:00:00             # time limit in d-hh:mm:ss

# Load required modules 
module load hisat2
module load samtools

# Set paths to variables
genome_index="/blue/eny2890/share/03_index_genome/index" #this needs to be the path to the name of indexed files in a folder with the genom
e fasta
reads_dir="/blue/eny2890/{path_to_your_curated_reads}" #file path to the location where all of your final reads are (trimmed or not)
output_dir="/blue/eny2890/{path_to_mapped_reads_directy_you_created}" #use the mkdir command to make a directory for the mapped files, .bam
 and .bai, to be deposited in

# Change to the reads directory
cd "$reads_dir"

# List all R1 files in the directory
read_files_R1=$(ls *R1_clean.fastq.gz)

# Loop through R1 files and perform HISAT2 mapping
for r1 in $read_files_R1; do
    # Extract sample name from R1 filename (assuming the filename format is something like: SampleName_R1_100.fastq.gz)
    sample_name=$(basename "$r1" _R1_clean.fastq.gz)
    # Define R2 filename based on R1 filename
    r2=$(echo "$r1" | sed 's/_R1_clean/_R2_clean/')
    # Run HISAT2
    hisat2 -p 8 -x "$genome_index" -1 "$r1" -2 "$r2" | samtools sort -@ 8 -o "$output_dir/$sample_name.sorted.bam"

# Create index file for each bam
samtools index $output_dir/$sample_name.sorted.bam

# Remove the intermediate SAM file
rm $output_dir/$sample_name.sam

done