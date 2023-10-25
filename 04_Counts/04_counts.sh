#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --output=counts_%j.out  # running + error log output
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={your.name}@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 4-12:00:00           # time limit in d-hh:mm:ss

# Load required modules
module load htseq

# Set paths and variables
input_dir="/blue/eny2890/{path_to_your_reads_folder}" # add your reads file directory
gtf_file="/blue/eny2890/share/04_counts/augustus.hints.gtf.gff" # this is the annotated genome; you don't have to change this one!
output_dir="/blue/eny2890/{path_to_your_count_files_folder}" # use the mkdir command to make a folder for counts and list it here

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Get the list of BAM files
bam_files=$(ls $input_dir/*.bam)

# Iterate over each BAM file
for bam_file in $bam_files; do
    # Extract the sample name from the BAM file name
    sample_name=$(basename $bam_file .bam)
    # Run HTSeq to calculate read counts
    htseq-count -f bam -s no -i Parent -t exon $bam_file $gtf_file > $output_dir/$sample_name.counts.tsv
done
