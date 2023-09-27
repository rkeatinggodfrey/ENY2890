#!/bin/sh
#SBATCH --job-name=fastqc_rna	# Job name that will shop up in queue
#SBATCH --output=fastqc_%j.out	# Running + error log output
#SBATCH --mail-type=END,FAIL	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=your.email@ufl.edu>	# Your email address for notifications
#SBATCH --nodes=1		# Use one node
#SBATCH --ntasks=1		# Run a single task
#SBATCH --cpus-per-task=4	# number of cpus
#SBATCH --mem=4gb		# memory
#SBTACH -t 10:00:00 		# time limit in d-hh:mm:ss


module load fastqc
fastqc file.R1.fq file.R2.fq    #You can specify as many files to process in a single run as you like. 