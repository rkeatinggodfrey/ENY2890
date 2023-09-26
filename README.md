# ENY2980
Course in Undergraduate Research Education in bioinformatics and differential gene expression 😁

In this course you will work in small teams to compare gene expression in the whiteline sphinx moth (*Hyles lineata*)


#### Pipeline overview

1. Run Quality check on reads
2. Trim adapters 
3. Index the genome 
4. Map reads to indexed genome
5. Count reads mapped to genome



### 1. Run Quality check on reads using fastqc

#### Resources:
+ https://github.com/s-andrews/FastQC

Navigate to the folder where your fastq.gz files are located, then adapt the fastqc.sh script to run quality check on the reads
you have chosen to use in your analysis

```bash
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
fastqc Hl22042008-f-P_S10_L001_R1_001.fastq.gz Hl22042008-f-P_S10_L001_R2_001.fastq.gz

```



