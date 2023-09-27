# ENY2980
Course-based Undergraduate Research Experience in bioinformatics and differential gene expression 🌙🦋🧬😁

<img width="30px" src="./ENY2890/images/Hyles_lineata_on_sheet_220826.jpg" alt="Hyles lineata Finger Rock 2022" />

## Our Data
+ You have access to gene expression data from the **mouthparts**, **genetalia**, and **legs** from **males** and **females**
+ In your group, you will come up with a hypothesis and a general experimental design using this data



## Broad Analysis Goals
In this course you will work in small teams to compare gene expression in the whiteline sphinx moth (*Hyles lineata*). We will focus on two analyses:
+ Differential Gene Expression (DGE) wherein we will count the number of times a gene shows up in each sample and compare across samples using statistics to determine if samples or groups of samples have show greater or fewer transcripts of a gene than others.
+ Presence/absence wherein we will simply look at what genes show up in each sample or group of samples and report the presence or absence of certain genes.




## Get your files!
+ Once you know what hypothesis you'd like to test, copy the fastq.gz files you need to your data folder (*call it what ya want, it's where you will do your analyes*)
+ Each pair of reads a sample are located in their own folder and can be distinguished by R1 (forward) or R2 (reverse) in their name.
+ To copy the sample folders you need, navigate to your data folder and use the cp command with the recursive flag -r so it copies the contents of the sample folder. 
```
cp -r /blue/eny2890/share/RNAseq/folder_of_reads_you_need .
```




## Pipeline overview
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

To submit a slurm script to the scheduler, use 
```sbatch scriptname.sh```
+ the submission script must start with #!/bin/sh
+ the script filename extension should be .sh



