# ENY2980
Course-based Undergraduate Research Experience in bioinformatics and differential gene expression in sphinx moths (Lepidoptera: Sphingidae)   🌙🦋🧬😁

<p align="center">
<img width="450px" src="./Images/Hyles_lineata_on_sheet_220826.png">
</p>

## Our Data
+ You have access to gene expression data from the **mouthparts**, **claspers/ovipositor**, and **legs** from **males** and **females**.
+ In your group, you will come up with a hypothesis and a general experimental design using this data.
+ Data for this project is in ```/blue/eny2890/share/RNASeq/```
+ There are subfolders that contain the data for each sex x body part combination.


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
4. Map reads to indexed genome & count
5. Run differential gene expression analysis



<p align="center">
<img width="300px" src="./Images/DGE_Pipeline.png">
</p>



### 1. Run Quality check on reads using fastqc

#### Resources:
+ https://github.com/s-andrews/FastQC

#### (A) Navigate to the folder where your fastq.gz files are located, then adapt the fastqc.sh script to run quality check on the reads you have chosen to use in your analysis

```cd /blue/eny2890/your_name/your_data_folder```

#### (B) Now, create a submission or bash script (.sh) to submit your fastqc job to the HiPerGator scheduler.

+ the submission script must start with #!/bin/sh
+ the script filename extension should be .sh

```nano fastqc.sh``` to create a new, blank file, then write the following (adapted for your data) in the file


```bash
#!/bin/sh
#SBATCH --job-name=fastqc_rna	# Job name that will show up in queue
#SBATCH --output=fastqc_%j.out	# Running + error log output
#SBATCH --mail-type=END,FAIL	 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=your.email@ufl.edu>	# Your email address for notifications
#SBATCH --nodes=1		# Use one node
#SBATCH --ntasks=1		# Run a single task
#SBATCH --cpus-per-task=4	# number of cpus
#SBATCH --mem=4gb		# memory
#SBATCH -t 10:00:00 		# time limit in d-hh:mm:ss


module load fastqc
fastqc FileName_R1_001.fastq.gz FileName_R2_001.fastq.gz 

```

#### (C) To submit a slurm script to the scheduler, use the ```sbatch``` command:
```sbatch fastqc.sh```

#### (D) You can check if your job is running
```squeue -u your.name```



### 2. Trimming

If there were adapters in your reads, you should run a trimming step to remove them, then re-run fastqc to make sure the adapters were removed properly. 

To do this you will use a program called [**Trimmomatic**](http://www.usadellab.org/cms/index.php?page=trimmomatic).

In the 02_Trim folder at the top of this page, you will find an example submission script and the adapter sequence file you need

#### Resources
+ https://github.com/timflutre/trimmomatic/tree/master
+ Trimming and filtering lesson on carpentries lab github: https://carpentries-lab.github.io/metagenomics-analysis/03-trimming-filtering/index.html




### 3. part a: Indexing the genome 

Once you are satisfied with the quality of your reads, whether you trimmed them or not, put all of the fastq.gz files (one read pair with *R1 and *R2 for each sample) into one folder.   

To make sure you have all the read files you need for the differential gene expression anaylsis, you can use the following command to count the number of files in the folder: ```ls -1 | wc -l``` which says "list the files in here one per line and pipe list to a word count of the list, which gives you the total number of files in the folder.  

You can now ```rm``` any untrimmed fastq.gz files you do not need anymore.  

We only need one version of the indexed genome that everyone can call upon in their mapping script. We indexed this in class on 10/11/2023 and the submission script to do this is located in the 03_Index_Genome folder at the top of this page. This script created a folder that contains the full genome in fasta format and a series of .ht2 files with indexing information.  

The file path on HiPerGator is: ```/blue/eny2890/share/03_index_genome/index```



### 3. part b: Mapping reads 

Now that we have our RNA reads and our indexed genome, we can figure out where in the genome the RNA reads are located. To do this we use a program called [**Hisat2**](https://daehwankimlab.github.io/hisat2/).   

This is actually the same program we used to index the genome, but here we will use a more complicated submission script to get *both sets of reads mapped to the genome at once* and output a single mapped reads file for both of them.   


<p align="center">
<img width="400px" src="./Images/231014_ReadMappingwHisat2_sh01.png">
</p>  




<p align="center">
<img width="400px" src="./Images/231014_ReadMappingwHisat2_sh02.png">
</p>



