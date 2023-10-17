## Mapping RNAreads to a reference genome using Hisat2

In order to know what genes our RNA reads (sequences) came from, we need to map them to a reference genome. Luckily, we have a high quality reference genome for *Hyles lineata* [Godfrey et al. 2022](https://academic.oup.com/g3journal/article/13/6/jkad090/7147210)!  

We can use the program [Hisat2](https://daehwankimlab.github.io/hisat2/) to do this.  

## Hisat2 instructions  

This is actually the same program we used to index the genome, but here we will use a more complicated submission script to get *both sets of reads mapped to the genome at once* and output a single mapped reads file for both of them.   

You should create the directory where you want the mapped reads deposited. Something like:  
```mkdir mapped```

Let's break down the script so you can modify and use it to run the mapping step:

### (1) Submission script heading   

The script needs to start with #! (she-bang) and needs to be a ".sh" script. **CHANGE** the email address to your own.  

The .out file will be the log of the program running. You can call it a .out file or a .log file! The ```%j``` part of the name means "job number".  

```
#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --output=hisat2_%j.out  # running & error log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your.email@ufl.edu> # CHANGE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 72:00:00             # time limit in d-hh:mm:ss

# Load required modules         # load the programs you need
module load hisat2
module load samtools
```

### (2) Defining directories  

In order to map to the genome, we need to tell Hisat2 where the indexed genome is located, where our RNA reads are located, and where you want the mapped read files (in .bam and .bai format) to be deposited.  `  

Let's look at the script to see what directories we need to define.  

This line defines the genome index. It needs to include the prefix of the indexed file, 'Hlineata'  
```genome_index="/blue/eny2890/share/03_index_genome/index/Hlineata" ```  

**CHANGE THIS** to the directory where your reads are located  
```reads_dir="/blue/eny2890/{path_to_your_curated_reads}"```

**CHANGE THIS** is the directory you made for the mapped .bam and .bai files  
```output_dir="/blue/eny2890/{path_to_mapped_reads_directy_you_created}"``` 


### (3) Defining R1 file names

this line navigates into the reads directory you defined above as "reads_dir"  
```cd "$reads_dir"```  

This line lists all R1 files in the directory. This depends on the pattern of you file names.  
```read_files_R1=$(ls *R1_clean.fastq.gz)```  

For example, if your file names end with R1_001.fastq.gz, you can change the ```*R1_clean.fastq.gz``` to ```*R1_001.fastq.gz```  

### (4) Loop through R1 and R2 files 

Let's look a what the loop says in plain English.  

For the read files in this folder, do the following  
```for r1 in $read_files_R1; do```  

Extract sample names. **if** your names end with R1_001.fastq.gz, you can **CHANGE** this to "_R1_001.fastq.gz"  
```sample_name=$(basename "$r1" _R1_clean.fastq.gz)```  

Define R2 based on R1 file name so you do not have to enter it again. **If** your file name ends with _001.fastq.gz you can **CHANGE** ```_clean``` to ```_001```  
```r2=$(echo "$r1" | sed 's/_R1_clean/_R2_clean/')```

### (5) Perform HISAT2 mapping  

FINALLY we run Hisat2! You don't need to change anything here. This calls on the genome index and the r1 and r2 file names we already defined above.  
```hisat2 -p 8 -x "$genome_index" -1 "$r1" -2 "$r2" | samtools sort -@ 8 -o "$output_dir/$sample_name.sorted.bam"```

We also need a mapped reads index for counting genes, so we create that here. It will be in .bai format.    
```samtools index $output_dir/$sample_name.sorted.bam```

Remove any large sam files created in this process  
```rm $output_dir/$sample_name.sam```

And we are  
```done```
