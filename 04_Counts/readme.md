### Gene counts with HTSeq

In this step of the pipeline we will use our aligned sequences reads (.bam files) and a genome file (the .gtf file in the counts folder) that includes genomic features (e.g., exons) to count how many reads map to each feature.  

We will use a module on HiPerGator called ```htseq-count``` to do this. You can find information about the settings for this program (here)[https://htseq.readthedocs.io/en/release_0.11.1/count.html].  

### Submitting a job to run HTSeq using a loop script


#### Setting up the directories  

Set the filepath to your bam files
```input_dir="/blue/eny2890/{path_to_your_mapped_reads_folder}"```  

This is the file path to the genome feature file (.gtf) you can use, so you don't need to change this.  
```gtf_file="/blue/eny2890/share/04_counts/augustus.hints.gtf.gff"```   


Create a directory for your count files to be stored and set the filepath here.
```output_dir="/blue/eny2890/{path_to_your_count_files_folder}"```


#### Running htseq-count
In our .sh script, we run htseq-count using the following parameters:  

```htseq-count -f bam -s no -i Parent -t exon $bam_file $gtf_file > $output_dir/$sample_name.counts.tsv```  

+ This line of code tells htseq-count that we are providing .bam files (```-f bam```). 
+ The ```-t exon``` flag indicates we want to count reads that mapped to exons and not any of the other features in the third column of our .gtf file. 
+ The ```-i Parent``` flag in the code says we want reads that map to features with this same id as the same. In our genome gtf feature file, "Parent" indicates the primary gene ID.


### Interpreting your count data  

Once your count files are finished, you can use the ```tail``` command to look at some statistics on how many reads mapped to features ("Parent" gene IDs).  

Navigate to the counts folder and take a look at any file using 

```tail {filename}.sorted.counts.csv```  

You can read definitions for all of these classifications on the [htseq-count page](https://htseq.readthedocs.io/en/release_0.11.1/count.html#0)  

The classifications do not include how many reads actually did map. We can calculate that by adding up the second column of the counts file using this line of code:   
```awk '{s+=$2}END{print s}' Hl22042001-f-G_S1_L001.sorted.counts.csv```