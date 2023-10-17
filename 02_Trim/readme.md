#### Step by step trimmomatic loop script instructions 

(1) Make a folder for trimming

```mkdir trim```

(2) Navigate to the folder where your .fastq.gz RNA read files are stored and move all of the .fastq.gz files to the trim folder

```mv *fastq.gz trim/```

(3) Navigate into that folder to work from there

```cd trim/``` 

(4) Copy the TruSeq3-PE.fa into your trim folder (where you should currently be located)

```cp /blue/eny2890/share/02_trimmomatic/TruSeq3-PE.fa .```

(you can look at it to make sure it matches what is on the git)
```cat TruSeq3-PE.fa``` to view file

(5) Copy the trimmomatic_loop.sh into your trim folder

```cp /blue/eny2890/share/02_trimmomatic/trimmomatic_loop.sh .```

nano the file to edit it and change the email address to your own. 

```nano trimmomatic_loop.sh```

As long as the *.fastq.gz file have their original names (e.g., Hl22042008-f-G_S5_L001_R1_001.fastq.gz) the script will be able to run on them as it is written.

(6) Run the loop script

```sbatch trimmomatic.sh```


(7) Make sure your job is running

```-squeue -u {your.name}```

(8) Look at the log file to see it running or any errors

```cat trimmo_{your_job_number}.log```

