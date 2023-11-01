## A useful reference for working on HiperGator

### A silly set of gifs to help you remember

[Linux commands](https://docs.google.com/presentation/d/1rGAfFQNdW3bNnVxPWDT8Dgj-MlUr970G-IZ-5kykWHU/edit#slide=id.p)


### Navigating into your folder no matter where you are
```cd /blue/eny2890/{your.folder.name}```

### Listing the contents of a director
```ls {directory name}```  
```ls -l {directory name}```

### Checking what jobs you are running
```squeue -u {your.name}```

### Cancelling a job 
To cancel a job, run the squeue command and look copy-paste the JOBID  
```scancel {JOBID}```

### Removing a file
```rm filename```

### Removing a directory and all of its contents
```rm -r /directory ```

### Moving a file 
To move a file from one directory to another  
``` mv {path/where/the/file/is/filename.ext} {path/where/it/should/be/moved/}```

### Moving a file from some other place to your current location
```mv {path/where/the/file/is/filename.ext} .```

### Copying a file from some other place to your current location
```cp {path/where/the/file/is/filename.ext} .```

### Copying a directory from some other place to your current location
```cp -r {path/to/the/directory} .```

### Copying all of the fastq.gz files from a directory to your current location 

```cp -path/to/the/directory/*.fastq.gz .```

### Copying all files with "cleaned" in the name from a directory to your current location 

```cp -path/to/the/directory/*cleaned* .```
Here the asterisks indicate any file that has cleaned in the name flanked by any other characters

### Starting a development node for running a small job or testing a sample analysis

To request 4gb of memory and 4 cpus for 1 hour, use

```module load ufrc 
srundev --mem=4gb --ntasks=1 --cpus-per-task=4 --time=01:00:00 
```

### Working with FASTA files

To retain the entire header  

```grep ">" {name_of_file}.fasta | sed 's,>,,g' > ors_prot_ms_acc.txt```

To retain only the >IDs  

```grep "^>" {name_of_file}.fasta | cut -d' ' -f1 > test```