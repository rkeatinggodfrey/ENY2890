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
```scancel {JOBID{}```

### Moving a file 
To move a file from one directory to another 
``` mv {path/where/the/file/is/filename.ext} {path/where/it/should/be/moved/}```

### Moving a file from some other place to your current location
```mv {path/where/the/file/is/filename.ext} .```

### Copying a file from some other place to your current location
```cp {path/where/the/file/is/filename.ext} .```

### Copying a directory from some other place to your current location
```cp -r {path/to/the/directory} .```


### Starting a development node for running a small job or testing a sample analysis

To request 4gb of memory and 4 cpus for 1 hour, use

```module load ufrc 
srundev --mem=4gb --ntasks=1 --cpus-per-task=4 --time=01:00:00 
```