## A useful reference for working on HiperGator

### Navigating into your folder no matter where you are
```cd /blue/eny2890/your.folder.name```

### Listing the contents of a director
```ls directory name```

### Checking what jobs you are running
```squeue -u your.name```

### Cancelling a job 
To cancel a job, run the squeue command and look copy-paste the JOBID
```scancel JOBID```

### Moving a file 
To move a file from one directory to another use
``` mv where/the/file/is/filename.ext where/it/should/be/moved/```



### Starting a development node for running a small job or testing a sample analysis

To request 4gb of memory and 4 cpus for 1 hour, use

```module load ufrc 
srundev --mem=4gb --ntasks=1 --cpus-per-task=4 --time=01:00:00 
```