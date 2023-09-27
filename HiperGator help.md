## A useful reference for working on HiperGator

### Checking what jobs you are running
```squeue -u your.name```

### Cancelling a job 
To cancel a job, run the squeue command and look copy-paste the JOBID
```scancel JOBID```



### Starting a development node for running a small job or testing a sample analysis

To request 4gb of memory and 4 cpus for 1 hour, use

```module load ufrc 
srundev --mem=4gb --ntasks=1 --cpus-per-task=4 --time=01:00:00 
```