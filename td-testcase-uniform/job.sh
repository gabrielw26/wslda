#!/bin/bash
##        <-- SET name 
#PBS -N td-wslda
##        <-- SET number of nodes, ppn=40 to lock the whole node or ppn=8 to allow sharing with others
#PBS -l nodes=2:ppn=40
##        <-- SET memory request, increase it if you encounter problem of type: `PBS: job killed: mem job total exceeded limit`
#PBS -l mem=128gb
##        <-- SET walltime, correlate it with `walltime` input file variable
#PBS -l walltime=12:00:00  
#PBS -j oe
#PBS -q long

## ------ QUEUE SYSTEM ------
## For submission use:
##      qsub job.sh
## For checking status use:
##      qstat
## To kill job use:
##      qdel jobid
##       
## ------ COMPUTATION -------     
## For computation you must use /home2/scratch folder
##      cd /home2/scratch
## NOTE: files older than 100days are automatically removed from  /home2/scratch
## For storing results use location:
##      cd /home2/archive
##      
## NOTE: dwarf is not heterogeneous (different nodes have different number of GPUs of a different type),
## method of distributing tasks across nodes is provided in machine.h. 

## execute code
cd $PBS_O_WORKDIR
source ./env.sh

##         <--- SET n = 8*nodes, do NOT modify `-npernode 8`
mpirun -n 16 -npernode 8 ./td-wslda-2d input.txt


