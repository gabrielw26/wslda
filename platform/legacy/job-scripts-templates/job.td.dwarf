#!/bin/bash
#PBS -N td-wslda
#PBS -l nodes=7:ppn=20
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -q long

## -------- WIKI INFO -------
## http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/dwarf
## 
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
## NOTE: dwarf is not heterogeneous (different nodes have different number of GPUs of different type),
## you need to deliver proper information about structure of the system in dwarfnodes.txt file

module load cuda/9.0

mpirun -n 40 -hostfile dwarfnodes.txt ./td-wslda-2d input.txt


