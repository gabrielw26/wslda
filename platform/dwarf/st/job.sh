#!/bin/bash
#PBS -N st-wslda
#PBS -l nodes=1:ppn=40
## If you want to submit to specific node use this: 
##PBS -l nodes=node2067.grid4cern.if.pw.edu.pl:ppn=40
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -q long

## -------- WIKI INFO -------
## http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/dwarf
## 
## ------ QUEUE SYSTEM ------
## to submit use:
##     qsub job.st.dwarf
## to check status use:
##     qstat 
## to delete job use:
##     qdel jobid
##       
## ------ COMPUTATION -------    
## For computation with ELPA single node runs are the most efficient
## Recommended node is 2067.    
## For computation use /home2/scratch
##      cd /home2/scratch
## NOTE: files older than 100days are automatically removed from  /home2/scratch
## For storing results use location:
##      cd /home2/archive

# execute code
cd $PBS_O_WORKDIR
source env.sh

mpirun -n 20 ./st-wslda-2d input.txt

