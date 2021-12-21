#!/bin/bash -l
#SBATCH -J compile # name that will be displayed in queue system
#SBATCH --output="compile.out" # stdout will be saved here
#SBATCH -N 1 # here set number of nodes 
#SBATCH -n 1 #  
#SBATCH --mem-per-cpu=1GB  
#SBATCH --time=00:05:00 # here set walltime hh:mm:ss
#SBATCH -p standard

## ------ QUEUE SYSTEM ------
## For submission use:
##      sbatch job.sh
## For checking status use:
##      squeue
## To kill job use:
##      scancel pid

cd $SLURM_SUBMIT_DIR 
source ./env.sh

# compile your code
make
