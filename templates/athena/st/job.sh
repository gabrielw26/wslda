#!/bin/bash -l
#SBATCH -J test # <-- SET: name that will be displayed in queue system
#SBATCH --output="test.%J.out" # <-- SET: stdout will be saved here
#SBATCH -N 2 # <-- SET: number of nodes 
#SBATCH -n 256 # <-- SET: must be N*128
#SBATCH --gres=gpu:8 # DO NOT MODIFY
#SBATCH --time=02:00:00 # <-- SET: walltime hh:mm:ss
#SBATCH --mail-type=BEGIN,END,FAIL # notifications for job started, done & fail
#SBATCH --mail-user=your@email # <-- SET: send-to address
#SBATCH -p plgrid-gpu-a100
#SBATCH -A plginhsf-gpu-a100

## ------ QUEUE SYSTEM ------
## For submission use:
##      sbatch job.sh
## For checking status use:
##      squeue
## To kill job use:
##      scancel pid
## 
## ------ COMPUTATION -------     
## For computation it is recomended to use location:
## cd $SCRATCH

cd $SLURM_SUBMIT_DIR 
source ./env.sh

## NOTE: YOU MUST ASSURE: np=N*npernode
## RECOMMENDED: 64 processes per node, ie. 8 processes per GPU
mpirun -np 128 -npernode 64 ./st-wslda-2d input.txt 

