#!/bin/bash -l
#SBATCH -J test # name that will be displayed in queue system
#SBATCH --output="test.out" # stdout will be saved here
#SBATCH -N 1 # here set number of nodes 
#SBATCH -n 32 # must be N*32 
#SBATCH --mem-per-cpu=2GB  
#SBATCH --gpus-per-node=8 # use 8 if you want to use Altair nodes
#SBATCH --time=24:00:00 # here set walltime hh:mm:ss
#SBATCH --mail-type=BEGIN,END,FAIL # notifications for job started, done & fail
#SBATCH --mail-user=your@email # send-to address
#SBATCH -p tesla 

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
## cd ~/grant_518/scratch/

cd $SLURM_SUBMIT_DIR 
source ./env.sh

## NOTE: YOU MUST ASSURE: -np=N*gpus-per-node
##                        -npernode=gpus-per-node
mpirun -np 8 -npernode 8 ./td-wslda-2d input.txt 


