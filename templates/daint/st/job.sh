#!/bin/bash -l
#SBATCH --job-name=NAME
#SBATCH --output="NAME.out"
#SBATCH --error="NAME.err"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriel.wlazlowski@pw.edu.pl
#SBATCH --time=24:00:00
#SBATCH --nodes=288
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH -A pr125

## -------- WIKI INFO -------
## https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Piz%20Daint
## 
## ------ QUEUE SYSTEM ------
## For submission use:
##      sbatch job.sh
## For checking status use:
##      squeue -u gabrielw
## To kill job use:
##      scancel pid
##       
## ------ COMPUTATION -------     
## For computation you must use $SCRATCH folder
##      cd $SCRATCH
## NOTE: files older than 30days are automatically removed from  $SCRATCH
## For storing results use location:
##      cd /project/pr125/

# Load modules
source env.sh

# execute (with GPU-ELPA)
## YOU MUST ASSURE: -n=--nodes 
srun -n 288 --ntasks-per-node=1 ./st-wslda-3d input.txt

