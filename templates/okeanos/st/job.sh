#!/bin/bash -l
#SBATCH -J test # name that will be displayed in queue system
#SBATCH --output="test.out" # stdout will be saved here
#SBATCH -N 2 # here set number of nodes
#SBATCH -n 48 # must be N*24
#SBATCH -A g87-1072
#SBATCH --mem=124000
#SBATCH --hint=nomultithread
#SBATCH --time=00:30:00 # here set walltime hh:mm:ss
#SBATCH --mail-type=BEGIN,END,FAIL # notifications for job started, done & fail
#SBATCH --mail-user=your@email # send-to address

## -------- WIKI INFO -------
## https://gitlab.fizyka.pw.edu.pl/wtools/w-bsk/-/wikis/Setting%20up%20calculations
## https://gitlab.fizyka.pw.edu.pl/wtools/w-bsk/-/tree/devel/templates/okeanos
##
## ------ QUEUE SYSTEM ------
## For submission use:
##      sbatch job.sh
## For checking status use:
##      squeue
## To kill job use:
##      scancel pid


cd $SLURM_SUBMIT_DIR

setenv LD_LIBRARY_PATH /opt/cray/diag/lib
module load cray-fftw

## NOTE: YOU MUST ASSURE: -n=-N*24
srun -n 48 ./st-wbsk-2d input.txt
