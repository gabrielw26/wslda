#!/bin/bash

#SBATCH --job-name=l48          # <--- SET 
#SBATCH --output="l48.%J.out"   # <--- SET
#SBATCH --error="l48.%J.err"    # <--- SET
#SBATCH --nodes=32                # <--- SET: Number of nodes, each noode has 8 GPUs
#SBATCH --ntasks=256              # <--- SET: Number of processes you want to use, MUST be nodes*8 !!!
#SBATCH --gpus=256                # <--- SET: MUST be the same as ntasks !!!
#SBATCH --time=02:00:00             # <--- SET: Walltime HH:MM:SS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriel.wlazlowski@pw.edu.pl   # <--- SET: if you want to get e-mail notification
#SBATCH --partition=pilot 
#SBATCH --account=project_465000150 
#SBATCH --cpus-per-task=1        # Do not modify
#SBATCH --ntasks-per-node=8      # Do not modify

## ------ QUEUE SYSTEM ------
## For submission use:
##      sbatch job.sh
## For checking status use:
##      squeue -u username
## To kill job use:
##      scancel pid
##       
## ------ COMPUTATION -------     
## For computation you must use SCRATCH folder
##      cd /scratch/project_465000150
## or project FAST SCRATCH
##      cd /flash/project_465000150
## For more info see: https://docs.lumi-supercomputer.eu/storage/
## 
## For storing results use location:
##      cd /project/project_465000150/

# Set environment
export MPICH_GPU_SUPPORT_ENABLED=1 
source ./env.sh
make 3d
cp ./td-wslda-3d ./td-wslda-3d-48

# Execute the code
#        <--- NOTE: MUST be the same as ntasks !!!
srun -n 256 ./td-wslda-3d-48 input.txt



