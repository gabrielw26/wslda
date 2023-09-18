#!/bin/bash

#SBATCH --job-name=NAME          # <--- SET 
#SBATCH --output="NAME.%J.out"   # <--- SET
#SBATCH --error="NAME.%J.err"    # <--- SET
#SBATCH --nodes=2                # <--- SET: Number of nodes, each noode has 8 GPUs
#SBATCH --ntasks=16              # <--- SET: Number of processes you want to use, MUST be nodes*8 !!!
#SBATCH --gpus=16                # <--- SET: MUST be the same as ntasks !!!
#SBATCH --time=02:15:00             # <--- SET: Walltime HH:MM:SS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your@email   # <--- SET: if you want to get e-mail notification
#SBATCH --partition=standard-g
#SBATCH --account=project_465000307 
#SBATCH --cpus-per-task=1        # Do not modify
#SBATCH --ntasks-per-node=8      # Do not modify
#SBATCH --gpus-per-node=8        # Do not modify
#SBATCH --exclude=nid005917

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
##      cd /scratch/project_465000307
## or project FAST SCRATCH
##      cd /flash/project_465000307
## For more info see: https://docs.lumi-supercomputer.eu/storage/
## 
## For storing results use location:
##      cd /project/project_465000307/

# Set environment
# export MPICH_GPU_SUPPORT_ENABLED=1 
source ./env.sh

# Execute the code
#        <--- NOTE: MUST be the same as ntasks !!!
srun -n 16 ./st-wslda-3d input.txt



