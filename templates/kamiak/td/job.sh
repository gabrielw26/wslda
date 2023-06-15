#!/bin/bash -l
#SBATCH --partition=kamiak         # SET: Partition/Queue to use
#SBATCH --job-name=myJob           # SET: Job name
#SBATCH --output=myJob_%j.out      # SET: Output file (stdout)
#SBATCH --error=myJob_%j.err       # SET: Error file (stderr)
#SBATCH --mail-type=ALL            # SET: Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=your.name@wsu.edu  # SET: Email address for notifications
#SBATCH --time=01:00:00             # SET: Wall clock time limit Days-HH:MM:SS
#SBATCH --nodes=2                  # SET: Number of nodes (min-max)
#SBATCH --ntasks-per-node=4        # SET: Number of tasks per node (max)
#SBATCH --gres=gpu:tesla:4           # SET: number of gpus pre node
#SBATCH --cpus-per-task=1         # DO NOT MODIF: Number of cores per task (threads)

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

# to check GPU's
nvidia-smi 

# Execute with nodes*ntasks-per-node processes
# srun ./td-wslda-1d input.txt
# srun ./td-wslda-2d input.txt
srun ./td-wslda-3d input.txt

