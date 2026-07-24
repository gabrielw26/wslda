#!/bin/sh
#$ -cwd
#$ -l node_f=4      # <- SET: number of node
#$ -l h_rt=01:00:00 # <- SET: walltime hh:mm:ss
#$ -N td-wslda      # <- SET: name of your job
#$ -m abe           # activate e-mail notification
#$ -M your.mail@pw.edu.pl # <- SET: your e-mail
#$ -p -5            # <- SET: priority, -5 (standard), -4, -3 (highest) correspond to the priority 0, 1, 2 of the charging rule.

# Load modules
source ./env.sh

#           <-- SET to node_f*4 (each node has 4 GPUs), do NOT modify npernode=4 value 
mpirun -n 16 -npernode 4 -x LD_LIBRARY_PATH ./td-wslda-3d input.txt


## SUBMISSION:
## qsub -g hp190063
