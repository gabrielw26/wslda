#!/bin/sh
#$ -cwd
#$ -l node_f=4      # <- SET: number of nodes
#$ -l h_rt=01:00:00 # <- SET: walltime hh:mm:ss
#$ -N st-wslda      # <- SET: name of your job
#$ -m abe           # activate e-mail notification
#$ -M your.mail@pw.edu.pl # <- SET: your e-mail
#$ -p -5           # <- SET: priority, -5 (standard), -4, -3 (highest) correspond to the priority 0, 1, 2 of the charging rule.

# Load modules
source ./env.sh

## NOTE:
## When running and an hcoll-related error or segmentation fault occur, uncomment option below 
# export I_MPI_COLL_EXTERNAL=0

#                 <-- SET to node_f*4 (each node has 4 GPUs), do NOT change ppn=4 value
mpiexec.hydra -n 16 -ppn 4 ./st-wslda-3d input.txt 


## SUBMISSION:
## qsub -g hp190063
