#!/bin/bash
#SBATCH -A NFI125
#SBATCH -J NAME
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 1

# Set environment
source ./env.sh

# Srun command.
#             <---N*8       <---N*8
time srun -n 8 -c 1 --gpus 8 --gpus-per-node 8 --gpu-bind=closest ./td-wslda-3d input.txt

